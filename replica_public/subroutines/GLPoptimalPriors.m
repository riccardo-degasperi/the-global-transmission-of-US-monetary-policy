function [out,csmin,postMode] = GLPoptimalPriors(Y,X,Ybar,Xbar,hp,dims)

    hyperpriors = hp.hyperpriors;

    % Retrieve dimensions:
    m = dims.m;
    [T0,N] = size(Y);
    [~,Npm] = size(X);
    p = (Npm - m)/N;
    T = T0 + p;

    % Retrieve list of priors to optimise
    priorlist = {'niw','lag','soc','cop','std'};    %possible options
    iPrior    = [hp.niw, hp.lag, hp.soc, hp.cop, hp.std] == 1;
    priors    = priorlist(iPrior);                  %selected options
    
    % Bounds for maximization
    MIN.l1 = 0.0001;        %lambda in GLP
    MIN.l2 = 0.1;           %alpha in GLP
    MIN.l3 = 0.0001;        %theta in GLP
    MIN.l4 = 0.0001;        %miu in GLP
    MIN.SS2 = hp.SS2./100;  %psi in GLP
    MAX.l1 = 5;
    MAX.l2 = 5;
    MAX.l3 = 50;
    MAX.l4 = 50;
    MAX.SS2 = hp.SS2.*100;

    
    % Initial values for maximisation of log-posterior
    l10 = -log((MAX.l1 - hp.l1)/(hp.l1 - MIN.l1));
    
    if hp.lag
    l20 = -log((MAX.l2 - hp.l2)/(hp.l2 - MIN.l2));
    else
    l20 = [];
    end
    
    if hp.soc
    l30 = -log((MAX.l3 - hp.l3)/(hp.l3 - MIN.l3));
    else
    l30 = [];
    end
    
    if hp.cop
    l40 = -log((MAX.l4 - hp.l4)/(hp.l4 - MIN.l4));
    else
    l40 = [];
    end
    
    if hp.std
    s20 = -log((MAX.SS2 - hp.SS2)./(hp.SS2 - MIN.SS2));
    else
    s20 = [];
    end
    
    % Initial value of the parameter vector 
    X0 = [l10; l20; l30; l40; s20];

    % Initial guess for the inverse Hessian
    H0 = 10*eye(length(X0));    %identity
    
    % Minnesota prior mean                      
    b              = zeros(N*p+m,N);
    diagb          = ones(N,1);
    diagb(~hp.iRW) = 0;
    b(1:N,:)       = diag(diagb);    %constant and exogenous in last columns

 
    % Parameters of the hyperpriors, if selected
    if hyperpriors

        % Define modes of the hyperpriors
        mode.l1 = 0.2;                      %0.4 Miranda & Ricco; 0.2 Sims
        mode.l3 = 1;
        mode.l4 = 1;

        % Define standard deviations of the hyperpriors
        sd.l1   = 0.4;                      %0.2 Miranda & Ricco; 0.4 Sims
        sd.l3   = 1;
        sd.l4   = 1;

        % Define scale and shape of the IG on SS2/(d-n-1)
        scaleSS2 = 0.02^2;        

        % Compute scale and shape of the G on l1, l3 and l4
        priorcoef.l1 = GammaCoef(mode.l1,sd.l1,0);
        priorcoef.l3 = GammaCoef(mode.l3,sd.l3,0);
        priorcoef.l4 = GammaCoef(mode.l4,sd.l4,0);
        
        % Set scale (beta) and shape (alpha) of the G on SS2
        priorcoef.alpha.SS2 = scaleSS2;
        priorcoef.beta.SS2  = scaleSS2;
        
    else
        priorcoef = [];
    end

    % Maximise posterior of hyperparameters
    [csmin.fh,xh,csmin.gh,csmin.H,csmin.itct,csmin.fcount,csmin.retcodeh] = csminwel(...
        'logMLVAR_formin_mod',...   %objective function
        X0,...                  %initial values
        H0,...                  %initial value for PD inverse hessian
        [],...                  %opt for numerical gradient
        1e-16,...               %convergence criterion
        1000,...                %max number of iterations
        ...
        ...                     %------- Options for obj. function -------%
        Y,...                   %Y
        X,...                   %X
        p,...                   %lags
        m,...                   %number of exogenous variables
        T-p,...                 %length of Y
        N,...                   %number of endogenous vars
        b,...                   %Minnesota prior mean
        MIN,...                 %lower bounds
        MAX,...                 %upper bounds
        hp.SS2,...              %residual variance of univariate AR(p)
        hp.l2,...               %default for lag-decay
        hp.l5,...               %prior variance for constant and exogenous
        hp.iRW,...              %position of ~iRW
        iPrior,...              %position of selected priors
        Ybar,...                %mean of first p lags
        Xbar,...                %mean of first p lags (exogenous)
        hyperpriors,...         %switch: priors on hyperparameters
        priorcoef);             %if hyperpriors is on

    % Position of parameters in par
    p1 = strcmp(priorlist(1),priors);
    p2 = strcmp(priorlist(2),priors);
    p3 = strcmp(priorlist(3),priors);
    p4 = strcmp(priorlist(4),priors);
    p5 = strcmp(priorlist(5),priors);
    
    % Posterior mode of priors
    out.l1 = MIN.l1+(MAX.l1-MIN.l1)./(1+exp(-xh(p1)));
    
    if ~isempty(find(p2,1))
    out.l2 = MIN.l2+(MAX.l2-MIN.l2)./(1+exp(-xh(p2)));
    end
        
    if ~isempty(find(p3,1))
    out.l3 = MIN.l3+(MAX.l3-MIN.l3)./(1+exp(-xh(p3)));
    end
    
    if ~isempty(find(p4,1))
    out.l4 = MIN.l4+(MAX.l4-MIN.l4)./(1+exp(-xh(p4)));
    end
    
    if ~isempty(find(p5,1))
    out.SS2 = MIN.SS2+(MAX.SS2-MIN.SS2)./(1+exp(-xh(find(p5,1):end)));
    out.SS  = sqrt(out.SS2);
    end
    
    % VAR coefficients and residuals at the posterior mode
    [out.fh,postMode.beta,postMode.sigma] = logMLVAR_formin_mod(xh,Y,X,p,m,T-p,N,b,MIN,MAX,hp.SS2,hp.l2,hp.l5,hp.iRW,iPrior,Ybar,Xbar,hyperpriors,priorcoef);

    % NOTE: this line can be removed. It is just to check that logMLVAR_formin_mod2
    % gives the same results as the two BVAR codes that follows.


