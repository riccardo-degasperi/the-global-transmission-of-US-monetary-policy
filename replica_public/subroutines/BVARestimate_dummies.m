function out = BVARestimate_dummies(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt)
% Estimates the reduced-form BVAR using dummy observations.

% INPUT:
% - data     : output from makeXY;
% - dims     : output from makeXY;
% - hp       : hyperparameters and related options;
% - modelOpt : structure containing model options;
% - dataOpt  : structure containing data options; 
% - miscOpt  : structure containing misc options;
% - plotOpt  : structure containing plot options;

% OUTPUT:
% - out      : structure containing estimation results and IRFs;

% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%
disp('------------------------------------------')
disp('Estimating Bayesian VAR (dummy obs)')

% Unpacking
X               = data.X;
Y               = data.Y;
endo            = data.endo;
Xex             = data.Xex;
T               = dims.T;
N               = dims.N;
p               = dims.p;
pexo            = dims.pexo;
m               = dims.m;
mc              = dims.mc;
h               = dims.h;
BVARcase        = miscOpt.BVARcase;
dropExplosive   = miscOpt.dropExplosive;
sim             = modelOpt.sim;
burnin          = modelOpt.burnin;
jump            = modelOpt.jump;
draws           = (sim-burnin)/jump;
GLP             = hp.GLP;
plotReducedForm = plotOpt.plotReducedForm;
plotType        = plotOpt.plotType;
folder          = plotOpt.folder;
plotOpt.ylab    = dataOpt.ylab;
varendo         = dataOpt.varendo;
varendolong     = dataOpt.varendolong;
shockVar        = dataOpt.shockVar;
shockSign       = dataOpt.shockSign;
units           = dataOpt.units;

% Patch for chartNames when units are more than 1
if numel(units) > 1
    units0 = join(string(units),'_');
else
    units0 = units;
end

%-------------------------------------------------------------------------%
% Estimate univariate AR(order) for each endogenous variable
hp.SS  = NaN(N,1);         %container for residual standard deviations
hp.SS2 = NaN(N,1);         %container for residual variances
order  = 1;                %we fit an AR(order) for the whole sample T

for i = 1:N

    % Build X and Y matrices
    X2  = [ones(T-order,1) lagX(endo(:,i),order)]; % (T-p)x(p+1) matrix (with constant)
    Y2  = endo(order+1:end,i);                     % (T-p)x1 matrix
    
    % Intermediate step
    XX2 = X2'*X2;
    XY2 = X2'*Y2;

    % Compute beta and std of residuals
    b2        = XX2\XY2;
    hp.SS(i)  = std(Y2-X2*b2);
    hp.SS2(i) = hp.SS(i).^2;
    
end

% Compute mean across first p observations for endogenous
Ybar = mean(endo(1:p,:),1)';

% Compute mean across first p observations for exogenous
if ~isempty(Xex)
    Xbar = mean(Xex(1:p-pexo,:),1);
else
    Xbar = [];
end


%-------------------------------------------------------------------------%
% Optimal priors

hpVAR = hp;   %initialise hpVAR

if GLP
disp('--> Computing optimal priors...')
    
    [hpglp,csmin,postMode] = GLPoptimalPriors(Y,X,Ybar,Xbar,hp,dims);

    % Update hyperparameters
    f     = fieldnames(hpglp);
    for i = 1:length(f)
        hpVAR.(f{i}) = hpglp.(f{i});
    end

end


%-------------------------------------------------------------------------%
% Estimate reduced-form BVAR

% Generate dummy observations
[Yd,Xd] = dummyObs(Ybar,Xbar,hpVAR,dims);

% Combined X and Y matrices
Yd = [Y; Yd];
Xd = [X; Xd];

% Estimate VAR
beta    = Xd\Yd;
ed      = Yd-Xd*beta;            %including dummies
[TT,NN] = size(Xd);
d       = N+2;
sigma   = (ed'*ed)/(TT-NN+d+1);  %same Dof correction as GLP


%-------------------------------------------------------------------------%
% Generate confidence bands for reduced-form BVAR
disp('--> Generating confidence bands for BVAR')

betas   = NaN(N*p+m+mc,N,draws);
sigmas  = NaN(N,N,draws);
es      = NaN(T-p,N,draws);
[TT,NN] = size(Xd);                             %(Td + T); N*p+m
neig    = 0;
maxeig  = NaN(sim,1);

k = 1;                                          %initialise index for storage matrices
s = 1;                                          %initialise count
w = 1;                                          %initialise count

CinvS  = chol(inv((TT-NN+d+1)*nspd(sigma)));    %removing DoF correction from sigma
CinvXX = chol(inv(Xd'*Xd));

while s <= sim

    % Display advancement
    if s == 100*w
        if exist('txt','var')
            fprintf(repmat('\b', 1, numel(txt)));
        end
        disp(' ')
        txt = ['    Gibbs Sampler #: ',num2str(s),'/',num2str(sim),'\n'];
        fprintf(txt);
        w = w+1;
    end

    %Get Psi
    Xg = randn(TT+2-(N*p+m+mc),N)*CinvS;           % vec(Xg) ~ N(vec(0),kron(inv(TT*sigma),I))
    sigmas_tmp = inv(Xg'*Xg);                   % (Xg'Xg) ~ W(inv(TT*sigma),TT+2-(N*p+1))
    % NOTE: Dof same as Banbura, Giannone & Reichlin (2010).
    %       TT+2-(N*p+m) = T+N+2 (= df in GLP)
    
    if mod(s,jump) == 0 && s > burnin 
    sigmas(:,:,k) = sigmas_tmp;                 %store only after jump
    end
    
    %Get Cholesky factor to use in MN below
    Csigmas = chol(sigmas_tmp);

    %Draw from posterior of beta (matrixvariate)
    tmp  = randn(size(beta));
    betas_tmp = beta + CinvXX'*tmp*Csigmas;     % vec(betag) ~ N(vec(beta),kron(Sg,invXX))
    
    % Check maximum eigenvalue
    A = NaN(N*p);                               %create companion matrix
    A(1:N,:) = betas_tmp(1:end-m-mc,:)';
    A(N+1:end,:) = [eye(N*(p-1)) zeros(N*(p-1),N)];
    maxeig(s,1) = max(abs(eig(A)));             %compute maximum eigenvalue
    if maxeig(s,1) > 1
        neig = neig+1;
    end
    
    % Burnin phase
    if s <= burnin
        s = s+1;
    
    % If we keep explosive draws
    elseif ~dropExplosive && s > burnin
        
        if mod(s,jump) == 0
        betas(:,:,k) = betas_tmp;               %VAR coefficients
        es(:,:,k)    = Y - X*betas_tmp;         %residuals for structural identification
        k = k+1;                                %update index for storage matrices
        end
        
        s = s+1;

    % If we discard explosive draws
    elseif dropExplosive && s > burnin && maxeig(s,1) < 1
    
        if mod(s,jump) == 0
        betas(:,:,k) = betas_tmp;               %VAR coefficients
        es(:,:,k)    = Y - X*betas_tmp;         %residuals for structural identification
        k = k+1;                                %update index for storage matrices
        end
        
        s = s+1;
        
    end
end

% Display number of explosive draws
if ~dropExplosive
disp(['    Number of explosive draws: ' num2str(neig)])
end
disp(' ')
    

%-------------------------------------------------------------------------%
% Tabulate VAR coefficients
disp('--> Tabulate VAR coefficients...')

% Get confidence regions
betas_l = quantile(betas(1:end-m-mc,:,:),0.05,3); %5% lower bound
betas_m = quantile(betas(1:end-m-mc,:,:),0.5,3);  %median
betas_u = quantile(betas(1:end-m-mc,:,:),0.95,3); %95% upper bound

% Build labels
tmp1 = char(repmat(varendo',p,1));
tmp2 = repmat(1:p,N,1); tmp2 = tmp2(:);
tmp3 = repmat(' p = ',N*p,1);
Regressors = tmp1;
Lags = [tmp3 num2str(tmp2)];

% Initialise table
TT = table(Regressors,Lags);

% Populate table iteratively
for z = 1:numel(varendo)
    TT.(['Eq. ',varendo{z}]) = [betas_l(:,z) betas_m(:,z) betas_u(:,z)];
end

% Display table
disp(TT)

% Write table as Excel file
folder_ = [pwd '/plots/',folder];
if ~isempty(folder)
    if not(isfolder(folder_))
        mkdir(folder_)
    end
    folder_ = [folder_,'/'];
end
fileName = [char(units0),'_DummyObs_coeffs.xlsx'];
writetable(TT,[folder_,fileName]);

%-------------------------------------------------------------------------%
% Generate Reduced-form Impulse Response Functions
disp('--> Generating reduced-form IRFs...')

IRF = IRFbuild(betas,dims,shockSign);

%-------------------------------------------------------------------------%
% Pack output
out.betas   = betas;
out.sigmas  = sigmas;
out.es      = es;
out.hpVAR   = hpVAR;
out.maxeig  = maxeig;
out.neig    = neig;
out.IRF     = IRF;
out.tabulatedCoeffs = TT;
