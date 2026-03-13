function out = DrawConditionalForecast_irf_mod(bs,Ms,SSArest,iImpulse,dims)


% This function draws a conditional forcast from the truncated distribution
% according the variables conditioned upon. The function encompassed three
% types of conditioning: Conditioning on a path of one or more observables;
% Conditioning on a path of one or more shocks; or Conditioning on a path
% of both shocks and obserables.



% Adapted from Antolin-Diaz, Petrella, Rubio-Ramirez (2021) and from
% Breitenlechner, Georgiadis, Schumann (2022)
%
% Last modified: 29/02/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Additional options
version = 'ADPRR'; %either 'ADPRR' or 'BGS'


%-------------------------------------------------------------------------%

% Unpack
yCondition   = SSArest.yCondition;
Omega_f      = SSArest.Omega_f;
eCondition   = SSArest.eCondition;
Omega_g      = SSArest.Omega_g;
h            = dims.h+1;             %+1 for consistency with IRFbuild function
N            = dims.N;
draws        = size(bs,2);

% Initialise
y_cf         = nan(h,N,draws);
y_bl         = nan(h,N,draws);
stoKL        = nan(draws,1);
stoKLc       = nan(draws,1);
stoShock     = nan(h,N,draws);
w = 1;

for s = 1:draws

    % Display advancement
    if s == 100*w
        if exist('txt','var')
            fprintf(repmat('\b', 1, numel(txt)-1));
        end
        txt = ['       Drawing sequence: ',num2str(s),'/',num2str(draws),'\n'];
        fprintf(txt);
        w = w+1;
    end

    % Select objects for current draw
    b = bs(:,s);
    M = Ms(:,:,s);
 

    %---------------------------------------------------------------------%
    % Objects conditional shock only

    eCondition_bl = zeros(h,N);                % all shocks are restricted to unconditional distribution
    eCondition_bl(1,iImpulse) = 1;             % only shock of interest has impulse
    shocks2 = vec(eCondition_bl');             %nh x 1
    Xi2 = diag(isfinite(shocks2));
    Xi2(isnan(shocks2),:) = [];                %ks x nh
    C_lb2 = Xi2/(M');                          %unused
    shocks2(isnan(shocks2)) = 0;
    
    % Keep only the stochastic part
    %f_tplus1_tplush2 = C_lb2*b_tplus1_tplush + Xi2*e_tplus1_tplush2; %ks x 1
    f_tplus1_tplush2 = Xi2*shocks2;
    
    C2 = C_lb2;                                %unused
    D2 = double(Xi2);
     
    Omega_f_final2 = zeros(N*h);               %equal to Omega_g

    %---------------------------------------------------------------------%
    % Calculate Mean and Variance of unconditional forecasts

    switch version

        case 'ADPRR'

        % BASELINE SCENARIO (IMPULSE ONLY)
        Dstar2 = pinv(D2);   % Dstar is the generalised inverse of D
        Dhat2  = null(D2)';  % The rows of Dhat form an orthonormal basis for the null space of D

        % ---Mean under baseline scenario---
        mu_y_unc = M'*Dstar2*f_tplus1_tplush2+M'*(Dhat2'*Dhat2)/(M')*b;
        
        % ---Variance under baseline scenario---
        Omega_y_unc = M'*(Dstar2*Omega_f_final2*Dstar2'+(Dhat2'*Dhat2))*M;
        Omega_y_unc = (Omega_y_unc+Omega_y_unc')*.5;                              %for numerical precision

        case 'BGS' % Breitenlechner, Georgiadis, Schumann (2022)

        % BASELINE SCENARIO (IMPULSE ONLY)
        Dstar2 = pinv(D2);   % Dstar is the generalised inverse of D

        % Mean under baseline scenario
        mu_y_unc = b + M'*Dstar2*(f_tplus1_tplush2-C2*b);

        % Variance under baseline scenario (with correct sign!)
        Omega_y_unc = M'*M + M'*Dstar2*(Omega_f_final2-D2*D2')*Dstar2'*M;
        Omega_y_unc = nspd(Omega_y_unc);
            
    end

    %---------------------------------------------------------------------%
    % Draw objects from unconditional Posterior Distribution 
    irf_bl = mvnrnd(mu_y_unc,Omega_y_unc)';
    irf_bl = reshape(irf_bl,[N,h])';


    %---------------------------------------------------------------------%
    % IF XXXX, select trajectory for policy variable
    
    traj = irf_bl(:,iImpulse);


    %---------------------------------------------------------------------%
    % Objects conditional on impulse + scenario

    % Restrictions on observables
    irf_cf = vec(yCondition');                        %nh x 1
    C_ub = diag(isfinite(irf_cf));                    %nh x nh
    C_ub(isnan(irf_cf),:) = [];                       %ko x nh
    f_ub = irf_cf(isfinite(irf_cf));                  %ko x 1
    
    % Restrictions on shocks
    shocks_cf = vec(eCondition');                     %nh x 1
    Xi = diag(isfinite(shocks_cf));                   %nh x nh
    Xi(isnan(shocks_cf),:) = [];                      %ks x nh
    C_lb = Xi/(M');                                   %ks x nh
    shocks_cf(isnan(shocks_cf)) = 0; %also the unrestricted shocks have mean zero

    % Keep only the stochastic part (this is g_tplus1_tplush in Jorgo's code)
    %f_lb = C_lb*b + Xi*shocks_cf;                     %ks x 1
    f_lb = Xi*shocks_cf;
    
    %RD: it's not necessary to do this for f_tplus1_tplush_ub because we
    %specify the path of the chosen observable to zero (i.e. to steady state).

    % Stack restrictions on observables and shocks
    C = [C_ub ; C_lb];
    D = [C_ub*M' ; Xi];
    
    % If you choose "loose" restrictions on oservables:
    if isempty(Omega_f)
        D_upper_bar = C_ub*M';
        Omega_f = D_upper_bar*D_upper_bar';
    end  
    
    f_tplus1_tplush = [f_ub ; f_lb];
    Omega_f_final = blkdiag(Omega_f,Omega_g);


    %---------------------------------------------------------------------%
    % Calculate Mean and Variance of conditional and unconditional forecasts

    switch version

        case 'ADPRR'

        % % BASELINE SCENARIO (IMPULSE ONLY)
        % Dstar2 = pinv(D2);   % Dstar is the generalised inverse of D
        % Dhat2  = null(D2)';  % The rows of Dhat form an orthonormal basis for the null space of D
        % 
        % % ---Mean under baseline scenario---
        % % ADPRR original
        % mu_y_unc = M'*Dstar2*f_tplus1_tplush2+M'*(Dhat2'*Dhat2)/(M')*b;
        % % Jorgo Eq D.15
        % %mu_y_unc_alt = b + M'*Dstar2*(f_tplus1_tplush2-C2*b);
        % % Note: they are identical.
        % 
        % % ---Variance under baseline scenario---
        % % ADPRR original
        % Omega_y_unc = M'*(Dstar2*Omega_f_final2*Dstar2'+(Dhat2'*Dhat2))*M;
        % Omega_y_unc = (Omega_y_unc+Omega_y_unc')*.5;                              %for numerical precision
        % % Jorgo Eq D.16 (with correct sign!)
        % %Omega_y_unc_alt = M'*M + M'*Dstar2*(Omega_f_final2-D2*D2')*Dstar2'*M;
        % %Omega_y_unc_alt = nspd(Omega_y_unc_alt); %(Omega_y_unc_alt+Omega_y_unc_alt')*.5;                  %for numerical precision
        % % Note: they are virtually identical
        % 
        % % % Same in ADPRR and Jorgo Eq D.12
        % % mu_e_unc = Dstar2*f_tplus1_tplush2 - Dstar2*C2*b_tplus1_tplush;
        % % % ADPRR original
        % % Omega_e_unc = Dstar2*Omega_f_final2*Dstar2'+Dhat2'*Dhat2;
        % % Omega_e_unc = (Omega_e_unc+Omega_e_unc')*.5; 
        % % % Jorgo Eq D.13
        % % %Omega_e_unc_alt = Dstar2*Omega_f_final2*Dstar2'+(eye(n*h)-(Dstar2*D2*D2'*Dstar2'));
        % % %Omega_e_unc_alt = (Omega_e_unc_alt+Omega_e_unc_alt')*.5; 
        % % % Note: they are identical
        
    
        % COUNTERFACTUAL SCENARIO
        Dstar  = pinv(D);   % Dstar is the generalised inverse of D
        Dhat   = null(D)';  % The rows of Dhat form an orthonormal basis for the null space of D
        
        % ---Mean under structural scenario---
        % ADPRR original
        mu_y = M'*Dstar*f_tplus1_tplush+M'*(Dhat'*Dhat)/(M')*b;
        % Jorgo Eq D.15
        %mu_y_alt = b + M'*Dstar*(f_tplus1_tplush-C*b);
        % Note: they are virtually identical
        
        % ---Variance under structural scenario---
        % ADPRR original
        Omega_y = M'*(Dstar*Omega_f_final*Dstar'+(Dhat'*Dhat))*M;
        Omega_y = (Omega_y+Omega_y')*.5;                                          %for numerical precision
        % Jorgo Eq D.16 (with correct sign!)
        %Omega_y_alt = M'*M + M'*Dstar*(Omega_f_final-D*D')*Dstar'*M;
        %Omega_y_alt = nspd(Omega_y_alt);% (Omega_y_alt+Omega_y_alt')*.5;                              %for numerical precision
        % Note: they are virtually identical
        
        % % Same in ADPRR and Jorgo Eq D.12
        % mu_e = Dstar*f_tplus1_tplush - Dstar*C*b_tplus1_tplush;
        % % ADPRR original
        % Omega_e = Dstar*Omega_f_final*Dstar'+Dhat'*Dhat;
        % Omega_e = (Omega_e+Omega_e')*.5; 
        % % Jorgo Eq D.13
        % %Omega_e_alt = Dstar*Omega_f_final*Dstar'+(eye(n*h)-(Dstar*D*D'*Dstar'));
        % %Omega_e_alt = (Omega_e_alt+Omega_e_alt')*.5; 
        % % Note: they are virtually identical

        case 'BGS' % Breitenlechner, Georgiadis, Schumann (2022)

        % % BASELINE SCENARIO (IMPULSE ONLY)
        % Dstar2 = pinv(D2);   % Dstar is the generalised inverse of D
        % 
        % % Mean under baseline scenario
        % mu_y_unc = b + M'*Dstar2*(f_tplus1_tplush2-C2*b);
        % 
        % % Variance under baseline scenario (with correct sign!)
        % Omega_y_unc = M'*M + M'*Dstar2*(Omega_f_final2-D2*D2')*Dstar2'*M;
        % Omega_y_unc = nspd(Omega_y_unc);
    
        % COUNTERFACTUAL SCENARIO
        Dstar  = pinv(D);   % Dstar is the generalised inverse of D

        % Mean under structural scenario
        mu_y = b + M'*Dstar*(f_tplus1_tplush-C*b);
        
        % Variance under structural scenario (with correct sign!)
        Omega_y = M'*M + M'*Dstar*(Omega_f_final-D*D')*Dstar'*M;
        Omega_y = nspd(Omega_y);
        
    end
    
    
    %---------------------------------------------------------------------%
    % Draw objects from conditional Posterior Distribution 
    irf_cf = mvnrnd(mu_y,Omega_y)';
    % irf_bl = mvnrnd(mu_y_unc,Omega_y_unc)';
    
    shocks = (irf_cf'/M)';
    shocks = reshape(shocks,[N,h])';
    
    irf_cf = reshape(irf_cf,[N,h])';
    % irf_bl = reshape(irf_bl,[N,h])';
    
    
    %---------------------------------------------------------------------%
    % Kullback Leibler Divergence
    % RD: KL divergence is not deﬁned in case the baseline and the
    % counterfactual do not feature any uncertainty. In that case we need to
    % use the modified version proposed in Breitenlechner, Georgiadis & Schumann (2022).
    % See the Appendix of that paper.

    %KL  = 0.5*(trace(Omega_e) + mu_e'*mu_e - n*h + log(1/det(Omega_e)));
    %KLc = (1 + (1 - exp(-2*KL2/(n*h)))^(0.5))/2;
    
    
    %---------------------------------------------------------------------%
    % IMPLEMENT JORGO's VERSION OF KL (appendix of What goes around comes
    % around)
    
    % Impose a path on all variables at all horizons
    C_ub_kl = eye(N*h);   
    
    % Average paths for observables under bl and cf
    % f_bl = mu_y_unc; % = M'*e_tplus1_tplush2;
    % f_cf = mu_y;
    f_bl = mu_y_unc; % = M'*e_tplus1_tplush2;
    f_cf = mu_y;
    
    % Applying "loose" restrictions on oservables
    Omega_f_kl = M'*M;
    % NOTE: D = M', since Xi = 0 and C_lb = 0
    %       Xi = 0 means that all shocks can deviate from their unconditional distribution
    %       C_lb = Xi/M' = 0 means that no observable is driven by a specific path of the structural shocks
    
    % As a consequence:
    Sigma_e_kl = eye(N*h);
    mu_e_bl = inv(M')*f_bl;
    mu_e_cf = inv(M')*f_cf;
    Sigma_y_kl = M'*M;
    
    % Simplified Eq. D.32 in Jorgo's paper
    KL = 0.5*( (mu_e_cf-mu_e_bl)'*(mu_e_cf-mu_e_bl) );
    
    % Calibrate it so that is bound between 0.5 and 1
    KLc = (1 + (1 - exp(-2*KL/(N*h)))^(0.5))/2;
    % NOTE: closer to 0.5 means more plausible


    %---------------------------------------------------------------------%
    % Store results
    y_cf(:,:,s)     = irf_cf;
    y_bl(:,:,s)     = irf_bl;
    stoKL(s,1)      = KL;
    stoKLc(s,1)     = KLc;
    stoShock(:,:,s) = shocks;

end

% Pack output
out.y_cf   = y_cf;
out.y_bl   = y_bl;
out.KL     = stoKL;
out.KLc    = stoKLc;
out.shocks = stoShock;

