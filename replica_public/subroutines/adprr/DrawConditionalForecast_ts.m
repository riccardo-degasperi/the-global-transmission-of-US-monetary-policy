function out = DrawConditionalForecast_ts(bs,Ms,SSArest,dims)


% Adapted from Antolin-Diaz, Petrella, Rubio-Ramirez (2021)
%
% Last modified: 29/02/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpack
yCondition   = SSArest.yCondition;
Omega_f      = SSArest.Omega_f;
eCondition   = SSArest.eCondition;
Omega_g      = SSArest.Omega_g;
h            = dims.h+1;
N            = dims.N;
draws        = size(bs,2);

% Initialise
stoY_cf      = nan(h,N,draws);
stoY_bl      = nan(h,N,draws);
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
    % Find type of conditional forecast
    if ~isempty(yCondition) && isempty(eCondition)
        type = 'condition-on-observables';
    elseif isempty(yCondition) && ~isempty(eCondition)
        type = 'condition-on-shock';   
    elseif ~isempty(yCondition) && ~isempty(eCondition)
        type = 'condition-on-structural-scenario';
    else
        error('ERROR: customise SSArestrictions.m')
    end

    %---------------------------------------------------------------------%
    % Calculate objects depending on conditions
    switch type
        case 'condition-on-observables'
    
            y_cf = vec(yCondition');
            C = diag(isfinite(y_cf));
            C(isnan(y_cf),:) = [];
            
            D = C*M';
            if isempty(Omega_f) 
                Omega_f = D*D';
            end
            
            f = y_cf(isfinite(y_cf));
            Omega_f_final = Omega_f;
            
        case 'condition-on-shock'
            shocks = vec(eCondition');
            Xi = diag(isfinite(shocks));
            Xi(isnan(shocks),:) = [];
            C = Xi/M';
            D = C*M';
            
            shocks(isnan(shocks))=0;
            f = C*b + Xi*shocks;
            Omega_f_final = Omega_g;
            
        case 'condition-on-structural-scenario'
            
            y_cf = vec(yCondition');
            f_ub = y_cf(isfinite(y_cf));
            C_ub = diag(isfinite(y_cf));
            C_ub(isnan(y_cf),:) = [];
            
            D_ub = C_ub*M';
            if isempty(Omega_f)
                Omega_f = D_ub*D_ub';
            end
            
            
            shocks = vec(eCondition');
            Xi = diag(isfinite(shocks));
            Xi(isnan(shocks),:) = [];
            
            C_lb = Xi/(M');
            C = [C_ub ; C_lb];
            D = C*M';
            
            f = [f_ub ; C_lb*b];
            Omega_f_final = blkdiag(Omega_f,Omega_g);
            
    end

    %---------------------------------------------------------------------%
    % Calculate Mean and Variance of conditional and unconditional forecasts
    Dstar = pinv(D);  % Dstar is the generalised inverse of D
    Dhat = null(D)';  % The rows of Dhat form an orthonormal basis for the null space of D
    
    mu_y_bl    = b;
    Omega_y_bl = (M'*M);
    
    mu_y_cf    = M'*Dstar*f+M'*(Dhat'*Dhat)/(M')*b; % Equation (10)
    Omega_y_cf = M'*(Dstar*Omega_f_final*Dstar'+(Dhat'*Dhat))*M; % Equation (10)
    Omega_y_cf = nspd(Omega_y_cf); %(Omega_y_cf+Omega_y_cf')*.5; 
    
    mu_e       = Dstar*f - Dstar*C*b;
    Omega_e    = Dstar*Omega_f_final*Dstar'+Dhat'*Dhat;
    Omega_e    = nspd(Omega_e); %(Omega_e+Omega_e')*.5;

    %---------------------------------------------------------------------%
    % Draw objects from conditional Posterior Distribution 
    y_cf = mvnrnd(mu_y_cf,Omega_y_cf)';
    y_bl = mvnrnd(mu_y_bl,Omega_y_bl)';
    
    shocks = ((y_cf' - b')/(M))';
    shocks = reshape(shocks,[N,h])';

    y_cf = reshape(y_cf,[N,h])';
    y_bl = reshape(y_bl,[N,h])';

    %---------------------------------------------------------------------%
    % Kullback Leibler Divergence
    KL  = 0.5*(trace(Omega_e) + mu_e'*mu_e - N*h + log(1/det(Omega_e)));
    KLc = (1 + (1 - exp(-2*KL/(N*h)))^(0.5))/2;

    %---------------------------------------------------------------------%
    % Store results
    stoY_cf(:,:,s)  = y_cf;
    stoY_bl(:,:,s)  = y_bl;
    stoKL(s,1)      = KL;
    stoKLc(s,1)     = KLc;
    stoShock(:,:,s) = shocks;

end

% Pack output
out.y_cf   = stoY_cf;
out.y_bl   = stoY_bl;
out.KL     = stoKL;
out.KLc    = stoKLc;
out.shocks = stoShock;

