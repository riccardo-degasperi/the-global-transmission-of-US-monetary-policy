% This script obtains the optimal kappa1 and kappa2 by maximizing the
% marginal likelihood
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

% MODIFIED: r.degasperi.ac.uk (10/12/2019)

function [ml_opt,hp] = get_OptKappa_mod(Y,X,p,m,hp)

    if hp.cross %if optimisation of cross-variable shrinkage is active

        f = @(k) -ml_VAR_ACP_mod(p,m,Y,X,hp,[k(1) k(2)]);
        [k_opt,nml] = fminsearch(f,[log(hp.l1) log(hp.l6)]);
        ml_opt = -nml;

        hp.l1 = exp(k_opt(1));
        hp.l6 = exp(k_opt(2));

    else

        f = @(k) -ml_VAR_SCP_mod(p,m,Y,X,hp,k);
        [k_opt,nml] = fminsearch(f,log(hp.l1));
        ml_opt = -nml;

        hp.l1 = exp(k_opt(1));
        
    end
end