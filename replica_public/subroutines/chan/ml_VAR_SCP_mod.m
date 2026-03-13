% This script computes the log marginal likelihood of the BVAR with the
% asymmetric conjugate prior
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

% MODIFIED: r.degasperi.ac.uk (10/12/2019)

function lml = ml_VAR_SCP_mod(p,m,Y,Z,hp,k)
    
    hp.l1 = exp(k);  %tightness

    [T,n] = size(Y);
    lml = -n*T/2*log(2*pi);
    for i = 1:n    
        
        yi = Y(:,i);
        ki = n*p+i +m-1;
        
        % Prior
        [mi,Vi,nui,Si,~] = prior_ACPi_mod(n,p,m,i,hp);
        Xi = [Z -Y(:,1:i-1)];
        
        iVi = Vi\speye(ki);                                                       %inverse prior variance for equation i
        Kthetai = iVi + Xi'*Xi;
        CKthetai = chol(Kthetai,'lower');    
        thetai_hat = CKthetai'\(CKthetai\(iVi*mi + Xi'*yi));                      %draw from posterior of beta for eq i
        Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;   %draw from posterior of sigma for eq i
        
        lml = lml -1/2*(sum(log(diag(Vi))) + 2*sum(log(diag(CKthetai)))) ...
            + nui*log(Si) - (nui+T/2)*log(Si_hat) + gammaln(nui+T/2) - gammaln(nui);        
    end    
end