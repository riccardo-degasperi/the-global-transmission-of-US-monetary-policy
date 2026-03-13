% This script constructs the asymmetric conjugate prior given the
% hyperparameter values
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

% MODIFIED: r.degasperi.ac.uk (10/12/2019)

function [mi,Vi,nui,Si,sto_Vi] = prior_ACPi_mod(N,p,m,var_i,hp)

%unpack
sig2   = hp.SS2;
iRW    = hp.iRW;
iEq    = hp.iEq;
iParam = hp.iParam;

l1 = hp.l1;  %own-lags
l2 = hp.l2;  %lag-decay
l5 = hp.l5;  %intercept and exogenous
l6 = hp.l6;  %cross-variable extra shrinkage

eps = hp.eps;  %dogmatic shrinkage parameter

kappa3 = 1;  %not optimised yet

ki = var_i + N*p +m-1;
mi = zeros(ki,1); if iRW(var_i) == 1; mi(var_i) = 1; end
Vi = zeros(ki,1);

% construct Vi
for j=1:ki
    if j <= N*p+m
        l = ceil(j/N);   %lag length
        idx = mod(j,N);  %variable index
        if idx==0
            idx = N;
        end    
    else    
        idx = j - (N*p+1);
    end
    
    if j>=N*p+1 && j<=N*p+m % intercept
        Vi(j) = (l1*l5)^2;
    elseif j > N*p+m    % alpha_i
        Vi(j) = kappa3/sig2(idx);   
    elseif idx == var_i % own lag
        Vi(j) = (1/sig2(idx))*(l1/(l^l2))^2;
    else % lag of other variables
        
        if ismember(var_i,iEq) && ismember(idx,iParam)
            Vi(j) = eps;
        else
            Vi(j) = (1/sig2(idx))*(l1*l6/(l^l2))^2;
        end
        
    end    
end
Si = sig2(var_i)/2; 
nui = 1 + var_i/2;
sto_Vi = Vi;
Vi = sparse(1:ki,1:ki,Vi);
end
