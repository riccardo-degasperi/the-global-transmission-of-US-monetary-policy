function out = loglik(Y,X,beta,sigma)
%loglik computes the log-likelihood of a VAR.

% INPUT:
% Y     : TxN matrix
% X     : Tx(N*p+m) matrix
% beta  : Nx(N*p+m) matrix
% sigma : NxN matrix

%-------------------------------------------------------------------------%

e     = Y - X*beta;
iS    = inv(sigma);
dS    = log(det(iS));
[T,N] = size(Y);
sterm = 0;

for i = 1:T
    sterm = sterm + (e(i,:)*iS*e(i,:)');
end

%out = (-(T*N)/2)*log(2*pi) + (T/2)*dS - 0.5*sterm;
out = (T/2)*dS - 0.5*sterm;