function [logML,beta,sigma]=logMLVAR_formin_mod(par,Y,X,p,m,T,N,b,MIN,MAX,s20,l20,l5,iRW,iPrior,Ybar,Xbar,hyperpriors,priorcoef)

% This function computes the log-posterior (or the logML if hyperpriors=0), 
% the posterior mode of the coefficients and the covariance matrix of the residuals of the BVAR of 
% Giannone, Lenza and Primiceri (2012)
%
% Last modified: Riccardo Degasperi 28/11/2018



% SHOULD WORK WITH m > 1 BUT CHECK THE DoF!!!!!!!!

%-------------------------------------------------------------------------%

% List of possible priors
priorlist = {'niw','lag','soc','cop','std'};
priors    = priorlist(iPrior);                 %selected options

% Position of parameters in par
p1 = strcmp(priorlist(1),priors);
p2 = strcmp(priorlist(2),priors);
p3 = strcmp(priorlist(3),priors);
p4 = strcmp(priorlist(4),priors);
p5 = strcmp(priorlist(5),priors);

% DoF for prior IW distribution for sigma 
d = N+2;

% Update priors
l1 = MIN.l1+(MAX.l1-MIN.l1)./(1+exp(-par(p1)));
l2 = MIN.l2+(MAX.l2-MIN.l2)./(1+exp(-par(p2)));
l3 = MIN.l3+(MAX.l3-MIN.l3)./(1+exp(-par(p3)));
l4 = MIN.l4+(MAX.l4-MIN.l4)./(1+exp(-par(p4)));

if isempty(find(p5,1))
    s2 = s20*(d-N-1);     %default for s2
else
    s2 = MIN.SS2+(MAX.SS2-MIN.SS2)./(1+exp(-par(find(p5,1):end)));
end

if isempty(l2)
    l2 = l20;            %default for l2 (lag-decay)
end


%-------------------------------------------------------------------------%
% Estimate BVAR using GLP

% Minnesota covariance matrix
omega = zeros(N*p+m,1);
omega(N*p+1:end) = (l1*l5)^2;        %constant and exogenous on bottom
for i = 1:p
    %omega((i-1)*N+1:i*N)=(d-N-1)*(l1^2)*(1/(i^l2))./s2; % IN GLP THE DEFAULT l2 = 2
    omega((i-1)*N+1:i*N)=(d-N-1)*((l1/(i^l2))^2)./s2;
end

% prior scale matrix for the covariance of the shocks
PSI = diag(s2);

% Defaults if no soc or cop
Td = 0;
Xsoc = [];
Ysoc = [];
Xcop = [];
Ycop = [];
Yd = Y;
Xd = X;

% Generate dummy observations (constant and exogenous on bottom)
if ~isempty(l3)     %Sum-of-coefficients

    Ysoc = diag(Ybar.*iRW)/l3;
    Xsoc = [kron(ones(1,p),Ysoc) zeros(N,m)];

    Yd = [Yd;Ysoc];
    Xd = [Xd;Xsoc];
    Td = Td + N;
end

if ~isempty(l4)     %Copersistence

    Ycop = Ybar'/l4;
    Xcop = [kron(ones(1,p),Ycop) 1/l4 Xbar/l4];

    Yd = [Yd;Ycop];
    Xd = [Xd;Xcop];
    Td = Td + 1;
end

% Total number of observations
T = T + Td;

% Estimate BVAR
beta  = (Xd'*Xd + diag(1./omega))\(Xd'*Yd + diag(1./omega)*b);
e     = Y - X*beta;
sigma = (e'*e + PSI + (beta-b)'*diag(1./omega)*(beta-b))/(T+d+N+1);

%-------------------------------------------------------------------------%
% Compute log-posterior

aaa = diag(sqrt(omega))*(X'*X)*diag(sqrt(omega));
bbb = diag(1./sqrt(s2))*(e'*e + (beta-b)'*diag(1./omega)*(beta-b))*diag(1./sqrt(s2));

eigaaa = real(eig(aaa));
eigaaa(eigaaa<1e-12) = 0;
eigaaa = eigaaa+1;

eigbbb = real(eig(bbb));
eigbbb(eigbbb<1e-12) = 0;
eigbbb = eigbbb+1;

logML = - N*T*log(pi)/2 + sum(gammaln((T+d-[0:N-1])/2)-gammaln((d-[0:N-1])/2)) +...
        - T*sum(log(s2))/2 - N*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;

if ~isempty(l3) || ~isempty(l4)
    yd = [Ycop;Ysoc];
    xd = [Xcop;Xsoc];
    
    % prior mode of the VAR coefficients
    % betahatd=(xd'*xd+diag(1./omega))\(xd'*yd+diag(1./omega)*b);
    betad = b;     % this is the case for our priors (the line above delivers the same but is numerically not very stable)
    
    % VAR residuals at the prior mode
    ed = yd-xd*betad;
    
    aaa = diag(sqrt(omega))*(xd'*xd)*diag(sqrt(omega));
    bbb = diag(1./sqrt(s2))*(ed'*ed + (betad-b)'*diag(1./omega)*(betad-b))*diag(1./sqrt(s2));
    
    eigaaa = real(eig(aaa));
    eigaaa(eigaaa<1e-12) = 0;
    eigaaa = eigaaa+1;
    
    eigbbb = real(eig(bbb));
    eigbbb(eigbbb<1e-12) = 0;
    eigbbb = eigbbb+1;
    
    % normalizing constant
    norm = - N*Td*log(pi)/2 + sum(gammaln((Td+d-[0:N-1])/2)-gammaln((d-[0:N-1])/2)) +...
           - Td*sum(log(s2))/2 - N*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
    
    logML = logML-norm;
end


if hyperpriors
    
    logML = logML+logGammapdf(l1,priorcoef.l1.k,priorcoef.l1.theta);
    
    if ~isempty(l3)
        logML = logML+logGammapdf(l3,priorcoef.l3.k,priorcoef.l3.theta);
    end
    
    if ~isempty(l4)
        logML = logML+logGammapdf(l4,priorcoef.l4.k,priorcoef.l4.theta);
    end
    
    if ~isempty(find(p5,1))
        logML = logML+sum(logIG2pdf(l5/(d-N-1),priorcoef.alpha.SS2,priorcoef.beta.SS2));
    end
    
end

logML = -logML;

%-------------------------------------------------------------------------%
function r=logGammapdf(x,k,theta)                  %k:shape; theta:scale
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);    %log of gamma pdf

function r=logIG2pdf(x,alpha,beta)                          %alpha:shape; beta:scale
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-gammaln(alpha);  %log of IG pdf