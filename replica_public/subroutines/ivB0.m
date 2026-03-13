function [B0,nnrB0,ncoB0,out] = ivB0(res,iv0,sigmas,dims,dataOpt)

% Unpack
T         = dims.T-dims.p;  %effective sample size
p         = dims.p;
m         = dims.m;
shockVar  = dataOpt.shockVar;
shockSize = dataOpt.shockSize;
shockIV   = dataOpt.shockIV;
varendo   = dataOpt.varendo;

% Retrieve dimensions
[Tiv,N,draws] = size(res);  %number of observations in the proxy
if isempty(draws)           %patch for when draws = 1
    draws = 1;
end
[~,Niv,~] = size(iv0);        %number of instruments
if ~isempty(shockIV)
    Nsiv = numel(shockIV);    %number of proxy-identified shocks
    Ns = numel(shockVar);     %total number of shocks
    iPiv = ismember(varendo,shockVar(shockIV));
else
    Nsiv = numel(shockVar);   %total number of shocks
    Ns = numel(shockVar);     %total number of shocks
    iPiv = ismember(varendo,shockVar);
end
iP   = ismember(varendo,shockVar);
if Nsiv>Niv
    error('ERROR: you need at least one proxy for each shock you are trying to identify.')
end

% Patch to compute normalised B0s even when shockSize = '1sd' 
if strcmp('1sd',shockSize)
    shockSize = 1;
end

% Place policy indicators first (works also for Niv>1)
res2 = [res(:,iPiv,:), res(:,~iPiv,:)];

% NOTE:
% - (A'A)^-1 * A'B = A\B when A is a M-by-N matrix with M>N.
% - Kronecker allows to run the regressions of each of the N
%   sets of residuals on the IV simultaneously.

% 2SLS
B0          = NaN(N,N,draws);            %normalised impact matrix
ncoB0       = NaN(draws,1);              %normalisation constant
nnrB0       = zeros(N,N,draws);          %non-normalised impact matrix
out.Fstat   = zeros(draws,Nsiv);         %first-stage F statistic
out.Frobust = zeros(draws,Nsiv);         %heteroscedasticity-robust F
out.RM      = zeros(Niv,Niv,draws);      %reliability matrix
out.corrIV  = zeros(draws,Niv);          %corr(IVs,shocks)
%out.sigma   = NaN(N,N,draws);            %proxy sample covariance matrix

for s = 1:draws

    % Patch for wild bootstrap (might not work with k>1)
    if size(iv0,3) == draws && draws > 1
    z = [ones(size(iv0(:,s))) iv0(:,:,s)];             %add constant
    iv = iv0(:,:,s);
    elseif size(iv0,3) == 1
    z = [ones(size(iv0,1),1) iv0];                   %add constant
    iv = iv0;
    else
        error('ERROR: Check the proxy in ivB0.')
    end

    % First stage (regress policy equations residuals on IVs)
    uhat = nan(Tiv,Nsiv);
    for is = 1:Nsiv

        % Regress each residual on all IVs and store fitted values
        firstStage = z\res2(:,is,s);
        uhat(:,is) = z*firstStage;

        % F stat (regression of relevant innovations on IVs)
        RSS = (res2(:,is,s)-uhat(:,is))'*(res2(:,is,s)-uhat(:,is));
        TSS = (res2(:,is,s)-repmat(mean(res2(:,is,s)),Tiv,1))'*(res2(:,is,s)-repmat(mean(res2(:,is,s)),Tiv,1));
        R2  = 1-RSS/TSS;
        Fstat = (R2/Niv) / ((1-R2)/(Tiv-Niv-1));
        %Fstat = ((TSS-RSS)/Niv)/(RSS/(Tiv-Niv-1)); %alternative
        out.Fstat(s,is) = Fstat;
        
        % Heteroscedasticity-robust F (see Wooldridge 2002 p. 55)
        R = [zeros(Niv,1) eye(Niv)];  %all coefficients are
        q = zeros(Niv,1);             %jointly equal to zero
        iXX = inv(z'*z);
        Shat = zeros(Niv+1);
        rr = res2(:,is,s)-uhat(:,is);
        for ii = 1:Tiv
            Shat = Shat+rr(ii)^2*z(ii,:)'*z(ii,:);
        end
        Frobust = 1/Niv*(R*firstStage-q)'*inv(Tiv/(Tiv-Niv-1)*R*iXX*Shat*iXX*R')*(R*firstStage-q);    % large sample (Wald test)
        out.Frobust(s,is) = Frobust;

    end

    % Second stage (regress remaining residuals on fitted)
    secondStage = [ones(Tiv,1) uhat]\res2(:,Nsiv+1:end,s);
    
    % Retrieve normalisation constant (Mertens & Ravn 2013)
    %sigma0 = (res2(:,:,s)'*res2(:,:,s))/(Tiv-N*p-m);       %problematic if Tiv < N*p+m 
    sigma0 = sigmas(:,:,s);                                 %for these sigmas T-p < N*p+m never binds

    % Reorder sigma
    sigma = [sigma0(iPiv,iPiv) sigma0(iPiv,~iPiv); sigma0(~iPiv,iPiv) sigma0(~iPiv,~iPiv)];
    
    sig11   = sigma(1:Nsiv,1:Nsiv);
    sig21   = sigma(Nsiv+1:N,1:Nsiv);
    sig22   = sigma(Nsiv+1:N,Nsiv+1:N);
    b21ib11 = secondStage(2:end,:)';
    Z       = b21ib11*sig11*b21ib11' - (sig21*b21ib11'+b21ib11*sig21') + sig22;
    b12b12p = (sig21 - b21ib11*sig11)'/Z*(sig21 - b21ib11*sig11);
    b11b11p = sig11 - b12b12p;
    
    % If you are identifying only one shock
    if Nsiv == 1

        b11 = sqrt(b11b11p);
        b1 = [b11; b21ib11*b11];
        b1unit = [1; b21ib11]*shockSize;

        % Store normalisation constant
        ncoB0(s) = b11;

    % If instead you have k>1 instruments k>1 shocks
    else

        b22b22p = sig22 + b21ib11*(b12b12p-sig11)*b21ib11';
        b12ib22 = (b12b12p*b21ib11' + (sig21-b21ib11*sig11)')/b22b22p;
        ib11iS1 = eye(Nsiv) - b12ib22*b21ib11;
        b21iS1  = b21ib11/ib11iS1;
        S1S1p   = ib11iS1*b11b11p*ib11iS1';

        % Additional restriction to separate the identified shocks
        S1 = chol(S1S1p,'lower');
        % NOTE: you can use sign restrictions instead, or impose any
        % other restriction.

        % Obtain first k columns of impact matrix
        b1 = [inv(ib11iS1); b21iS1]*S1;

        % Normalize b1 relative to the variable in ShockVar that is placed first in varendo
        b1unit = b1./b1(1,:)*shockSize;

    end

    % Proxy-sample covariance matrix
    %out.sigma(:,:,s) = sigma0;

    % Generate impact matrix (reordering the variables)
    B0(:,:,s)      = zeros(N,N);
    B0(iPiv,iPiv,s)  = b1unit(1:Nsiv,:);
    B0(~iPiv,iPiv,s) = b1unit(Nsiv+1:end,:);

    % Store non-normalised impact matrix
    nnrB0(:,:,s)      = zeros(N,N);
    nnrB0(iPiv,iPiv,s)  = b1(1:Nsiv,:);
    nnrB0(~iPiv,iPiv,s) = b1(Nsiv+1:end,:);
    
            
    % Reliability matrix (Mertens & Ravn eqs. A4 and A5)
    if Nsiv == 1 && Niv == 1

        d = eye(Nsiv)*sum(sum(iv,2)~=0)/Tiv;
        Bi = [1/b11+(sig21-sig11*b21ib11)'*inv(Z)/b11*b21ib11 -(sig21-sig11*b21ib11)'*inv(Z)/b11];
        et = (Bi*res2(:,:,s)')';
        sigmu = cov(iv,res2(:,1,s),1);
        PHI = sigmu(1,2)/b11;
        Gamma = d\PHI;
        E  = Gamma*et(sum(iv,2)~=0);
        V  = iv(sum(iv,2)~=0)-E;
        RM = inv(E'*E+V'*V)*E'*E;
        corrIV = sqrt(RM); %correlation between IV and true shock

    else

        d = eye(Nsiv)*sum(sum(iv,2)~=0)/Tiv;
        sigmm = cov(iv,1);
        sigmu = cov([iv,res2(:,1:Nsiv,s)],1);
        RM = inv(sigmm)*sigmu(1:Niv,Niv+1:end)*inv(b11b11p)*sigmu(1:Niv,Niv+1:end)'/d; %Lambda
        RMeigs = sort(eig(RM)); %if below unity ==> measurement error
        corrIV = sqrt(RMeigs); %correlation between PC of IVs and true shock

    end
    out.RM(:,:,s) = RM;
    out.corrIV(s,:) = real(corrIV);

end


