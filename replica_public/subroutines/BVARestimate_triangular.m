function out = BVARestimate_triangular(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt)
% Estimates the reduced-form BVAR using a triangularisation algorithm.

% INPUT:
% - data     : output from makeXY;
% - dims     : output from makeXY;
% - hp       : hyperparameters and related options;
% - modelOpt : structure containing model options;
% - dataOpt  : structure containing data options; 
% - miscOpt  : structure containing misc options;

% OUTPUT:
% - out      : structure containing estimation results and IRFs;

% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%
disp('------------------------------------------')
disp('Estimating Bayesian VAR (triangular)')

% Error if name of routine is misspecified
if ~strcmp(miscOpt.caseTriangular,'Chan') 
error('ERROR: check miscOpt.caseTriangular in the initial options.')
end

% Unpacking
X               = data.X;
Y               = data.Y;
endo            = data.endo;
T               = dims.T;
N               = dims.N;
p               = dims.p;
m               = dims.m;
dropExplosive   = miscOpt.dropExplosive;
sim             = modelOpt.sim;
burnin          = modelOpt.burnin;
jump            = modelOpt.jump;
draws           = (sim-burnin)/jump;
shutVar         = modelOpt.shutVar;
shutEq          = modelOpt.shutEq;
GLP             = hp.GLP;
folder          = plotOpt.folder;
plotOpt.ylab    = dataOpt.ylab;
varendo         = dataOpt.varendo;
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
order  = 1;                %we fit an AR(order) for the whole sample T
hp.SS  = NaN(N,1);         %container for residual standard deviations
hp.SS2 = NaN(N,1);         %container for residual variances
hp.e2  = NaN(T-order,N);   %container for residuals


for i = 1:N

    % Build X and Y matrices
    X2  = [ones(T-order,1) lagX(endo(:,i),order)]; % (T-p)x(p+1) matrix (with constant)
    Y2  = endo(order+1:end,i);                     % (T-p)x1 matrix
    
    % Intermediate step
    XX2 = X2'*X2;
    XY2 = X2'*Y2;

    % Compute beta and std of residuals
    b2         = XX2\XY2;
    hp.e2(:,i) = Y2-X2*b2;
    hp.SS(i)   = std(Y2-X2*b2);
    hp.SS2(i)  = hp.SS(i).^2;
    
end

hpVAR = hp;   %initialise hpVAR

%-------------------------------------------------------------------------%
% Patch: reorder equations moving shutEq on top

% NOTE: this is necessary to cast priors that are not block-diagonal.

%Find location of equations to shut
lEq1 = ismember(varendo,shutEq);               %for reordering
lEq2 = false(1,N);  lEq2(1:sum(lEq1)) = true;  %to go back to previous order

%Reorder X and Y
Y = [Y(:,lEq1) Y(:,~lEq1)];

for i = 1:p
    X_tmp = X(:,(i-1)*N+1:i*N);
    X(:,(i-1)*N+1:i*N) = [X_tmp(:,lEq1) X_tmp(:,~lEq1)];
end

% Reorder position of coefficients to shrink
lCoef = ismember(varendo,shutVar);
lCoef = [lCoef(lEq1) lCoef(~lEq1)];

% Reorder other relevant items
hpVAR.iRW = [hpVAR.iRW(lEq1);hpVAR.iRW(~lEq1)];
hpVAR.SS  = [hpVAR.SS(lEq1);hpVAR.SS(~lEq1)];
hpVAR.SS2 = [hpVAR.SS2(lEq1);hpVAR.SS2(~lEq1)];
hpVAR.e2  = [hpVAR.e2(:,lEq1) hpVAR.e2(:,~lEq1)];

% Update positions for casting dogmatic prior
% hpVAR.iPol   = find(strcmp(varendo,shockVar));      %position of policy variable
hpVAR.iParam = find(lCoef);                         %variables to shrink to zero
hpVAR.iEq    = 1:sum(lEq1);                         %equations to shrink


%-------------------------------------------------------------------------%
% EQUATION-BY-EQUATION HOMOSCEDASTIC Chan

% Check
if dropExplosive
error('ERROR: dropExplosive not implemented for this case.');
end

% Optimal priors
if GLP
[ml_opt,hpVAR] = get_OptKappa_mod(Y,X,p,m,hpVAR);
end

% Containers
Beta = zeros(sim,N*(N*p+m));
Alp  = zeros(sim/jump,N*(N-1)/2);
Sig  = zeros(sim/jump,N);

count_alp = 0;

for ii = 1:N
    
    yi = Y(:,ii);
    ki = N*p+ii +m-1;
    
    [mi,Vi,nui,Si,~] = prior_ACPi_mod(N,p,m,ii,hpVAR);
    Xi = [X -Y(:,1:ii-1)];                                                      %adding -Y to the regressors (we are in structural form)
    
    % compute the parameters of the posterior distribution
    iVi = Vi\speye(ki);
    Kthetai = iVi + Xi'*Xi;
    CKthetai = chol(Kthetai,'lower');    
    thetai_hat = (CKthetai')\(CKthetai\(iVi*mi + Xi'*yi));                      %sigma is not there in Xi'yi because Sigma = eye(N) in structural form
    Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;
    
    % sample sig and theta
    Sigi = 1./gamrnd(nui+(T-p)/2,1./Si_hat,sim,1);
    U = randn(sim,ki).*repmat(sqrt(Sigi),1,ki);
    Thetai = repmat(thetai_hat',sim,1) + U/CKthetai;
    
    Sig(:,ii) = Sigi;
    Beta(:,(ii-1)*(N*p+m)+1:ii*(N*p+m)) =  Thetai(:,1:N*p+m);
    Alp(:,count_alp+1:count_alp+ii-1) =  Thetai(:,N*p+m+1:end);                  %estimates of the free elements in A
    count_alp = count_alp + ii -1;
end   

store_alp  = Alp(burnin+1:jump:end,:);
store_beta = Beta(burnin+1:jump:end,:);
store_Sig  = Sig(burnin+1:jump:end,:);

[store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);

beta_tmp = permute(store_Btilde,[2,1]);
betas = reshape(beta_tmp,N*p+m,N,draws);
sigmas = permute(store_Sigtilde,[2,3,1]);

es = NaN(T-p,N,draws);
for s = 1:draws
    es(:,:,s) = Y - X*betas(:,:,s);
end

% Check maximum eigenvalue
neig = 0;

A = NaN(N*p,N*p,draws);                                   %create companion matrix
A(1:N,:,:) = permute(betas(1:end-m,:,:),[2,1,3]);
A(N+1:end,:,:) = repmat([eye(N*(p-1)) zeros(N*(p-1),N)],1,1,draws);
for s = 1:size(A,3)
    maxeig(s,1) = max(abs(eig(A(:,:,s))));
    if maxeig(s,1) >= 1
        neig = neig+1;
    end
end


% Patch: reorder equations

% Reorder betas
betas_tmp3 = [];
for i = 1:p
    betas_tmp2 = NaN(N,N,draws);
    betas_tmp1 = betas((i-1)*N+1:i*N,:,:);
    betas_tmp2(lEq1,:,:) = betas_tmp1(lEq2,:,:);
    betas_tmp2(~lEq1,:,:) = betas_tmp1(~lEq2,:,:);
    betas_tmp3 = [betas_tmp3;betas_tmp2];
end
betas_tmp3 = [betas_tmp3;betas(N*p+1:N*p+m,:,:)];

betas0 = NaN(size(betas));
betas0(:,lEq1,:) = betas_tmp3(:,lEq2,:);
betas0(:,~lEq1,:) = betas_tmp3(:,~lEq2,:);

% Reorder residuals
es_tmp = NaN(size(es));
es_tmp(:,lEq1,:) = es(:,lEq2,:);
es_tmp(:,~lEq1,:) = es(:,~lEq2,:);
es0 = es_tmp;

% Reorder sigmas
sigmas0 = nan(size(sigmas));
sigmas0(lEq1,lEq1,:) = sigmas(lEq2,lEq2,:);
sigmas0(~lEq1,~lEq1,:) = sigmas(~lEq2,~lEq2,:);
sigmas0(lEq1,~lEq1,:) = sigmas(lEq2,~lEq2,:);
sigmas0(~lEq1,lEq1,:) = sigmas(~lEq2,lEq2,:);

% Pack
es = es0;
betas = betas0;
sigmas = sigmas0;

%-------------------------------------------------------------------------%
% Tabulate VAR coefficients
disp('--> Tabulate VAR coefficients...')

% Get confidence regions
betas_l = quantile(betas(1:end-m,:,:),0.05,3); %5% lower bound
betas_m = quantile(betas(1:end-m,:,:),0.5,3);  %median
betas_u = quantile(betas(1:end-m,:,:),0.95,3); %95% upper bound

% N.B. We are dropping the coeffs on constant and exogenous

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
fileName = [char(units0),'_Triangular_coeffs.xlsx'];
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
out.hpVAR   = hp;
out.maxeig  = maxeig;
out.neig    = neig;
out.IRF     = IRF;