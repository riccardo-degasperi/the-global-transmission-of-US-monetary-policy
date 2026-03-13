function HD = getHD(in,data,dims,modelOpt)
% See BEAR toolboox compendium Ch.5
% Credit: Ambrogio Cesa-Bianchi's Toolboox v3.0

% Display message
disp('--> Estimate historical decomposition...')

% Unpack
sim             = modelOpt.sim;
burnin          = modelOpt.burnin;
jump            = modelOpt.jump;
const           = modelOpt.constant;
betas           = in.betas;
es              = in.es;
B0              = in.B0nnr;  %1sd
maxeig          = in.maxeig(burnin+1:jump:end);
N               = dims.N;
Nx              = dims.Nx;
T               = dims.T - dims.p;
m               = dims.m;
mc              = dims.mc;
p               = dims.p;
pexo            = dims.pexo;
draws           = (sim-burnin)/jump;
X               = data.X;
Nex             = Nx*pexo+1;

% Check
if max(maxeig)>1
    warning('max(eig(A)) > 1. If sum(shocks) diverge, try to use Sum of Coefficients and/or Dummy Initial Observation priors. For a quick check, set miscOpt.dropExplosive = 1. The explosive draws are likely the cause for this issue.')
end

% Initialise containers
HDendo  = zeros(N,T+1,draws);
HDshock = zeros(N,T+1,N,draws);
HDinit  = zeros(N,T+1,draws);
HDconst = zeros(N,T+1,draws);
HDtrend = zeros(N,T+1,draws);
HDexo   = zeros(N,T+1,Nex,draws);

for s = 1:draws

    % Get companion form
    A = NaN(N*p);
    A(1:N,:) = betas(1:end-m-mc,:,s)';
    A(N+1:end,:) = [eye(N*(p-1)) zeros(N*(p-1),N)];

    % Contribution of each shock
    eps = B0(:,:,s)\es(:,:,s)'; % structural errors 
    B0_big = zeros(N*p,N);
    B0_big(1:N,:) = B0(:,:,s);
    Icomp = [eye(N) zeros(N,(p-1)*N)];
    HDshock_big = zeros(N*p,T+1,N);
    for i = 1:N % for each variable
        eps_big = zeros(N,T+1); % matrix of shocks conformable with companion
        eps_big(i,2:end) = eps(i,:);
        for t = 2:T+1
            HDshock_big(:,t,i) = B0_big*eps_big(:,t) + A*HDshock_big(:,t-1,i);
            HDshock(:,t,i,s) =  Icomp*HDshock_big(:,t,i);
        end
    end

    % Initial value
    HDinit_big = zeros(N*p,T+1);
    HDinit_big(:,1) = X(1,1:N*p)';
    HDinit(:,1) = Icomp*HDinit_big(:,1);
    for t = 2:T+1
        HDinit_big(:,t) = A*HDinit_big(:,t-1);
        HDinit(:,t,s) = Icomp*HDinit_big(:,t);
    end
    
    % Constant
    HDconst_big = zeros(N*p,T+1);
    CC = zeros(N*p,1);
    if const>0
        CC(1:N,:) = betas(N*p+1,:,s);
        for t = 2:T+1
            HDconst_big(:,t) = CC + A*HDconst_big(:,t-1);
            HDconst(:,t,s) = Icomp * HDconst_big(:,t);
        end
    end
    
    % Linear trend
    HDtrend_big = zeros(p*N,T+1);
    TT = zeros(N*p,1);
    if const>1
        TT(1:N,:) = betas(N*p+2,:,s);
        for t = 2:T+1
            HDtrend_big(:,t) = TT*(t-1) + A*HDtrend_big(:,t-1);
            HDtrend(:,t,s) = Icomp * HDtrend_big(:,t);
        end
    end
    
    % Exogenous
    HDexo_big = zeros(N*p,T+1);
    EXO = zeros(N*p,Nx*(pexo+1));
    if Nx>0
        for i = 1:Nex
            VARexo = X(:,N*p+const+i);
            EXO(1:N,i) = betas(N*p+const+i,:,s);
            for t = 2:T+1
                HDexo_big(:,t) = EXO(:,i)*VARexo(t-1,:)' + A*HDexo_big(:,t-1);
                HDexo(:,t,i,s) = Icomp * HDexo_big(:,t);
            end
        end
    end
    
    % All decompositions must add up to the original data
    HDendo(:,:,s) = HDinit(:,:,s) + HDconst(:,:,s) + HDtrend(:,:,s) + sum(HDexo(:,:,:,s),3) + sum(HDshock(:,:,:,s),3);
        
end

% Get means
meanHDendo  = mean(HDendo,3);
meanHDinit  = mean(HDinit,3);
meanHDconst = mean(HDconst,3);
meanHDtrend = mean(HDtrend,3);
meanHDexo   = mean(HDexo,4);
meanHDshock = mean(HDshock,4);

HDshockp95  = prctile(HDshock,95,4);
HDshockp84  = prctile(HDshock,84,4);
HDshockp50  = prctile(HDshock,50,4);
HDshockp16  = prctile(HDshock,16,4);
HDshockp5   = prctile(HDshock,5,4);


% Save and reshape all HDs
HD.shock    = zeros(T+p,N,N);  % [nobs x shock x var]
HD.shockp95 = zeros(T+p,N,N);  % [nobs x shock x var]
HD.shockp84 = zeros(T+p,N,N);  % [nobs x shock x var]
HD.shockp50 = zeros(T+p,N,N);  % [nobs x shock x var]
HD.shockp16 = zeros(T+p,N,N);  % [nobs x shock x var]
HD.shockp5  = zeros(T+p,N,N);  % [nobs x shock x var]
    for i = 1:N
        for ii = 1:N
            HD.shock(:,ii,i)    = [NaN(p,1); meanHDshock(i,2:end,ii)'];
            HD.shockp95(:,ii,i) = [NaN(p,1); HDshockp95(i,2:end,ii)'];
            HD.shockp84(:,ii,i) = [NaN(p,1); HDshockp84(i,2:end,ii)'];
            HD.shockp50(:,ii,i) = [NaN(p,1); HDshockp50(i,2:end,ii)'];
            HD.shockp16(:,ii,i) = [NaN(p,1); HDshockp16(i,2:end,ii)'];
            HD.shockp5(:,ii,i)  = [NaN(p,1); HDshockp5(i,2:end,ii)'];
        end
    end
HD.init  = [nan(p-1,N); meanHDinit(:,1:end)'];    % [nobs x var]
HD.const = [nan(p,N);   meanHDconst(:,2:end)'];   % [nobs x var]
HD.trend = [nan(p,N);   meanHDtrend(:,2:end)'];   % [nobs x var]
HD.exo   = zeros(T+p,N,Nx);                       % [nobs x var x var_ex]
    for i = 1:Nx
        HD.exo(:,:,i) = [nan(p,N);   meanHDexo(:,2:end,i)'];
    end
HD.endo  = [nan(p,N);   meanHDendo(:,2:end)'];    % [nobs x var]



end