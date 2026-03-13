function out = BVARidentification(in,data,dims,modelOpt,dataOpt,miscOpt)
% Identifies a Bayesian VAR

% INPUT:
% - in       : output from BLPestimate.m
% - data     : output from makeXY;
% - dims     : output from makeXY;
% - modelOpt : structure containing model options;
% - dataOpt  : structure containing data options;
% - miscOpt  : structure containing misc options; 

% OUTPUT:
% - out      : structure containing estimation results and IRFs;


% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%

% Options
precision = 1.0000e-05;

%-------------------------------------------------------------------------%
% Unpacking
betas           = in.betas;
sigmas          = in.sigmas;
es              = in.es;
hpVAR           = in.hpVAR;
maxeig          = in.maxeig;
neig            = in.neig;
IRF             = in.IRF;
dates           = data.dates;
iFD             = logical(data.iFD);
X               = data.X;
Y               = data.Y;
T               = dims.T;
N               = dims.N;
p               = dims.p;
m               = dims.m;
mc              = dims.mc;
h               = dims.h;
firstDifference = miscOpt.firstDifference;
sim             = modelOpt.sim;
burnin          = modelOpt.burnin;
jump            = modelOpt.jump;
draws           = (sim-burnin)/jump;
identification  = modelOpt.identification;
varendo         = dataOpt.varendo;
shockVar        = dataOpt.shockVar;
shockSize       = dataOpt.shockSize;
shockSign       = dataOpt.shockSign;
shockIV         = dataOpt.shockIV;

% Check whether the shock is normalised or one SD
check1sd = strcmp('1sd',shockSize);


%-------------------------------------------------------------------------%
% Structural identification
disp('--> Recovering structural form...')

% Estimate impact matrix
switch identification
    
    case 'cholesky' %-----------------------------------------------------%
        
        B0nnr = NaN(N,N,draws);
        B0    = NaN(N,N,draws);
        
        for s = 1:draws
            
            % Cholesky decomposition
            B0nnr(:,:,s) = chol(nspd(sigmas(:,:,s)),'lower');

            % Unit normalisation
            B0(:,:,s) = B0nnr(:,:,s)*diag(1./diag(B0nnr(:,:,s)));
            
        end
        
        
    case 'iv' %-----------------------------------------------------------%
        
        % Initialise containers
        B0nnr = NaN(N,N,draws);
        
        % Match length of residual series to IV
        [es2,iv2,~,~,~,~] = ivMatch(p,dates,ivDates,es,iv);

        % Estimate identified columns of B0
        [B0_nr_part,B0_nnr_part,~,res] = ivB0(es2,iv2,sigmas,dims,dataOpt);

        % First-stage F-stat: median and bands
        Fdistr = [quantile(res.Fstat,0.05)' quantile(res.Fstat,0.5)' quantile(res.Fstat,0.95)'];
        Frdistr = [quantile(res.Frobust,0.05)' quantile(res.Frobust,0.5)' quantile(res.Frobust,0.95)'];

        % Correlation between (principal components of) proxies and shocks (see Mertens & Ravn, 2013, p.8)
        RMdistr = [quantile(res.corrIV,0.05)' quantile(res.corrIV,0.5)' quantile(res.corrIV,0.95)']; 

        % Get frequentist 1st stage F statistic
        if T-p > N*p+m+mc
            resid = Y-X*(X\Y);
            sigmaf = resid'*resid/(T-p-N*p-m-mc);
            [resid1s,iv1s,~,~,~,~] = ivMatch(p,dates,ivDates,resid,iv);
            [~,~,~,res_freq] = ivB0(resid1s,iv1s,sigmaf,dims,dataOpt);
            Fstat = res_freq.Fstat;
            Frobust = res_freq.Frobust;
        else
            Fstat = 'NA';
            Frobust = 'NA';
            warning('Frequentist F-stat cannot be computed because T-p < N*p+m')
        end
        
        % Obtain remaining columns of B0 (without structural meaning)
        % (Venditti & Veronese, 2021)

        % Specify number of shocks
        Ns = numel(shockVar);
        
        for s = 1:draws

            %1) Reorder policy variables first
            sigma = sigmas(:,:,s);
            iP = ismember(varendo,shockVar); %position of policy variables
            b1 = [B0_nnr_part(iP,iP,s);B0_nnr_part(~iP,iP,s)];
            Sig = [sigma(iP,iP) sigma(iP,~iP); sigma(~iP,iP) sigma(~iP,~iP)];
            
            %2) Compute Chol(Sigma) and omega1hat:
            Sigtr = chol(0.5*(Sig + Sig'))';
            Omega1hat = Sigtr\b1;
            % NOTE: We use the non-normalized matrix of impacts

            %3) Draw Omega (nxn):
            x = randn(N);
            Omega_draw = getqr(x);

            %4) Replace omega1hat in Omega:
            Omega_draw(:,1:Ns) = Omega1hat;

            %5) Gram-Schmidt orthogonalization
            Omegastar = gramschmidt(Omega_draw);
            % NOTE: if the additional restrictions used to distinguish the
            % shocks when you identify more than one shock with more than
            % one proxy are such that the identified shocks are not
            % orthogonal, this procedure will ortogonalise them,
            % potentially altering the impact coefficients on the shocks
            % ordered after the first one. If you use Choleski, the
            % identified columns will remain the same.
            if Ns>1 && sum(sum(Omega1hat(:,2:Ns)-Omegastar(:,2:Ns))) > precision
                error('ERROR: Gram-Schmidt is changing the proxy-identified columns.')
            end

            %6) Compute impact matrix B
            tmpB0nnr = Sigtr*Omegastar;

            %7) Check that the covariance matrix is maintained
            Sig_check = (Sigtr*Omegastar)*(Sigtr*Omegastar)';
            test = abs(Sig-Sig_check) > precision;
            if sum(sum(test))>0
                error('ERROR: The orthogonalized covariance differs from the input covariance.')
            end

            %8) Reorder columns
            B0nnr(iP,iP,s)   = tmpB0nnr(1:Ns,1:Ns);
            B0nnr(~iP,iP,s)  = tmpB0nnr(Ns+1:end,1:Ns);
            B0nnr(iP,~iP,s)  = tmpB0nnr(1:Ns,Ns+1:end);
            B0nnr(~iP,~iP,s) = tmpB0nnr(Ns+1:end,Ns+1:end);

        end
        
        % Store normalised partial impact matrix for later use
        B0 = B0_nr_part;

end


%-------------------------------------------------------------------------%
% Compute shock series
shock = zeros(T-p,numel(shockVar),draws);
iP = ismember(varendo,shockVar);
for s = 1:draws

    if check1sd

        % One SD normalisation
        shock(:,:,s) = (B0nnr(:,iP,s)'/sigmas(:,:,s)*es(:,:,s)')'; 

    else

        % Unit normalisation
        shock(:,:,s) = (B0(:,iP,s)'/sigmas(:,:,s)*es(:,:,s)')'/(B0(:,iP,s)'/sigmas(:,:,s)*B0(:,iP,s));
    end
end


%-------------------------------------------------------------------------%
% Generate Structural-form Impulse Response Functions
disp('--> Generating structural-form IRFs...')

% If shockSize = '1sd', use non-normalised impact matrix to generate IRFs
if check1sd
    B0_ = B0nnr;
else
    B0_ = B0;
end

% Generate structural impulse responses
IRFs = IRFsbuild(betas,IRF,B0_,dims,shockSign);


%-------------------------------------------------------------------------%
% Get IRFs median, upper and lower bounds
IRFsulm = IRFbands(IRFs,N,h);


%-------------------------------------------------------------------------%
% Pack output
out.betas        = betas;
out.sigmas       = sigmas;
out.es           = es;
out.hpVAR        = hpVAR;
out.maxeig       = maxeig;
out.neig         = neig;
out.IRF          = IRF;
out.B0nnr        = B0nnr;
out.B0           = B0;
out.IRFs         = IRFs;
out.IRFsulm      = IRFsulm;
out.shock        = shock;
if strcmp(identification,'iv')
    out.Fstat    = Fstat;
    out.Frobust  = Frobust;
    out.Fdistr   = Fdistr;
    out.Frdistr  = Frdistr;
    out.RMdistr  = RMdistr;
end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% QR decomposition
function out = getqr(a)
    [q,r] = qr(a,0);
    out = q*diag(sign(diag(r)));

% Alternative version
% % (See Rubio-Ramirez, Waggoner & Zha, 2010)
% for i = 1:size(q,1)
%     if r(i,i) < 0
%         q(:,i) = -q(:,i);
%     end
% end
% out = q;

end
