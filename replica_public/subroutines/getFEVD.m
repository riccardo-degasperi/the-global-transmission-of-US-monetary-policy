function out = getFEVD(in,dims,modelOpt)
% See BEAR toolboox compendium Ch.5
% Credit: Ambrogio Cesa-Bianchi's Toolboox v3.0

% Display message
disp('--> Estimate forecast error variance decomposition...')

% Unpack
IRF             = in.IRF;             %reduced form IRFs
sigmas          = in.sigmas;
B0              = in.B0nnr;           %non normalised impact matrix
sim             = modelOpt.sim;
burnin          = modelOpt.burnin;
jump            = modelOpt.jump;
h               = dims.h+1;           %+1 to have a round h horizons in the charts
N               = dims.N;
draws           = (sim-burnin)/jump;

% Initialise containers
PSI        = zeros(N,N,h,draws);
FEVD       = zeros(h,N,N,draws);
SE         = zeros(h,N,draws);

% Compute Wold representation
for s = 1:draws
    for i = 1:N
        for ii = 1:N
            for j = 1:h
                PSI(i,ii,j,s) = IRF{i,ii}(j,s);
            end
        end
    end
end

% Compute FEVD
for s = 1:draws

    % Initialise containers for iteration s
    MSE   = zeros(N,N,h);
    MSE_j = zeros(N,N,h);

    for i = 1:N
        
        % Compute the mean square error
        MSE(:,:,1) = sigmas(:,:,s);
        for j = 2:h
           MSE(:,:,j) = MSE(:,:,j-1) + PSI(:,:,j,s)*sigmas(:,:,s)*PSI(:,:,j,s)';
        end

        % Compute the i^th structural forecast error
        MSE_j(:,:,1) = B0(:,i,s)*B0(:,i,s)';
        for j = 2:h
            MSE_j(:,:,j) = MSE_j(:,:,j-1) + PSI(:,:,j,s)*MSE_j(:,:,1)*PSI(:,:,j,s)';   
        end
    
        % Compute the Forecast Error Covariance Decomposition
        FECD = MSE_j./MSE;
    
        % Select only the variance terms
        for j = 1:h
            for ii = 1:N
                FEVD(j,i,ii,s) = 100*FECD(ii,ii,j);
                SE(j,:,s) = sqrt(diag(MSE(:,:,j))');
            end
        end
        % FEVD(j,i,ii): matrix with 'j' horizons, the FEVD due to 'i' shock for 'ii' variable

    end
end

% Median, mean, and bands
out.FEVDmean = mean(FEVD,4);
out.FEVDp95  = prctile(FEVD,95,4);
out.FEVDp84  = prctile(FEVD,84,4);
out.FEVDp50  = prctile(FEVD,50,4);
out.FEVDp16  = prctile(FEVD,16,4);
out.FEVDp5   = prctile(FEVD,5,4);

end