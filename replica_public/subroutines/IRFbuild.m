function IRF = IRFbuild(betas,dims,shockSign)

h  = dims.h;
m  = dims.m;
mc = dims.mc;

IRFtype = 'hamilton';

if strcmp(shockSign,'positive')
    sign = 1;
elseif strcmp(shockSign,'negative')
    sign = -1;                        %set to -1 to flip responses to negative shocks
else
    error('ERROR: shockSign can be either "negative" or "positive".')
end

% Retrieve dimensions
[tmp,N,sim] = size(betas);
p = (tmp-m-mc)/N;

switch IRFtype

    case 'hamilton'

    IRF = cell(N,N);
    Yh  = NaN(h+1,N*p);

    for s = 1:sim  %loop over simulations
        
        for i = 1:N      %loop over shocks
        Yh(1,:) = zeros(1,N*p);
        Yh(1,i) = sign;

            for j = 1:h  %loop over horizons
            Yh(j+1,:) = [Yh(j,:)*betas(1:end-m-mc,:,s) Yh(j,1:N*(p-1))];
            end

            for ii = 1:N
            IRF{ii,i}(:,s) = Yh(:,ii);
            end
        end
    end

    case 'companion'
        
    % Construct selector
    J = [eye(N) zeros(N,N*(p-1))];

    % Initialise companion form
    A            = NaN(N*p);
    A(N+1:end,:) = [eye(N*(p-1)) zeros(N*(p-1),N)];
    
    % Initialise IRF container
    IRF_tmp = NaN(N,N,h,sim);
    
    for s = 1:sim

        % Complete companion form
        A(1:N,:) = betas(1:end-m-mc,:,s)'; %drop coeffs on constant and exogenous

        % Build IRF
        for j = 1:h+1  %loop over horizons
            IRF_tmp(:,:,j,s) = sign.*J*A^(j-1)*J';      
        end
    end

    % Make it compatible to case 'hamilton'
    IRF = cell(N,N);
    for i = 1:N
        for ii = 1:N
            for j = 1:h+1
                for s = 1:sim
                IRF{i,ii}(:,s) = IRF_tmp(i,ii,:,s);
                end
            end
        end
    end

end %switch
