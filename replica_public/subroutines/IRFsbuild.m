function IRFs = IRFsbuild(betas,IRF,B0,dims,shockSign)
% IRFsbuild constructs structural IRFs starting from the reduced-form ones.

% INPUT:
% betas : (Np+m) x N x sim matrix of coefficients
% IRF   : reduced-form IRFs generated with IRFbuild.m
% B0    : N x N x sim      impact matrix

%-------------------------------------------------------------------------%

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

switch IRFtype

    case 'hamilton'

    % Retrieve dimensions
    [tmp,N,sim] = size(betas);
    p = (tmp-m-mc)/N;

    IRFs = cell(N,N);
    Yh  = NaN(h+1,N*p);

    for s = 1:sim     %loop over simulations

        for i = 1:N   %loop over shocks
        Yh(1,:) = zeros(1,N*p);
        Yh(1,1:N) = sign*B0(:,i,s)';

            if h == 0

                for ii = 1:N
                IRFs{ii,i}(:,s) = Yh(:,ii);
                end

            else

                for j = 1:h       %loop over horizons
                Yh(j+1,:) = [Yh(j,:)*betas(1:end-m-mc,:,s) Yh(j,1:N*(p-1))];
                end

                for ii = 1:N
                IRFs{ii,i}(:,s) = Yh(:,ii);
                end

            end

        end
    end


    case 'VMA'


    % Retrieve dimensions
    [~,N,sim] = size(betas);

    % Initialise containers for IRFs and 
    sigmas  = NaN(N,N);
    IRFs = cell(N,N);

    for s = 1:sim      %loop over simulations
        for j = 1:h+1  %loop over horisons

            % Get the NxN psi matrix
            for i = 1:N
                for ii = 1:N
                    sigmas(i,ii) = IRF{i,ii}(j,s); 
                end
            end

            % Obtain the NxN psitilda
            psitilde = sigmas*B0(:,:,s);

            % Put it back into the cell
            for i=1:N
                for ii=1:N  
                    IRFs{i,ii}(j,s) = psitilde(i,ii);
                end
            end

        end %horizons
    end %simulations

end

