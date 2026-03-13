function [bb,MM] = compactForecastNotation(B0,betas,data,dims)
% Rewrite SVAR forecast from autoregressive form to moving average form.
% Adapted from Antolin-Diaz, Petrella, Rubio-Ramirez (2021).
% See Online Appendix, Section A
%
% This only works for modelOpt.constant = 1.
%
% Last modified: 06/03/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpack
m         = dims.m;
mc        = dims.mc;
N         = dims.N;
p         = dims.p;
h         = dims.h+1;             %+1 for consistency with IRFbuild function
T         = dims.T - dims.p;
Y         = data.Y;
D         = data.Xcovid;
draws     = size(betas,3);

% Loop over draws
for s = 1:draws

    B = betas(1:end-m-mc,:,s);
    c = betas(end-m-mc+1:end-mc,:,s);    %constant
    d = betas(end-mc+1:end,:,s);         %parameters on dummies

    % Reshape betas
    B = permute(reshape(B',N,N,p),[2,1,3]);
    
    %---------------------------------------------------------------------%
    % Create K matrices
    K_h = nan(N,N,h);
    K_h(:,:,1) = eye(N);
    
    for i = 1:h-1
       sumterm = zeros(N,N);
       for j = 1:i
           if j <= p
               sumterm = sumterm + K_h(:,:,i+1-j)*B(:,:,j);
           end
       end
       K_h(:,:,i+1) = eye(N) + sumterm;
    end

    %---------------------------------------------------------------------%
    % Create Delta matrices (dummies are zero for all h>T)
    D_h = nan(1,N,h);
    Dhat = zeros(1,mc);     %here you can change the forward path of the dummy
    D_h(:,:,1) = Dhat*d;

    for i = 2:h
        sumterm = Dhat*d;
        for j = 1:i-1
            if j <= p
                sumterm = sumterm + D_h(:,:,i-j)*B(:,:,j);
            end
        end
        D_h(:,:,i) = sumterm;
    end

    %---------------------------------------------------------------------%
    % Create N matrices
    N_h_l = nan(N,N,h,p);

    for i = 1
        for l = 1:p
            N_h_l(:,:,1,l) = B(:,:,l);
        end
    end
    
    for i = 2:h
       for l = 1:p
           sumterm = zeros(N,N);
           for j = 1:i-1
               if j <= p
               sumterm = sumterm + N_h_l(:,:,i-j,l)*B(:,:,j);
               end
           end
           % RD: Correct. E.g. B2 doesn't exist for a VAR(1).
           if i+l-1 > p
               N_h_l(:,:,i,l) = sumterm;
           else
               N_h_l(:,:,i,l) = B(:,:,i+l-1) + sumterm;
           end
       end
    end
    
    %---------------------------------------------------------------------%
    % Create b matrix
    b = [];
    
    for h = 1:h
        
        sumterm = zeros(1,N);
        
        for l = 1:p
           sumterm  =  sumterm + Y(T+1-l,:)*N_h_l(:,:,h,l);
        end
        
        b_tplush = c*K_h(:,:,h) + D_h(:,:,h) + sumterm;
            
        b = [b b_tplush];
    
    end
    
    b = b';
    
    %---------------------------------------------------------------------%
    % Create M matrix
    M_i = nan(N,N,h);
    
    %M_i(:,:,1) = inv2(A0);
    M_i(:,:,1) = B0(:,:,s)';  %RD: impact matrix is already inverse!
    M = zeros(N*h);
    
    for i = 1:h-1
       sumterm = zeros(N,N);
       for j = 1:i
           if j <= p
               sumterm = sumterm + M_i(:,:,i+1-j)*B(:,:,j);
           end
       end
       M_i(:,:,i+1) = sumterm;
       
    end
    
    
    Mi_reshaped = reshape(M_i,N,N*h);
    
    for i = 1:h
     M(N*(i-1)+1:N*(i-1)+N,N*(i-1)+1:end) = Mi_reshaped(:,1:N*h-N*(i-1));
    end

    %---------------------------------------------------------------------%
    % Store
    bb(:,s) = b;
    MM(:,:,s) = M;

end
