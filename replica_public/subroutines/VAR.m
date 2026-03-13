function out = VAR(Y,X)

% INPUT:
% - Y and X matrices


%-------------------------------------------------------------------------%
% OPTIONS
statistics = 0;
cols = 2;

%-------------------------------------------------------------------------%

[T,Np] = size(X);    %T = T-p; Np = N*p+m
[~,N]  = size(Y);

% Check degres of freedom
if Np>=T
    error('ERROR: N*p+m > T-p. Not enough degrees of freedom.')
end

XY    = X'*Y;
XX    = X'*X;

out.beta   = XX\XY;
fitted     = X*out.beta;
out.e      = Y - fitted;
out.sigma  = (out.e'*out.e)/(T-Np);
out.loglik = -(T*N/2)*log(2*pi)+T/2*log(det(inv(out.sigma)))-(T*N/2);

%-------------------------------------------------------------------------%
% STATISTICS

% Plot residuals vs. fitted values
if statistics

vars = size(Y,2);
rows = ceil(vars/2);
    
figure
    for i = 1:vars
    
    subplot(rows,cols,i)
    scatter(fitted(:,i),out.e(:,i),'x')
    hold on
    
    ylabel('Residuals') 
    xlabel('Fitted values')
    
    end

end

