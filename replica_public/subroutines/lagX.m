function out = lagX(X,p)
% Generates matrix of lagged values
%-------------------------------------------------------------------------%
% Inputs:                                                                 %
% X   : TxN matrix of data                                                %
% p   : number of lags                                                    %
%                                                                         %
% Outputs:                                                                %
% out : (T-p)x(Np) matrix of lagged values                                %
%                                                                         %
% NOTE: the correct matrix Y(t) of contemporaneous values is X(p+1:end,:) %
%                                                                         %
% r.degasperi@warwick.ac.uk                                               %
%-------------------------------------------------------------------------%

    [T,N] = size(X);

    out = NaN(T-p,N*p);

    for i = 1:p
        out(:,(i-1)*N+1:i*N) = X(p+1-i:end-i,:);
    end

end

