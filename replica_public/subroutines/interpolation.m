function data = interpolation(data,interpolCase,ws)

% Find series with NaNs
[~,iCol] = find(isnan(data));
iCol = unique(iCol);
X = data(:,iCol);
[T,N] = size(X);
iNaN = isnan(X);

switch interpolCase

    case 'MA'

    % Replace NaNs with centered MA predictions
    for k = 1:N
        x = X(:,k);
        x(iNaN(:,k)) = nanmedian(X(:,k),1);
        x_MA = filter(ones(2*ws+1,1)/(2*ws+1),1,[x(1)*ones(ws,1);x;x(end)*ones(ws,1)]);
        x_MA = x_MA(2*ws+1:end);
        x(iNaN(:,k)) = x_MA(iNaN(:,k));
        X(:,k) = x;
    end

    % Replace in data matrices
    data(:,iCol) = X;

    case 'spline'

    % Replace NaNs with cubic spline predictions
    for k = 1:N  
        x = X(:,k);
        x = spline(find(~iNaN(:,k)),x(~iNaN(:,k)),1:T)';
        X(:,k) = x;
    end

    % Replace in data matrices
    data(:,iCol) = X;

end %switch


