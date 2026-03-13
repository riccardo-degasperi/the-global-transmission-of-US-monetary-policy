function res = handleVARtrend(Y,X,vecB,h,p)
% Removes deterministic component from Y.
% Expects vecB = [N(N*p+m)x1]

% Last modified: 28/06/2019 r.degasperi@warwick.ac.uk

%-------------------------------------------------------------------------%
% OPTIONS
plotTrend = 1;
savePlot  = 0;
%-------------------------------------------------------------------------%

nB      = numel(vecB);
[T,N]   = size(Y);
[~,Npm] = size(X);
m = Npm-N*p;

% trend        = NaN(T+h+1,Npm);
% trend(:,end) = 1;
trend          = zeros(T+h+1,Npm);
trend(1:T,N*p+1:end) = X(:,N*p+1:end);
trend(1,:)     = X(1,:);

for j = 1:T+h
    
    trend(j+1,1:N)       = trend(j,:)*reshape(vecB,nB/N,N);
    trend(j+1,N+1:end-m) = trend(j,1:end-N-m);
    
end

trend       = trend(:,1:N);
Y_detrended = Y-trend(2:T+1,:);

%-------------------------------------------------------------------------%
if plotTrend
    
    hh = findobj('type','figure');
    nn = length(hh);
    figure(nn+1)
    plotRows = ceil(N/3);
    pln = 1;
    for j = 1:N
        subplot(plotRows,ceil(N/plotRows),pln)

        plot(Y(:,j))
        hold on
        plot(trend(2:T,j),'--r')
        axis tight
        
        pln = pln+1;
    end
    
    set(gcf,'PaperUnits','centimeters','PaperSize',[18 15]) %[x y]
    set(gcf,'PaperPosition',[-1 0 20 15]) %[left bottom width height]
    
    if savePlot
    print(gcf,'-dpdf','VARtrend.pdf');
    saveas(gcf,'VARtrend.fig');
    end
    
end

res.detrended = [zeros(p,N);Y_detrended];
res.trend     = trend(2:end,:);           %includes horizons

