function [] = seriesExplorer(data,dataOpt,plotOpt,lags)
% seriesExplorer plots series one by one. For each series plots ACF and
% PACF and runs an Augmented-Dickey-Fueller with lags up to "lags" (to be
% specified). Then, for each series, it takes first differences and plots
% the differenced series, ACF and PACF, and runs a new ADF.
% For the series in levels, the code runs an ADF with trend and drift if
% iRW = 1, while it runs an ADF with drift only if iRW = 0.
%
% last modified: 26/03/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%
disp('    Check stationarity')

% Unpack
endo          = [data.endo data.endo_global];
iRW_endo      = [data.iRW_endo; data.iRW_global];
dates         = data.dates;
varendo       = [dataOpt.varendo dataOpt.varendo_global];
varendolong   = [dataOpt.varendolong dataOpt.varendolong_global];
units         = dataOpt.units;
folder        = plotOpt.folder;
saveCharts    = plotOpt.saveCharts;
saveFig       = plotOpt.saveFig;

for u = 1:size(units,2)
    for i = 1:size(varendo,2)

    % Select one series at the time
    tmp = endo(:,i,u);
    disp('    ---------------------------------------')
    disp(['    Series:  ',varendolong{i}])
    disp('    ---------------------------------------')

    % Plot (levels)
    figure(1)
    subplot(2,3,1)
    plot(1:size(tmp,1),tmp)
    hold on
    plot(1:size(tmp,1),zeros(size(tmp)),'--')
    hold off

        if numel(units) == 1
        title([varendolong{i},' (levels)'])
        else
        title([varendolong{i},': ',units{u},' (levels)'])
        end

        % Set X-axis and Y-axis
        minband = min(tmp);                %minimum value of lower bound
        maxband = max(tmp);                %maximum value of upper bound
        space   = maxband-minband;
        Ymin    = minband-0.2*space;        %lower limit of the graph
        Ymax    = maxband+0.2*space;        %upper limit of the graph
        hh      = size(dates,1);
        set(gca,'XLim',[1 hh],'YLim',[Ymin Ymax],'FontSize',10,'FontName','Courier New','FontWeight','normal');

        % Set x-axis labels
        ax = gca;
        xticklabels(dates(ax.XTick));

    % ACF
    subplot(2,3,2)
    autocorr(tmp)
    set(gca,'FontSize',10,'FontName','Courier New','FontWeight','normal');

    % PACF
    subplot(2,3,3)
    parcorr(tmp)
    set(gca,'FontSize',10,'FontName','Courier New','FontWeight','normal');

    % ADF & KPSS
    if     iRW_endo(i) == 0    % stationarity (drift)

        [ADFreject(i,:,u),ADFpval(i,:,u)]   = adftest(tmp,'model','ARD','lags',1:lags);

    elseif iRW_endo(i) == 1    % trend-stationarity (drift and trend)

        [ADFreject(i,:,u),ADFpval(i,:,u)]   = adftest(tmp,'model','TS','lags',1:lags);
        
        t0 = round(sqrt(size(tmp,1)));   %optimal lag length (Kwiatkowski et al. Journal of Econometrics, 1992)
        [KPSSreject(i,:,u),KPSSpval(i,:,u),reg] = kpsstest(tmp,'lags',[(t0-3):(t0+3)],'trend',true);

    end

    disp('    ADF results: (H0: unit root; 0: fail to reject H0)')
    [ADFreject(i,:,u);ADFpval(i,:,u)]
    if iRW_endo(i) == 1
    disp('    KPSS results: (H0: trend-stationarity; 0: fail to reject H0)')
    [KPSSreject(i,:,u);KPSSpval(i,:,u)]
    end

    pause

    % H0: unit root
    % H1: trend-stationarity
    % 0: fail to reject H0

    
    % First differences
    tmp2 = tmp(1:end-1) - tmp(2:end);

    % Plot (differences)
    subplot(2,3,4)
    plot(1:size(tmp2,1),tmp2)
    hold on
    plot(1:size(tmp2,1),zeros(size(tmp2)),'--')
    hold off

        if numel(units) == 1
        title([varendolong{i},' (differenced)'])
        else
        title([varendolong{i},': ',units{u},' (differenced)'])
        end
    
        title([varendolong{i},': ',units{u},' (differenced)'])

        % Set X-axis and Y-axis
        minband = min(tmp2);                %minimum value of lower bound
        maxband = max(tmp2);                %maximum value of upper bound
        space   = maxband-minband;
        Ymin    = minband-0.2*space;        %lower limit of the graph
        Ymax    = maxband+0.2*space;        %upper limit of the graph
        hh      = size(dates,1);
        set(gca,'XLim',[1 hh],'YLim',[Ymin Ymax],'FontSize',10,'FontName','Courier New','FontWeight','normal');

        % Set x-axis labels
        ax = gca;
        xticklabels(dates(ax.XTick));

    % ACF
    subplot(2,3,5)
    autocorr(tmp2)
    set(gca,'FontSize',10,'FontName','Courier New','FontWeight','normal');

    % PACF
    subplot(2,3,6)
    parcorr(tmp2)
    set(gca,'FontSize',10,'FontName','Courier New','FontWeight','normal');

    % ADF (drift)
    [ADFreject_d(i,:,u),ADFpval_d(i,:,u)] = adftest(tmp2,'model','ARD','lags',1:lags);

    disp('    ---------------------------------------')
    disp('    ADF results (differences): (H0: unit root; 0: fail to reject H0)')
    [ADFreject_d(i,:,u);ADFpval_d(i,:,u)]
    
    
    %---------------------------------------------------------------------%
    % Save plot
    if saveCharts

        folder_ = [pwd '/plots/',folder,'/seriesExplorer'];
        if ~isempty(folder)
            if not(isfolder(folder_))
                mkdir(folder_)
            end
            folder_ = [folder_,'/'];
        end

        if numel(units) == 1
            outFileName = ['ACF_',varendo{i}];
        else
            outFileName = ['ACF_',units{u},'_',varendo{i}];
        end

        set(gcf,'PaperUnits','centimeters','PaperSize',[39 26]) %[x y]
        set(gcf,'PaperPosition',[1.5 0.25 34.5 23])  %[left bottom width height]
        print(gcf,'-dpdf',[folder_,outFileName,'.pdf']);

        if saveFig
        savefig(gcf,[folder_,outFileName,'.fig']);
        end

    end
    %---------------------------------------------------------------------%

    pause
        
    end
end

close all