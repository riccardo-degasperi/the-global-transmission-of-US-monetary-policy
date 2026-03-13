function [] = plotSeries(data,dataOpt,plotOpt)
% plotSeries plots input time series one by one.
%
% last modified: 26/03/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%
disp('    Plot time series')

% Unpacking
endo        = data.endo;
exo         = data.exo;
dates       = data.dates;
varendo     = dataOpt.varendo;
varendolong = dataOpt.varendolong;
varexo      = dataOpt.varexo;
varexolong  = dataOpt.varexolong;
units       = dataOpt.units;
folder      = plotOpt.folder;
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
font        = plotOpt.font;


% Plot endogenous series
Nfig   = size(units,2);
Nplots = size(varendo,2);
hh     = size(dates,1);
for t = 1:Nfig
figure(t)
    for i = 1:Nplots
        
        hold on
        subplot(Nplots,1,i)
        
        tmp = endo(:,i,t);
        
        plot(tmp,'LineWidth',2)
        
        %xticklabels(dates)
        title(varendolong{i})
        
        % Set X-axis and Y-axis
        minband = min(tmp);                %minimum value of lower bound
        maxband = max(tmp);                %maximum value of upper bound
        space   = maxband-minband;
        Ymin    = minband-0.2*space;       %lower limit of the graph
        Ymax    = maxband+0.2*space;       %upper limit of the graph
        set(gca,'XLim',[1 hh],'YLim',[Ymin Ymax],'FontSize',10,'FontName',font,'FontWeight','normal');
        
        % Set x-axis labels
        ax = gca;
        xticklabels(dates(ax.XTick));
        
    end

    % top super-title
    ax = axes('Units','Normal','Position',[0.10 .075 .85 .88],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title([units{t}],'FontSize',11,'FontName',font,'FontWeight','normal');

    %---------------------------------------------------------------------%
    if saveCharts

        folder_ = [pwd '/plots/',folder,'/seriesExplorer'];
        if ~isempty(folder)
            if not(isfolder(folder_))
                mkdir(folder_)
            end
            folder_ = [folder_,'/'];
        end

        if numel(units) == 1
            outFileName = 'endoSeries';
        else
            outFileName = ['endoSeries_',units{t}];
        end

        cols = 1;
        rows = Nplots;

        set(gcf,'PaperUnits','centimeters','PaperSize',[39*cols 6.5*rows]) %[x y]
        set(gcf,'PaperPosition',[1.5*cols 0.25*rows 34.5*cols 5.75*rows])  %[left bottom width height]
        print(gcf,'-dpdf',[folder_,outFileName,'.pdf']);

        if saveFig
        savefig(gcf,[folder_,outFileName,'.fig']);
        end

    end
    %---------------------------------------------------------------------%
end

% Plot exogenous series
if ~isempty(exo)
    Nplots = size(varexo,2);
    figure(Nfig+1)
    for i = 1:Nplots

        hold on
        subplot(Nplots,1,i)

        plot(exo(:,i),'LineWidth',2)
        xticklabels(dates)
        title(varexolong{i},'FontSize',10,'FontName',font,'FontWeight','normal')

    end

    %---------------------------------------------------------------------%
    if saveCharts

        folder_ = [pwd '/plots/',folder,'/seriesExplorer'];
        if ~isempty(folder)
            if not(isfolder(folder_))
                mkdir(folder_)
            end
            folder_ = [folder_,'/'];
        end

        outFileName = 'exoSeries';

        cols = 1;
        rows = Nplots;

        set(gcf,'PaperUnits','centimeters','PaperSize',[39*cols 6.5*rows]) %[x y]
        set(gcf,'PaperPosition',[1.5*cols 0.25*rows 34.5*cols 5.75*rows])  %[left bottom width height]
        print(gcf,'-dpdf',[folder_,outFileName,'.pdf']);

        if saveFig
        savefig(gcf,[folder_,outFileName,'.fig']);
        end

    end 
    %---------------------------------------------------------------------%

end
