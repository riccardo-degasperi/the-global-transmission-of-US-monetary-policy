function [ ] = SSAplotMultipleScenarios(Y,dates,varnames,chartName0,customLabels,plotOpt)
% Plot multiple SSA (time-series) scenarios

% INPUT:
% - Y : structure with many fields: cf (T x N x 5 containing the realised
%       data and conditional forecast) for various scenarios and bl (T x N
%       x 5 containing realised data and unconditional forecast). The 5
%       columns in dim = 3 are the median and credible regions.


% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpacking
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
folder      = plotOpt.folder;
cb90        = plotOpt.cb90;
font        = plotOpt.font;

% Colors & Styles
style = {'-','--',':'};
cbands = [1 0 0; 0 0 1; 0 1 0];

% Determine dimensions
names = fieldnames(Y);
Nres = numel(names);
[T,N,~] = size(Y.(names{1}));
rows = N;
cols = 1;

% Number of figures already plotted
nF = findobj('type','figure');   
nn = length(nF);

% Get elements
dates_ = datenum(dates);
bl0 = Y.(names{Nres});

% Plot
figure(nn+1)
for i = 1:N

    subplot(rows,cols,i)

    bl = squeeze(bl0(:,i,:));

    hold on
    go = [];    %initialize graphic object container for legend

    % Plot confidence bands (BL)
    if cb90
    Xpatch1=[dates_' fliplr(dates_')];
    Ypatch1=[bl(:,1)' fliplr(bl(:,5)')];
    IRFpatch1=patch(Xpatch1,Ypatch1,[.75 .75 .75]);  %build confidence bands 90%
    set(IRFpatch1,'facealpha',.5);                   %add transparency
    set(IRFpatch1,'edgecolor','none');               %remove edge color
    else
    Xpatch2=[dates_' fliplr(dates_')];
    Ypatch2=[bl(:,2)' fliplr(bl(:,4)')];
    IRFpatch2=patch(Xpatch2,Ypatch2,[.65 .65 .65]);  %build confidence bands 68%
    set(IRFpatch2,'facealpha',.5);                   %add transparency
    set(IRFpatch2,'edgecolor','none');               %remove edge color
    end

    for ii = 1:Nres-1

        cf0 = Y.(names{ii});
        cf = squeeze(cf0(:,i,:));
    
        % Plot confidence bands (CF)
        if cb90
        Xpatch1=[dates_' fliplr(dates_')];
        Ypatch1=[cf(:,1)' fliplr(cf(:,5)')];
        IRFpatch1=patch(Xpatch1,Ypatch1,cbands(ii,:));  %build confidence bands 90%
        set(IRFpatch1,'facealpha',0);                   %add transparency
        set(IRFpatch1,'edgecolor',cbands(ii,:));        %remove edge color
        set(IRFpatch1,'LineStyle',style{ii});           %set line style
        else
        Xpatch2=[dates_' fliplr(dates_')];
        Ypatch2=[cf(:,2)' fliplr(cf(:,4)')];
        IRFpatch2=patch(Xpatch2,Ypatch2,cbands(ii,:));  %build confidence bands 68%
        set(IRFpatch2,'facealpha',0);                   %add transparency
        set(IRFpatch2,'edgecolor',cbands(ii,:));        %remove edge color
        set(IRFpatch2,'LineStyle',style{ii});           %set line style
        end
    
        % Plot median (CF)
        go(ii) = plot(dates_,cf(:,3),'Color',cbands(ii,:),'LineStyle',style{ii},'LineWidth',2); %plot median

    end

    % Plot median (BL)
    go(Nres) = plot(dates_,bl(:,3),'Color','k','LineWidth',2); %plot median

    % Plot x-axis
    plot([dates_(1) dates_(end)],[0 0],'k--');                         %plot horizontal axis

    hold off

    % Set X-axis and Y-axis
    if cb90
    minband=min([min(bl(:,1)) min(cf(:,1))]);        %minimum value of lower bound
    maxband=max([max(bl(:,5)) max(cf(:,5))]);        %maximum value of upper bound
    else
    minband=min([min(bl(:,2)) min(cf(:,2))]);        %minimum value of lower bound
    maxband=max([max(bl(:,4)) max(cf(:,4))]);        %maximum value of upper bound
    end
    space=maxband-minband;
    Ymin=minband-0.2*space-realmin;       %lower limit of the graph
    Ymax=maxband+0.2*space+realmin;       %upper limit of the graph
    set(gca,'XLim',[dates_(1) dates_(end)],'YLim',[Ymin Ymax],'FontName',font);
    datetick('x',12,'keeplimits');
    
    ax = gca;
    ax.YAxis.Exponent = 0;
    xtickangle(0)
            
    % Create titles
    title(varnames{1,i},'FontWeight','normal');  %variable shocked
    box on

    % Create overall legend
    if i == 1
        hL = legend(go,customLabels,'Location','northwest');
        units = 'normalized';
        set(hL,'Units',units);
    end
        
end


% Save plot
if saveCharts

    folder_ = [pwd '/plots/',folder];
    if ~isempty(folder)
        if not(isfolder(folder_))
            mkdir(folder_)
        end
        folder_ = [folder_,'/'];
    end
    
    outFileName = chartName0;
    
    set(gcf,'PaperUnits','centimeters','PaperSize',[39*cols 6.5*rows]) %[x y]
    set(gcf,'PaperPosition',[1.5*cols 0.25*rows 34.5*cols 5.75*rows])  %[left bottom width height]
    print(gcf,'-dpdf',[folder_,outFileName,'.pdf']);
    
    if saveFig
    savefig(gcf,[folder_,outFileName,'.fig']);
    end

end


end
