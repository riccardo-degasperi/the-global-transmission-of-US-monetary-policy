function [ ] = IRFdisplay_across(IRFulm,varacross,units,ylab,chartName,plotOpt)

% Unpacking
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
folder      = plotOpt.folder;
horz        = plotOpt.frequency;
font        = plotOpt.font;

% Determine dimensions
U = size(units,1);
h = size(IRFulm{1,1},2);

% Number of figures already plotted
hh = findobj('type','figure');
nn = length(hh);

% r1 = ceil(N/cols);
% r2 = rows0;
% rows = min([r1,r2]);              %number of rows
% Nfig = ceil(N/(rows*cols));

% rows = ceil(sqrt(U));
% cols = ceil(U/rows);

rows = 3;
cols = 5;


irf = figure(1+nn);                 %open figure
set(irf,'name',[chartName,'_',varacross]);
count = 0;

for i = 1:U

    count = count + 1;

    subplot(rows,cols,count)

    tmp = IRFulm{i,1};

    hold on

    % Plot confidence bands
    Xpatch1=[(0:1:h-1) (h-1:-1:0)];
    Ypatch1=[tmp(1,:) fliplr(tmp(5,:))];
    IRFpatch1=patch(Xpatch1,Ypatch1,[.75 .75 .75]);  %build confidence bands 90%
    set(IRFpatch1,'facealpha',.5);                   %add transparency
    set(IRFpatch1,'edgecolor','none');               %remove edge color

    Xpatch2=[(0:1:h-1) (h-1:-1:0)];
    Ypatch2=[tmp(2,:) fliplr(tmp(4,:))];
    IRFpatch2=patch(Xpatch2,Ypatch2,[.65 .65 .65]);     %build confidence bands 68%
    set(IRFpatch2,'facealpha',.5);                   %add transparency
    set(IRFpatch2,'edgecolor','none');               %remove edge color

    % Plot median and X-axis
    plot(0:h-1,tmp(3,:),'Color',[.2 .4 .6],'LineWidth',2); %plot median
    plot([0 h],[0 0],'k--');                         %plot horizontal axis

    hold off

    % Set X-axis and Y-axis
    minband=min(tmp(1,:));        %minimum value of lower bound
    maxband=max(tmp(5,:));        %maximum value of upper bound
    space=maxband-minband;
    Ymin=minband-0.2*space-realmin;       %lower limit of the graph
    Ymax=maxband+0.2*space+realmin;       %upper limit of the graph
    set(gca,'XLim',[0 h-1],'YLim',[Ymin Ymax],'XTick',0:6:h-1,'FontName',font);

    ax = gca;
    ax.YAxis.Exponent = 0;
    xtickangle(0)
            
    % Create titles
    title(units{i},'FontWeight','bold');  %unit
    ylabel(ylab)
    
    % Create overall legend
%         if i == (t-1)*rows*cols+1
%             hL = legend('90%','68%','median');
%             position = [0.9 0.1 0.03 0.03];
%             units = 'normalized';
%             set(hL,'Position', position,'Units', units);
%         end

    box on

    if i == U
        xlabel(['Horizon (',horz,')'])
        break
    end

end

% Top super-title
    ax = axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title(['Response of: ',varacross],'FontSize',11,'FontName',font,'FontWeight','normal');

% Save plot
if saveCharts
    
folder_ = [pwd '/plots/',folder];
if ~isempty(folder)
    if not(isfolder(folder_))
        mkdir(folder_)
    end
    folder_ = [folder_,'/'];
end

outFileName = chartName;

    set(gcf,'PaperUnits','centimeters','PaperSize',[6.5*cols 6.5*rows]) %[x y]
    set(gcf,'PaperPosition',[0.25*cols 0.25*rows 5.75*cols 5.75*rows])  %[left bottom width height]
print(gcf,'-dpdf',[folder_,outFileName,'_',varacross,'.pdf']);

if saveFig
    savefig(gcf,[folder_,outFileName,'_',varacross,'.fig']);
end

end
    
