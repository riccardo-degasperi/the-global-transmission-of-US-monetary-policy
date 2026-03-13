function [ ] = IRFdisplay_all(IRFulm,longnames,chartName,plotOpt)
% Plots IRFs

% Unpack
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
folder      = plotOpt.folder;
font        = plotOpt.font;

% Determine dimensions
N = size(IRFulm,1);
h = size(IRFulm{1,1},2);

% Number of figures already plotted
hh = findobj('type','figure');
nn = length(hh);

% Open figure and title
figure(1+nn);

% Intitiate count
count = 0;

for i=1:N             %loop over variables

    for ii=1:N        %loop over shocks

    count = count+1;  %update count
    subplot(N,N,count)
    tmp = IRFulm{i,ii};

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
    plot([0,h-1],[0 0],'k--');                         %plot horizontal axis

    hold off

    % Set X-axis and Y-axis
    minband=min(tmp(1,:));        %minimum value of lower bound
    maxband=max(tmp(5,:));        %maximum value of upper bound
    space=maxband-minband;
    Ymin=minband-0.2*space-realmin;       %lower limit of the graph
    Ymax=maxband+0.2*space+realmin;       %upper limit of the graph
    set(gca,'XLim',[0 h-1],'YLim',[Ymin Ymax],'XTick',0:6:h-1,'FontName',font);

    % Create top and side titles
    if i==1
        title(longnames{1,ii},'FontWeight','normal');   %variable shocked
    end
    if ii==1
        ylabel(longnames{1,i},'FontWeight','normal');   %variable responding
    end 

    box on
    
    end
end

% top super-title
ax = axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Shock in position:','FontSize',11,'FontName',font,'FontWeight','normal');
% side super-title
ylabel('Response of:','FontSize',12,'FontName',font,'FontWeight','normal');
set(get(ax,'Ylabel'),'Visible','on');

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
cols = N;
rows = N;

    set(gcf,'PaperUnits','centimeters','PaperSize',[6.5*cols 6.5*rows]) %[x y]
    set(gcf,'PaperPosition',[0.25*cols 0.25*rows 5.75*cols 5.75*rows])  %[left bottom width height]
print(gcf,'-dpdf',[folder_,outFileName,'.pdf']);

if saveFig
    savefig(gcf,[folder_,outFileName,'.fig']);
end

end

    

