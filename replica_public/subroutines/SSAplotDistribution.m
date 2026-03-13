function [ ] = SSAplotDistribution(shockSequence,plotHorizon,varnames,chartName0,plotOpt)
% Plot distribution of SSA sequence of shocks vs. normal distribution

% IRFulm    : cell(N,N) containing median and upper/lower bounds of IRF;
% varlist   : short names of the N variables in IRFulm;
% longnames : long names of the N variables for subplot titles;
% shockVar  : short name of variable whose shock is normalized to 1;
% plotOpt   : structure with charts options;
% plotCase  : {'allSeries','varplot'}

% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

color0 = [0.8500 0.3250 0.0980];

% Unpacking
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
folder      = plotOpt.folder;
font        = plotOpt.font;
rows        = plotOpt.allSeries.plotRows;
cols        = plotOpt.allSeries.plotCols;

% Number of figures already plotted
nF = findobj('type','figure');   
nn = length(nF);

% Select horizon
hh = squeeze(shockSequence(plotHorizon,:,:));

% Plot
figure(nn+1)
for i = 1:size(hh,1)
    [f,xi] = ksdensity(hh(i,:),'kernel','epanechnikov');
    subplot(rows,cols,i)
    area([-4:0.1:4],normpdf([-4:0.1:4]),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
    hold on
    if abs(mean(hh(i,:)))<0.001
        xline(0,'Color',color0)
    elseif mean(hh(i,:)) > 0.99 && mean(hh(i,:)) < 1.01
        xline(1,'Color',color0)
    else
        plot(xi,f,'Color',color0);
    end
    ylim([0 2])
    xlim([-4, 4])
    title(varnames{i},'FontWeight','normal')
    set(gca,'FontName',font);
    box on
end
ax = axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Shock in position:','FontSize',11,'FontName',font,'FontWeight','normal');

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

set(gcf,'PaperUnits','centimeters','PaperSize',[6.5*cols 6.5*rows]) %[x y]
set(gcf,'PaperPosition',[0.25*cols 0.25*rows 5.75*cols 5.75*rows])  %[left bottom width height]
print(gcf,'-dpdf',[folder_,outFileName,'.pdf']);

if saveFig
savefig(gcf,[folder_,outFileName,'.fig']);
end

end


end
