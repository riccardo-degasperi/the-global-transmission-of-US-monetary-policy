function [ ] = IRFdisplay_homogeneity(IRFstruc,dataOpt,plotOpt,legNames,chartName)
% Plots IRFs



% Unpacking
varlist     = dataOpt.varendo;
longnames   = dataOpt.varendolong;
shockVar    = dataOpt.shockVar;
plotCase    = 'allSeries';

saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
rows0       = plotOpt.(plotCase).plotRows;
cols        = plotOpt.(plotCase).plotCols;
folder      = plotOpt.folder;
horz        = plotOpt.frequency;
ylab        = plotOpt.ylab;


style = [0.3010, 0.7450, 0.9330;
        0, 0.4470, 0.7410;
        0, 0.75, 0.75;
        0, 0, 1;
        0.8500, 0.3250, 0.0980;
        0.6350, 0.0780, 0.1840;
        1, 0, 0;
        1, 1, 0;
        0.9290, 0.6940, 0.1250;
        0.75, 0.75, 0;
        0.4940, 0.1840, 0.5560;
        0.75, 0, 0.75;
        1.0000, 0.0745, 0.6510;
        0.4660, 0.6740, 0.1880;
        0, 0.5, 0;];


% Determine dimensions
names = fieldnames(IRFstruc);
N = size(varlist,2);
h = size(IRFstruc.(names{1}){1,1},2);

% Legend names (if not defined)
if isempty(legNames)
    names0 = strrep(names,'_',', ');
    legNames = [cellstr(char(names0(1:end-1)));names0(end)];
end

% Number of figures already plotted
hh = findobj('type','figure');
nn = length(hh);

r1 = ceil(N/cols);
r2 = rows0;
rows = min([r1,r2]);              %number of rows
Nfig = ceil(N/(rows*cols));


% Modify default groot styles and colours
% set(groot,...
%     'defaultAxesColorOrder',[0 0 0],...
%     'defaultAxesLineStyleOrder','-|--|:|-.|:o|')

% Loop over figures
for t = 1:Nfig

    irf = figure(t+nn);      %open figure
    set(irf,'name',chartName);
    count = 0;
    iShock = strcmp(varlist,shockVar); %position of shock variable
    ShockName = longnames(iShock);     %name of the shock variable

    for i = (t-1)*rows*cols+1:t*rows*cols
        
        count = count + 1;
        tmp_margins = [];              %collate series to then pick minband and maxband

        subplot(rows,cols,count)

        hold on
        
        for j = 1:numel(names)
        
            tmp = IRFstruc.(names{j}){i,iShock};
            tmp_margins = [tmp_margins tmp];

            % Plot median and X-axis
            if strcmp(names{j},'median')
                
                Xpatch1=[(0:1:h-1) (h-1:-1:0)];
                Ypatch1=[tmp(1,:) fliplr(tmp(5,:))];
                IRFpatch1=patch(Xpatch1,Ypatch1,[.65 .65 .65]);  %build confidence bands 90%
                set(IRFpatch1,'facealpha',.5);                   %add transparency
                set(IRFpatch1,'edgecolor','none');               %remove edge color

%                 Xpatch2=[(0:1:h-1) (h-1:-1:0)];
%                 Ypatch2=[tmp(2,:) fliplr(tmp(4,:))];
%                 IRFpatch2=patch(Xpatch2,Ypatch2,[.65 .65 .65]);  %build confidence bands 68%
%                 set(IRFpatch2,'facealpha',.5);                   %add transparency
%                 set(IRFpatch2,'edgecolor','none');               %remove edge color
                
                %plot(0:h-1,tmp(3,:),'r-','LineWidth',2);
                
            else
                plot(0:h-1,tmp(3,:),'Color',style(j,:),'LineWidth',1.5,'MarkerSize',3);
            end


        end

        plot([0,h-1],[0 0],'k--');                         %plot horizontal axis
        
        hold off

        % Set X-axis and Y-axis
        minband=min(tmp_margins(3,:));        %minimum value of lower bound
        maxband=max(tmp_margins(3,:));        %maximum value of upper bound
        space=maxband-minband;
        Ymin=minband-0.2*space-realmin;       %lower limit of the graph
        Ymax=maxband+0.2*space+realmin;       %upper limit of the graph
        set(gca,'XLim',[0 h-1],'YLim',[Ymin Ymax],'XTick',0:6:h-1,'FontName','Times New Roman');

        % Create titles
        title(longnames{1,i},'FontWeight','normal');  %variable shocked
%         title(varlist{1,i},'FontWeight','normal');  %variable shocked
        ylabel(ylab(i))
        
        % Create overall legend
        if i == (t-1)*rows*cols+1
            hL = legend(legNames,'Orientation','horizontal','NumColumns',5);
            position = [0.45 0.01 0.1 0.04];
            units = 'normalized';
            set(hL,'Position',position,'Units',units,'FontSize',10);
        end
        
        box on
        
        if i == N
            xlabel(['Horizon (',horz,')'])
            break
        end
        
    end
    
    % Top super-title
%     ax = axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
%     set(get(ax,'Title'),'Visible','on')
%     title(['Shock to: ',ShockName],'FontSize',11,'FontName','Times New Roman','FontWeight','normal');

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
    print(gcf,'-dpdf',[folder_,outFileName,'_',num2str(t),'.pdf']);

    if saveFig
        savefig(gcf,[folder_,outFileName,'_',num2str(t),'.fig']);
    end

    end
    
end %loop over figures

% Restore default groot styles and colours
% set(groot,'defaultAxesLineStyleOrder','remove')
% set(groot,'defaultAxesColorOrder','remove')

