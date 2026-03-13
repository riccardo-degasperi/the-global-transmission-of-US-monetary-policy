function [ ] = IRFdisplay_overlay(IRFstruc,varlist,iskip,longnames,shockVar,legNames,plotCase,plotOpt,chartName)
%IRFdisplay_asymmetry plots the IRFs, overlaying for each variable the
%various cases specified in IRFstruc as elements of the structure. Each
%element of the structure is a NxN cell. Each cell contains a horizons
%by simulations matrices of IRFs.

% INPUT:
% IRFstruc  : structure with as many fields as overlays (E.g. IRFs identified using 'positive' and 'negative' MP shocks)
% varlist   : variables to plot (identifiers)
% longnames : names of variables to plot
% shockVar  : policy variable
% legNames  : legend entries
% plotCase  : selects an option out of plotOpt
% plotOpt   : structure containing plot options
% chartName : output file name

% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%


% Unpacking
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
rows0       = plotOpt.(plotCase).plotRows;
cols        = plotOpt.(plotCase).plotCols;
folder      = plotOpt.folder;
cb90        = plotOpt.cb90;
horz        = plotOpt.frequency;
ylab        = plotOpt.ylab;
font        = plotOpt.font;

% Colors & Styles
style = {'-','--',':'};

cbands = [.55 .55 .55 .45 .45 .45;.75 .75 .75 .65 .65 .65];
%cbands = [.4 .5 .5 .3 .5 .5;.8 .7 .7 .8 .6 .5];
cmedian = [0 .5 .5;.8 .3 .1];

% Position of shock & name
iShock = strcmp(varlist,shockVar); %position of shock variable
ShockName = longnames(iShock);     %name of the shock variable

% iShock for EM
iShock2 = iShock(iskip);

% Determine dimensions
names = fieldnames(IRFstruc);
N = size(varlist,2);
h = size(IRFstruc.(names{1}){1,iShock},2);

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
countEM = 0;
for t = 1:Nfig

    irf = figure(t+nn);      %open figure
    set(irf,'name',chartName);
    count = 0;
    

    for i = (t-1)*rows*cols+1:t*rows*cols
        
        count = count + 1;
        countEM = countEM +1;

        subplot(rows,cols,count)

        hold on
        
        tmp0 = [];  %initialize tmp0, usedd to compute max and min height below
        go = [];    %initialize graphic object container for legend
        
        for j = 1:numel(names)

            if j == 2 && iskip(i) == 0
                countEM = countEM - 1;
                break
            end

            if j == 1
            tmp = IRFstruc.(names{j}){i,iShock};
            elseif j == 2
            tmp = IRFstruc.(names{j}){countEM,iShock2};
            end
            tmp0 = [tmp0 tmp];

            % Plot confidence bands
            if cb90
            Xpatch1=[(0:1:h-1) (h-1:-1:0)];
            Ypatch1=[tmp(1,:) fliplr(tmp(5,:))];
            IRFpatch1=patch(Xpatch1,Ypatch1,cbands(j,1:3));  %build confidence bands 90%
            set(IRFpatch1,'facealpha',.5);                   %add transparency
            set(IRFpatch1,'edgecolor','none');               %remove edge color
            end
            
            Xpatch1=[(0:1:h-1) (h-1:-1:0)];
            Ypatch1=[tmp(2,:) fliplr(tmp(4,:))];
            IRFpatch1=patch(Xpatch1,Ypatch1,cbands(j,4:6));  %build confidence bands 68%
            set(IRFpatch1,'facealpha',.5);                   %add transparency
            set(IRFpatch1,'edgecolor','none');               %remove edge color

            % Plot median and X-axis
            go(j) = plot(0:h-1,tmp(3,:),style{j},'LineWidth',2,'Color',cmedian(j,:)); %plot median


        end

        plot([0,h-1],[0 0],'k-');                         %plot horizontal axis
        
        hold off
        
        % Set X-axis and Y-axis
        if cb90
        minband=min(tmp0(1,:));        %minimum value of lower bound
        maxband=max(tmp0(5,:));        %maximum value of upper bound
        else
        minband=min(tmp0(2,:));        %minimum value of lower bound
        maxband=max(tmp0(4,:));        %maximum value of upper bound
        end
        space=maxband-minband;
        Ymin=minband-0.2*space-realmin;       %lower limit of the graph
        Ymax=maxband+0.2*space+realmin;       %upper limit of the graph
        set(gca,'XLim',[0 h-1],'YLim',[Ymin Ymax],'XTick',0:6:h-1,'FontName',font);

        ax = gca;
        ax.YAxis.Exponent = 0;
        xtickangle(0)
        
        % Create titles
        title(longnames{1,i},'FontWeight','normal');  %variable shocked
%         title(varlist{1,i},'FontWeight','normal');  %variable shocked
        ylabel(ylab(i))

        % Create overall legend
        if i == (t-1)*rows*cols+1
            hL = legend(go,legNames,'Location','southeast');
%             position = [0.8 0.15 0.1 0.1];              %comment out to have legend inside first subplot
            units = 'normalized';
            set(hL,'Units',units);
%             set(hL,'Position',position,'Units',units);  %comment out to have legend inside first subplot
        end
        
        box on
        
        if i == N || i == t*rows*cols
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

