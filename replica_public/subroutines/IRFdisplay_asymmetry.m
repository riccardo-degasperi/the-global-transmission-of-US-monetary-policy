function [ ] = IRFdisplay_asymmetry(IRFstruc,varlist,longnames,shockVar,legNames,plotCase,plotOpt,chartName0)
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
horz        = plotOpt.frequency;
ylab        = plotOpt.ylab;
cb90        = plotOpt.cb90;
font        = plotOpt.font;

% Patch: compatibility with TVAR code
chartName = chartName0;

% Colors & Styles
style = {'-','--',':'};
cbands = [.4 .4 .4;.6 .6 .6;.8 .8 .8];
cmedian = [1 0.4 0.15;0 0.45 0.75;0.25 0.5 0.15];

% Position of shock & name
if ischar(shockVar)
    shockVar = cellstr(shockVar);
end

% Loop over shocks
for ii = 1:numel(shockVar)

    iShock = strcmp(varlist,shockVar{ii}); %position of shock variable
    ShockName = longnames{iShock};         %name of the shock variable
    
    if numel(shockVar) > 1
        chartName = [chartName0,'_',shockVar{ii}];
    end

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
    for t = 1:Nfig
    
        irf = figure(t+nn);      %open figure
        set(irf,'name',chartName);
        count = 0;
    
        for i = (t-1)*rows*cols+1:t*rows*cols
            
            count = count + 1;
    
            subplot(rows,cols,count)
    
            hold on
            
            tmp0 = [];  %initialize tmp0, usedd to compute max and min height below
            go = [];    %initialize graphic object container for legend
            
            for j = 1:numel(names)
            
                tmp = IRFstruc.(names{j}){i,iShock};
                tmp0 = [tmp0 tmp];
    
                % Plot confidence bands
                Xpatch1=[(0:1:h-1) (h-1:-1:0)];
                if cb90
                Ypatch1=[tmp(1,:) fliplr(tmp(5,:))];             %build confidence bands 90%
                else
                Ypatch1=[tmp(2,:) fliplr(tmp(4,:))];             %build confidence bands 68%
                end
                IRFpatch1=patch(Xpatch1,Ypatch1,cbands(j,:));    
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
    
            xtickangle(0)
    
            % Create titles
            title(longnames{1,i},'FontWeight','normal');  %variable shocked
    %         title(varlist{1,i},'FontWeight','normal');  %variable shocked
            ylabel(ylab(i))
    
            % Create overall legend
            if i == (t-1)*rows*cols+1
                hL = legend(go,legNames,'FontSize',5);
    %             position = [0.8 0.15 0.1 0.1];              %comment out to have legend inside first subplot
                units = 'normalized';
                set(hL,'Units',units);
    %             set(hL,'Position',position,'Units',units);  %comment out to have legend inside first subplot
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
    %     title(['Shock to: ',ShockName],'FontSize',11,'FontName',font,'FontWeight','normal');
    
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

end %loop over shocks


