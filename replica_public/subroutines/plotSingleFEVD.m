function [] = plotSingleFEVD(FEVD,dataOpt,plotOpt,chartName)

    % Unpack
    varendo     = dataOpt.varendo;
    shockVar    = dataOpt.shockVar;
    longnames   = dataOpt.varendolong;
    saveCharts  = plotOpt.saveCharts;
    saveFig     = plotOpt.saveFig;
    folder      = plotOpt.folder;
    horz        = plotOpt.frequency;
    font        = plotOpt.font;
    row         = plotOpt.allSeries.plotRows;
    col         = plotOpt.allSeries.plotCols;

    % Get dimensions and positions
    [h,N,N] = size(FEVD.FEVDp50);
    iShock = contains(varendo,shockVar); %logical position of shock variables
    iShock = find(iShock);               %position of shock variables
    shockName = varendo(iShock);         %name of the shock variables
    Ns = numel(iShock);                  %number of shocks
    

    % Loop over shocks
    for j = 1:Ns

        % Select shock
        iShockj = false(1,N);
        iShockj(iShock(j)) = 1;

        % Repack FEVD and confidence regions
        FEVDulm = cell(N,1);
        for i = 1:N
        FEVDulm{i} = [FEVD.FEVDp5(:,iShockj,i) FEVD.FEVDp16(:,iShockj,i) FEVD.FEVDp50(:,iShockj,i) FEVD.FEVDp84(:,iShockj,i) FEVD.FEVDp95(:,iShockj,i)]';
        end
        
        % Plot FEVD and regions
%         row = round(sqrt(N));
%         col = ceil(sqrt(N));
        nn = length(findobj('type','figure'));
        figure(nn+1)
        for ii = 1:N
            subplot(row,col,ii);
            tmp = FEVDulm{ii};
    
            hold on
    
            % Plot confidence bands
            Xpatch1=[(0:1:h-1) (h-1:-1:0)];
            Ypatch1=[tmp(1,:) fliplr(tmp(5,:))];
            IRFpatch1=patch(Xpatch1,Ypatch1,[0.3020 0.7451 0.9333]);  %build confidence bands 90%
            set(IRFpatch1,'facealpha',.5);                   %add transparency
            set(IRFpatch1,'edgecolor','none');               %remove edge color
    
            Xpatch2=[(0:1:h-1) (h-1:-1:0)];
            Ypatch2=[tmp(2,:) fliplr(tmp(4,:))];
            IRFpatch2=patch(Xpatch2,Ypatch2,[0 0.4471 0.7412]);  %build confidence bands 68%
            set(IRFpatch2,'facealpha',.5);                   %add transparency
            set(IRFpatch2,'edgecolor','none');               %remove edge color
    
            % Plot median and X-axis
            plot(0:h-1,tmp(3,:),'Color',[.2 .4 .6],'LineWidth',2); %plot median
            %plot([0 h-1],[0 0],'k--');                         %plot horizontal axis
    
            hold off
    
            % Set X-axis and Y-axis
            minband=min(tmp(1,:));        %minimum value of lower bound
            maxband=max(tmp(5,:));        %maximum value of upper bound
            space=maxband-minband;
            %Ymin=minband-0.2*space-realmin;       %lower limit of the graph
            Ymax=maxband+0.2*space+realmin;       %upper limit of the graph
            set(gca,'XLim',[0 h-1],'YLim',[0 Ymax],'XTick',0:6:h-1,'FontName',font);
            box on
            grid on
    
            % Create titles
            title(longnames{1,ii},'FontWeight','normal');  %variable shocked
            ylabel('%');
            if ii == N
                xlabel(['Horizon (',horz,')'])
                break
            end
    
            % Create overall legend
    %         if i == (t-1)*rows*cols+1
    %             hL = legend('90%','68%','median');
    %             position = [0.9 0.1 0.03 0.03];
    %             units = 'normalized';
    %             set(hL,'Position', position,'Units', units);
    %         end
    
        end
    
        % Top super-title
        ax = axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        title('Forecast Error Variance Decomposition','FontSize',11,'FontName',font,'FontWeight','bold');
        
        % Save plot
        if saveCharts
    
            folder_ = [pwd '/plots/',folder,'/FEVD'];
            if ~isempty(folder)
                if not(isfolder(folder_))
                    mkdir(folder_)
                end
                folder_ = [folder_,'/'];
            end
            
            set(gcf,'PaperUnits','centimeters','PaperSize',[col*6.5 row*6.5]) %[x y]
            set(gcf,'PaperPosition',[0.25*col 0.25*row 5.75*col 5.75*row])    %[left bottom width height]
            print(gcf,'-dpdf',[folder_,chartName,'_',shockName{j},'.pdf']);
            
            if saveFig
            savefig(gcf,[folder_,chartName,'_',shockName{j},'.fig']);
            end
        end

    end

end