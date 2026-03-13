function [] = plotSingleHD(HD,data,dataOpt,plotOpt,chartName)

% OPTIONS
type = 'shock-v-all';               %type = {'shock-v-all','shock-v-data'}. Plot shock contribution versus demeaned series, or versus sum of contributions of all shocks.
nticks = 6;                         %number of ticks to display
%-------------------------------------------------------------------------%

% Unpack
varendo     = dataOpt.varendo;
varendolong = dataOpt.varendolong;
shockVar    = dataOpt.shockVar;
frequency   = dataOpt.frequency;
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
folder      = plotOpt.folder;
font        = plotOpt.font;
dates       = data.dates;
Y           = data.Y;

% Get dimensions
[T,N] = size(Y);                     %get dimensions (excluding initial obs)
p = size(dates,1)-size(Y,1);         %get number of lags
iShock = contains(varendo,shockVar); %logical position of shock variables
iShock = find(iShock);               %position of shock variables
shockName = varendo(iShock);         %name of the shock variables
Ns = numel(iShock);                  %number of shocks

% Loop over shocks
for j = 1:Ns

    iShockj = false(1,N);
    iShockj(iShock(j)) = 1;

    % Initialise containers
    HDulm = cell(N,1);                   %for shock contributions
    HDsum = nan(size(Y));                %for sum of shocks contributions
    yplot = nan(size(Y));                %for demeaned policy variable
    
    % Repack HD and confidence regions
    for i = 1:N
    
        % Contribution of shock of interest
        HDulm{i} = [HD.shockp5(p+1:end,iShockj,i) HD.shockp16(p+1:end,iShockj,i) HD.shockp50(p+1:end,iShockj,i) HD.shockp84(p+1:end,iShockj,i) HD.shockp95(p+1:end,iShockj,i)]';
        
        % Contribution of all shocks
        HDsum(:,i) = sum(HD.shock(p+1:end,:,i),2);
    
        % Demeaned true series
        yplot(:,i) = Y(:,i) - mean(Y(:,i),1);
    end
    
    % Get ticks
    dates = dates(p+1:end);              %drop initial observations
    if strcmp(frequency,'yearly')
        [years,ia,ic] = unique(dates);
        if numel(years) <= nticks
            skip = 1;
        else
            skip = floor(numel(years)/nticks);
        end
        tick_pos = 1:skip:numel(years);
        ticklabels = years(tick_pos);
        ticks = ia(tick_pos);
    elseif strcmp(frequency,'monthly')
        [years,ia,ic] = unique(extractBefore(dates,5));
        if numel(years) <= nticks
            skip = 1;
        else
            skip = floor(numel(years)/nticks);
        end
        tick_pos = 1:skip:numel(years);
        ticklabels = years(tick_pos);
        ticks = ia(tick_pos);
    elseif strcmp(frequency,'weekly') || strcmp(frequency,'daily')
        [years,ia,ic] = unique(extractAfter(dates,6));
        if numel(years) <= nticks
            skip = 1;
        else
            skip = floor(numel(years)/nticks);
        end
        tick_pos = 1:skip:numel(years);
        ticklabels = years(tick_pos);
        ticks = ia(tick_pos);
    else
        error('ERROR: dataOpt.frequency can be "daily", "weekly", "monthly" or "yearly".')
    end
    
    % Plot
    %==========================================================================
    xdim = 16; %cm
    ydim = 9;  %cm
    nn = length(findobj('type','figure'));
    for f = 1:N
        figure(nn+f)
    
        tmp_hd = HDulm{f,1};
        switch type
            case 'shock-v-data'
            tmp_y = yplot(:,f);
            case 'shock-v-all'
            tmp_y = HDsum(:,f);
        end
    
        hold on
    
        % Plot confidence bands
        Xpatch1=[(0:1:T-1) (T-1:-1:0)];
        Ypatch1=[tmp_hd(1,:) fliplr(tmp_hd(5,:))];
        IRFpatch1=patch(Xpatch1,Ypatch1,[.75 .75 .75]);  %build confidence bands 90%
        set(IRFpatch1,'facealpha',.5);                   %add transparency
        set(IRFpatch1,'edgecolor','none');               %remove edge color
    
        % Plot confidence bands
        Xpatch2=[(0:1:T-1) (T-1:-1:0)];
        Ypatch2=[tmp_hd(2,:) fliplr(tmp_hd(4,:))];
        IRFpatch2=patch(Xpatch2,Ypatch2,[.65 .65 .65]);  %build confidence bands 68%
        set(IRFpatch2,'facealpha',.5);                   %add transparency
        set(IRFpatch2,'edgecolor','none');               %remove edge color
    
        % Plot median and X-axis
        plot(0:T-1,tmp_hd(3,:),'Color',[.2 .4 .6],'LineWidth',2); %plot median
        plot([0 T-1],[0 0],'k--');                                %plot horizontal axis
    
        % Plot true demeaned series / sum of shocks
        plot(0:T-1,tmp_y,'Color','k','LineWidth',2);
    
        hold off
    
        % Set X-axis and Y-axis
        minband=min(tmp_y);        %minimum value of lower bound
        maxband=max(tmp_y);        %maximum value of upper bound
        space=maxband-minband;
        Ymin=minband-0.2*space-realmin;       %lower limit of the graph
        Ymax=maxband+0.2*space+realmin;       %upper limit of the graph
        set(gca,'FontName',font);
        xticks(ticks)
        xticklabels(ticklabels)
        xlim([1 T]);
        ylim([Ymin Ymax]);
        title(['Historical decomposition: ',varendolong{f}], 'FontWeight','bold','FontSize',10);
        box on
        grid on
        
        % Save plot
        if saveCharts
            folder_ = [pwd '/plots/',folder,'/HD'];
            if ~isempty(folder)
                if not(isfolder(folder_))
                    mkdir(folder_)
                end
                folder_ = [folder_,'/'];
            end
            
            set(gcf,'PaperUnits','centimeters','PaperSize',[xdim ydim]) %[x y]
            set(gcf,'PaperPosition',[0 1 xdim ydim-2])    %[left bottom width height]
            print(gcf,'-dpdf',[folder_,chartName,'_of_',varendo{f},'_to_',shockName{j},'.pdf']);
            
            if saveFig
            savefig(gcf,[folder_,chartName,'_of_',varendo{f},'_to_',shockName{j},'.fig']);
            end
        end
    end

end

end