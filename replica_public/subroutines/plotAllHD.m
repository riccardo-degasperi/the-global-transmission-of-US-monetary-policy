function [] = plotAllHD(HD,data,dataOpt,plotOpt,legNames,chartName)

% OPTIONS
nticks = 6;                         %number of ticks to display


% Unpack
varendo     = dataOpt.varendo;
shockVar    = dataOpt.shockVar;
frequency   = dataOpt.frequency;
saveCharts  = plotOpt.saveCharts;
saveFig     = plotOpt.saveFig;
folder      = plotOpt.folder;
font        = plotOpt.font;
dates       = data.dates;
Y           = data.Y;

% Get dimensions
[T,N,N] = size(HD.shock);
p       = size(dates,1)-size(Y,1);        %get number of lags
iShock  = contains(varendo,shockVar);     %logical position of shock variables
Ns      = sum(iShock);

% Get legend entries
if isempty(legNames)
    legNames = varendo(iShock);
    if numel(legNames)<N
        legNames(end+1) = {'Not identified'};
    end
else
    if numel(legNames) ~= Ns
        error('ERROR: check custom legend for HD chart in file "main".')
    end
    if numel(legNames)<N
        legNames(end+1) = {'Not identified'};
    end
end

% Sum contributions of unidentified shocks into one category
HDcompressed = [HD.shock(p+1:end,iShock,:) sum(HD.shock(p+1:end,~iShock,:),2)];

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

%% Plot
%==========================================================================
xdim = 16; %cm
ydim = 9;  %cm
nn = length(findobj('type','figure'));
for f = 1:N
    figure(nn+f)
    AreaPlot(squeeze(HDcompressed(:,:,f))); hold on; 
    plot(sum(squeeze(HDcompressed(:,:,f)),2),'-k','LineWidth',2);
    xticks(ticks)
    xticklabels(ticklabels)
    xlim([1 T-p]);
    title([varendo{f}], 'FontWeight','bold','FontSize',10);
	set(gca,'Layer','top','FontName',font);
    grid on

    % Create overall legend
    hL = legend(legNames);
    set(hL,'Location','southoutside');
    
    % Save plot
    if saveCharts
        folder_ = [pwd '/plots/',folder,'/HD'];
        if ~isempty(folder)
            if not(isfolder(folder_))
                mkdir(folder_)
            end
            folder_ = [folder_,'/'];
        end
        
        leg = 0.5*numel(legNames);
        set(gcf,'PaperUnits','centimeters','PaperSize',[xdim ydim+leg]) %[x y]
        set(gcf,'PaperPosition',[0 1 xdim ydim-2+leg])    %[left bottom width height]
        print(gcf,'-dpdf',[folder_,chartName,'_',varendo{f},'.pdf']);
        
        if saveFig
        savefig(gcf,[folder_,chartName,'_',varendo{f},'.fig']);
        end
    end
end



end