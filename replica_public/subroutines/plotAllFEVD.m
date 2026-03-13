function [] = plotAllFEVD(meanFEVD,dataOpt,plotOpt,legNames,chartName)

    % Unpack
    varendo     = dataOpt.varendo;
    shockVar    = dataOpt.shockVar;
    saveCharts  = plotOpt.saveCharts;
    saveFig     = plotOpt.saveFig;
    folder      = plotOpt.folder;
    font        = plotOpt.font;

    % Get dimensions
    [h,N,N] = size(meanFEVD);
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
            error('ERROR: check custom legend for FEVD chart in file "main".')
        end
        if numel(legNames)<N
            legNames(end+1) = {'Not identified'};
        end
    end

    % Sum contributions of unidentified shocks into one category
    FEVDcompressed = [meanFEVD(:,iShock,:) sum(meanFEVD(:,~iShock,:),2)];

    % Plot FEVD
    row = round(sqrt(N));
    col = ceil(sqrt(N));
    nn = length(findobj('type','figure'));
    figure(nn+1)
    tcl = tiledlayout(row,col);
    for ii = 1:N
        nexttile(tcl);
        AreaPlot(squeeze(FEVDcompressed(:,:,ii)));
        xlim([1 h]); ylim([0 100]);
        title(varendo{ii}, 'FontWeight','bold','FontSize',10); 
        set(gca,'Layer','top','FontName',font);
    end

    % Create overall legend
    hL = legend(legNames);
    hL.Layout.Tile = 'south';

    % Save plot
    if saveCharts

        folder_ = [pwd '/plots/',folder,'/FEVD'];
        if ~isempty(folder)
            if not(isfolder(folder_))
                mkdir(folder_)
            end
            folder_ = [folder_,'/'];
        end
        
        leg = 0.5*numel(legNames);
        set(gcf,'PaperUnits','centimeters','PaperSize',[col*6.5 row*6.5+leg]) %[x y]
        set(gcf,'PaperPosition',[0.25*col 0.25*row 5.75*col 5.75*row+leg])    %[left bottom width height]
        print(gcf,'-dpdf',[folder_,chartName,'.pdf']);
        
        if saveFig
        savefig(gcf,[folder_,chartName,'.fig']);
        end
    end



end