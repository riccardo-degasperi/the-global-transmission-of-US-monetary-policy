function out = getSSA_irf(in,data,dims,dataOpt,plotOpt,ssaOpt)
%

%
% Last updated: 28/02/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Additional options:
plotIndividualChannels = 1; %plot conditional IRFs with bands
plotHorizon            = 1; %choose horizon for shock distribution chart


%-------------------------------------------------------------------------%
disp('--> Structural Scenario Analysis...')

% Unpack
betas           = in.betas;
B0nnr           = in.B0nnr ;
B0              = in.B0;
h               = dims.h;
N               = dims.N;
KLtresh         = ssaOpt.KLtresh;
varendo         = dataOpt.varendo;
varendolong     = dataOpt.varendolong;
varplot         = dataOpt.varplot;
shockVar        = dataOpt.shockVar;
units           = dataOpt.units;
plotOpt.ylab    = dataOpt.ylab;
ylab0           = dataOpt.ylab;
whichPlot       = plotOpt.whichPlot;
estimation      = plotOpt.estimation;
ssaOpt_         = ssaOpt;
impulse0         = ssaOpt.impulse;
if isfield(ssaOpt,'customLabels')
    customLabels = ssaOpt.customLabels;
else
    customLabels = {};
end

% Patch for chartNames when units are more than 1
if numel(units) > 1
    units0 = join(string(units),'_');
else
    units0 = units;
end

% Default chart labels
if isempty(customLabels)
    customLabels = cellstr([repmat('restriction',numel(ssaOpt.SSAy),1) num2str(vec(1:numel(ssaOpt.SSAy)))]);
end

%-------------------------------------------------------------------------%
% Loop over shocks
for k = 1:numel(impulse0)

    % Select shock
    impulse = impulse0(k);
    iImpulse = strcmp(varendo,impulse0(k));

    % Loop over scenarios
    for j = 1:numel(ssaOpt.SSAy)

        disp(' ')
        disp(['    Shock in position: ',impulse{:},'; Scenario: ',customLabels{j}])

        % Upload restrictions
        ssaOpt_.SSAy = ssaOpt.SSAy{j};
        ssaOpt_.SSAe = ssaOpt.SSAe{j};
        SSArest = feval(ssaOpt.SSAscript,impulse,dims,dataOpt,ssaOpt_);
        
        % Transform to compact forecast notation
        disp('    1. Transforming to compact forecast notation...')
        [b,M] = compactForecastNotation(B0nnr,betas,data,dims);
        
        % Draw Y | B,S,Q
        disp('    2. Drawing conditional forecasts...')
        [SSA,fail] = DrawConditionalForecast_irf(b,M,SSArest,iImpulse,dims);
        
        % Keep only draws that are plausible enough
        iKL = SSA.KLc <= KLtresh;
        y_cf2 = SSA.y_cf(:,:,iKL);
        y_bl2 = SSA.y_bl(:,:,iKL);
        effSample = sum(iKL);
        
        % Display effective sample
        disp(['       Effective sample respecting KL condition: ',num2str(effSample)])
        if effSample == 0
            error('ERROR: there are no valid draws respecting KL condition. Try with a looser "KLtresh" or with a less demanding scenario.')
        end
        
        %-----------------------------------------------------------------%
        % Get IRFs
        IRF_cf = squeeze(num2cell(zeros(h+1,effSample,N,N), [1 2]));
        IRF_bl = squeeze(num2cell(zeros(h+1,effSample,N,N), [1 2]));
        for i = 1:N
            for s = 1:effSample
            IRF_cf{i,iImpulse}(:,s) = y_cf2(:,i,s);
            IRF_bl{i,iImpulse}(:,s) = y_bl2(:,i,s);
            end
        end
        
        % Generate credibility regions
        IRFulm_cf = IRFbands(IRF_cf,N,h);
        IRFulm_bl = IRFbands(IRF_bl,N,h);
        
        % Store in structure
        IRFulm.IRFulm_bl = IRFulm_bl;
        IRFulm.IRFulm_cf = IRFulm_cf;
        
        
        %-----------------------------------------------------------------%
        % Plot charts
        disp('    3. Plotting scenarios...')
        
        % Save to folder SSA
        if j == 1
        folder0 = plotOpt.folder;
        plotOpt.folder = [plotOpt.folder,'/SSA/',SSArest.chooseSS];
        end
    
        % Define common chartname
        chartName0 = [char(units0),'_',estimation,'_',customLabels{j},'_SSA_',impulse{:}];
    
        switch whichPlot
    
            case 'none'
                
                %do nothing
        
            case 'all'        %plot all versions
        
            % Varendo
            chartName = [chartName0,'_varendo_overlaid'];
            legNames = {'baseline','counterfactual'};
            plotOpt.ylab = ylab0;
            IRFdisplay_SSA(IRFulm,varendo,varendolong,impulse,legNames,'allSeries',plotOpt,chartName)
            %IRFdisplay_stack(IRFstruc,varendo,varendolong,shockVar,legNames,'allSeries',plotOpt,chartName)
            % Use IRFdisplay_stack if you want bands for the conditional IRFs too.
        
            if plotIndividualChannels
            IRFdisplay_single(IRFulm_bl,varendo,varendolong,impulse,[chartName0,'_varendo_bl'],'allSeries',plotOpt)
            
            cb90_ = plotOpt.cb90;
            plotOpt.cb90 = 0;  %plot only 68% bands for conditional IRFs.
            IRFdisplay_single(IRFulm_cf,varendo,varendolong,impulse,[chartName0,'_varendo_cf'],'allSeries',plotOpt)
            plotOpt.cb90 = cb90_;
            end
        
            % Varplot
            if ~strcmp(varplot,impulse)
                varplot = [varplot,impulse];
            end
            iPlot = ismember(varendo,varplot);
            IRFulm_plot.IRFulm_bl = IRFulm_bl(iPlot,iPlot);
            IRFulm_plot.IRFulm_cf = IRFulm_cf(iPlot,iPlot);
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            plotOpt.ylab = ylab0(iPlot);
            chartName = [chartName0,'_varplot_overlaid'];
            IRFdisplay_SSA(IRFulm_plot,varplot,longnames,impulse,legNames,'varplot',plotOpt,chartName)
        
            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = ismember(varendo,impulse);
            iPlot2(posShock) = 1;        %reincluding policy variable
            IRFulm_plot.IRFulm_bl = IRFulm_bl(iPlot2,iPlot2);
            IRFulm_plot.IRFulm_cf = IRFulm_cf(iPlot2,iPlot2);
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            plotOpt.ylab = ylab0(iPlot2);
            chartName = [chartName0,'_nonvarplot_overlaid'];
            IRFdisplay_SSA(IRFulm_plot,nonvarplot,longnames,impulse,legNames,'varplot',plotOpt,chartName)
    
            case 'varendo'    %varendo only
        
            chartName = [chartName0,'_varendo_overlaid'];
            legNames = {'baseline','counterfactual'};
            plotOpt.ylab = ylab0;
            IRFdisplay_SSA(IRFulm,varendo,varendolong,impulse,legNames,'allSeries',plotOpt,chartName)
            %IRFdisplay_stack(IRFstruc,varendo,varendolong,shockVar,legNames,'allSeries',plotOpt,chartName)
            % Use IRFdisplay_stack if you want bands for the conditional IRFs too.
        
            if plotIndividualChannels
            IRFdisplay_single(IRFulm_bl,varendo,varendolong,impulse,[chartName0,'_varendo_bl'],'allSeries',plotOpt)
            
            cb90_ = plotOpt.cb90;
            plotOpt.cb90 = 0;  %plot only 68% bands for conditional IRFs.
            IRFdisplay_single(IRFulm_cf,varendo,varendolong,impulse,[chartName0,'_varendo_cf'],'allSeries',plotOpt)
            plotOpt.cb90 = cb90_;
            end
        
            case 'varplot1'    %varplot only
        
            % Varplot
            if ~strcmp(varplot,impulse)
                varplot = [varplot,impulse];
            end
            iPlot = ismember(varendo,varplot);
            IRFulm_plot.IRFulm_bl = IRFulm_bl(iPlot,iPlot);
            IRFulm_plot.IRFulm_cf = IRFulm_cf(iPlot,iPlot);
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            plotOpt.ylab = ylab0(iPlot);
            chartName = [chartName0,'_varplot_overlaid'];
            legNames = {'baseline','counterfactual'};
            IRFdisplay_SSA(IRFulm_plot,varplot,longnames,impulse,legNames,'varplot',plotOpt,chartName)
            
            if plotIndividualChannels
            plotOpt.ylab = ylab0;
            IRFdisplay_single(IRFulm_bl,varendo,varendolong,impulse,[chartName0,'_varendo_bl'],'allSeries',plotOpt)
            
            cb90_ = plotOpt.cb90;
            plotOpt.cb90 = 0;  %plot only 68% bands for conditional IRFs.
            IRFdisplay_single(IRFulm_cf,varendo,varendolong,impulse,[chartName0,'_varendo_cf'],'allSeries',plotOpt)
            plotOpt.cb90 = cb90_;
            end

            case 'varplot2'    %varplot only
        
            % Varplot
            if ~strcmp(varplot,impulse)
                varplot = [varplot,impulse];
            end
            iPlot = ismember(varendo,varplot);
            IRFulm_plot.IRFulm_bl = IRFulm_bl(iPlot,iPlot);
            IRFulm_plot.IRFulm_cf = IRFulm_cf(iPlot,iPlot);
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            plotOpt.ylab = ylab0(iPlot);
            chartName = [chartName0,'_varplot_overlaid'];
            legNames = {'baseline','counterfactual'};
            IRFdisplay_SSA(IRFulm_plot,varplot,longnames,impulse,legNames,'varplot',plotOpt,chartName)
        
            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = ismember(varendo,impulse);
            iPlot2(posShock) = 1;        %reincluding policy variable
            IRFulm_plot.IRFulm_bl = IRFulm_bl(iPlot2,iPlot2);
            IRFulm_plot.IRFulm_cf = IRFulm_cf(iPlot2,iPlot2);
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            plotOpt.ylab = ylab0(iPlot2);
            chartName = [chartName0,'_nonvarplot_overlaid'];
            IRFdisplay_SSA(IRFulm_plot,nonvarplot,longnames,impulse,legNames,'varplot',plotOpt,chartName)
    
            if plotIndividualChannels
            plotOpt.ylab = ylab0;
            IRFdisplay_single(IRFulm_bl,varendo,varendolong,impulse,[chartName0,'_varendo_bl'],'allSeries',plotOpt)
            
            cb90_ = plotOpt.cb90;
            plotOpt.cb90 = 0;  %plot only 68% bands for conditional IRFs.
            IRFdisplay_single(IRFulm_cf,varendo,varendolong,impulse,[chartName0,'_varendo_cf'],'allSeries',plotOpt)
            plotOpt.cb90 = cb90_;
            end

        end

        % Plot distribution of the shocks vs. normal
        chartName = [chartName0,'_shockDistribution'];
        SSAplotDistribution(SSA.shocks,plotHorizon,varendolong,chartName,plotOpt)
    
        % Pack output
        if numel(ssaOpt.SSAy) == 1
            out.(impulse{:}).IRFsulm       = IRFulm;
            out.(impulse{:}).effSample     = effSample;
            out.(impulse{:}).iKL           = iKL;
            out.(impulse{:}).KLc           = SSA.KLc;
            out.(impulse{:}).shockSequence = SSA.shocks;
        else
    
            % Remove whitespace
            customLabels_ = cellfun(@(x) strrep(x, ' ', '_'), customLabels, 'UniformOutput', false);
            customLabels_ = cellfun(@(x) strrep(x, '.', ''), customLabels_, 'UniformOutput', false);
            customLabels_ = cellfun(@(x) strrep(x, ',', ''), customLabels_, 'UniformOutput', false);
            customLabels_ = cellfun(@(x) strrep(x, '-', '_'), customLabels_, 'UniformOutput', false);

            out.(impulse{:}).(customLabels_{j}).IRFsulm       = IRFulm;
            out.(impulse{:}).(customLabels_{j}).effSample     = effSample;
            out.(impulse{:}).(customLabels_{j}).iKL           = iKL;
            out.(impulse{:}).(customLabels_{j}).KLc           = SSA.KLc;
            out.(impulse{:}).(customLabels_{j}).shockSequence = SSA.shocks;
    
            % Prepare structure for overlaid plot
            chIRFsulm.(customLabels_{j}) = IRFulm_cf;
            chIRFs.(customLabels_{j}) = IRF_cf;
            if j == numel(ssaOpt.SSAy)
                chIRFsulm.baseline = IRFulm_bl;
                chIRFs.baseline = IRF_bl;
            end
        end

    end %loop over scenarios

    % Restore folder path
    plotOpt.folder = folder0;
    
    % Overlaid scenarios plot
    if numel(ssaOpt.SSAy) > 1
        if numel(customLabels) == numel(ssaOpt.SSAy)
            customLabels = [customLabels;'baseline'];   %add label for baseline if not present
        end
        chartName = [char(units0),'_',estimation,'_SSA_',SSArest.chooseSS,'_',impulse{:},'_all'];
        plotOpt.ylab = ylab0;
        style = {'k-','k--','k-.','k:','ko','r-'};
        IRFdisplay_structure(chIRFsulm,varendo,varendolong,impulse,customLabels,style,'allSeries',plotOpt,chartName)
    
        % Pack output
        out.(impulse{:}).chIRFs = chIRFs;
    end

    % Pack output
    out.(impulse{:}).restrictions = SSArest;

end %loop over shocks

% Leave some space
disp(' ')

