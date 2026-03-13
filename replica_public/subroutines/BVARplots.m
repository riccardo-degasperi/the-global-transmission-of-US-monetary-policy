function [] = BVARplots(IRFsulm,chartName0,dataOpt,plotOpt)
% Plots impulse response functions from a VAR model

% INPUT:
% - IRFsulm    : structural IRFs (median and bands);
% - chartName0 : first part of figure file name;
% - dataOpt    : structure containing data options; 
% - plotOpt    : structure containing plot options;

% OUTPUT:
% - plots of responses;


% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%

% Unpacking
plotType       = plotOpt.plotType;
whichPlot      = plotOpt.whichPlot;
varendo        = dataOpt.varendo;
varendolong    = dataOpt.varendolong;
varplot        = dataOpt.varplot;
shockVar       = dataOpt.shockVar;
plotOpt.ylab   = dataOpt.ylab;
ylab0          = dataOpt.ylab;



%-------------------------------------------------------------------------%
% Plot Structural-form IRFs

switch plotType

%=========================================================================%
    case 'allShocks'  %plots shocks to all variables

        chartName = [chartName0,'_allShocks'];
        IRFdisplay_all(IRFsulm,varendolong,chartName,plotOpt)

        %In IRFdisplay_all [Ymin Ymax] must be Ymax > Ymin. For shocks to
        %non policy variables, Ymax = Ymin and the code crashes.
        
        %TODO
        %Drop plots if not identified.
        %If you allow for identification of m shocks, plot only those.
%=========================================================================%


    case 'singleShock'  %plots shocks to variable shockVar


        switch whichPlot
            
            case 'none'
                
            %do nothing
        
            case 'all'        %plot all versions
        
            % Varendo
            chartName = [chartName0,'_varendo'];
            IRFdisplay_single(IRFsulm,varendo,varendolong,shockVar,chartName,'allSeries',plotOpt)
        
            % Varplot
            if sum(ismember(shockVar,varplot)) < numel(shockVar)
                iS = ismember(shockVar,varplot)==0;
                varplot = [varplot,shockVar(iS)];
            end
            iPlot = ismember(varendo,varplot);
            IRFsulm_plot = IRFsulm(iPlot,iPlot);
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            plotOpt.ylab = ylab0(iPlot);
            chartName = [chartName0,'_varplot'];
            IRFdisplay_single(IRFsulm_plot,varplot,longnames,shockVar,chartName,'varplot',plotOpt)
        
            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = ismember(varendo,shockVar);
            iPlot2(posShock) = 1;        %reincluding policy variable
            IRFsulm_plot = IRFsulm(iPlot2,iPlot2);
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            plotOpt.ylab = ylab0(iPlot2);
            chartName = [chartName0,'_nonvarplot'];
            IRFdisplay_single(IRFsulm_plot,nonvarplot,longnames,shockVar,chartName,'varplot',plotOpt)
        
            case 'varendo'    %varendo only
        
            chartName = [chartName0,'_varendo'];
            IRFdisplay_single(IRFsulm,varendo,varendolong,shockVar,chartName,'allSeries',plotOpt)
        
            case 'varplot1'    %varplot only
        
            % Varplot
            if sum(ismember(shockVar,varplot)) < numel(shockVar)
                iS = ismember(shockVar,varplot)==0;
                varplot = [varplot,shockVar(iS)];
            end
            iPlot = ismember(varendo,varplot);
            IRFsulm_plot = IRFsulm(iPlot,iPlot);
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            plotOpt.ylab = ylab0(iPlot);
            chartName = [chartName0,'_varplot'];
            IRFdisplay_single(IRFsulm_plot,varplot,longnames,shockVar,chartName,'varplot',plotOpt)
        
            case 'varplot2'    %varplot only
        
            % Varplot
            if sum(ismember(shockVar,varplot)) < numel(shockVar)
                iS = ismember(shockVar,varplot)==0;
                varplot = [varplot,shockVar(iS)];
            end
            iPlot = ismember(varendo,varplot);
            IRFsulm_plot = IRFsulm(iPlot,iPlot);
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            plotOpt.ylab = ylab0(iPlot);
            chartName = [chartName0,'_varplot'];
            IRFdisplay_single(IRFsulm_plot,varplot,longnames,shockVar,chartName,'varplot',plotOpt)
        
            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = ismember(varendo,shockVar);
            iPlot2(posShock) = 1;        %reincluding policy variable
            IRFsulm_plot = IRFsulm(iPlot2,iPlot2);
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            plotOpt.ylab = ylab0(iPlot2);
            chartName = [chartName0,'_nonvarplot'];
            IRFdisplay_single(IRFsulm_plot,nonvarplot,longnames,shockVar,chartName,'varplot',plotOpt)
        
        end

        % Restore default
        plotOpt.ylab = ylab0;
end

