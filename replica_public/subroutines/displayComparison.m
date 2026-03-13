function [] = displayComparison(IRFsulm,dataOpt,plotOpt)
        
        %-------------------------------------------------------------------------%
        disp('--> Plotting comparison charts: symmetric vs. asymmetric priors...')

        % Unpacking
        shockVar        = dataOpt.shockVar;
        varendo         = dataOpt.varendo;
        varendolong     = dataOpt.varendolong;
        varplot         = dataOpt.varplot;
        units           = dataOpt.units;
        plotOpt.ylab    = dataOpt.ylab;
        whichPlot       = plotOpt.whichPlot;
        estimation      = plotOpt.estimation;

        % Patch for chartNames when units are more than 1
        if numel(units) > 1
            units0 = join(string(units),'_');
        else
            units0 = units;
        end
        
        % Get struc names
        names = fieldnames(IRFsulm);

        switch whichPlot

            case 'none'

            case 'all'        %plot all versions

            % Varendo
            chartName = [char(units0),'_DummyObs-v-',estimation,'_comparison_varendo'];
            legNames = {'Symmetric Prior','Asymmetric Prior'};
            IRFdisplay_asymmetry(IRFsulm,varendo,varendolong,shockVar,legNames,'allSeries',plotOpt,chartName);

            % Varplot
            if sum(ismember(shockVar,varplot)) < numel(shockVar)
                iS = ismember(shockVar,varplot)==0;
                varplot = [varplot,shockVar(iS)];
            end
            iPlot = ismember(varendo,varplot);
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot,iPlot);
            end
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            chartName = [char(units0),'_DummyObs-v-',estimation,'_comparison_varplot'];
            IRFdisplay_asymmetry(IRFsulm_plot,varplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = ismember(varendo,shockVar);
            iPlot2(posShock) = 1;        %reincluding policy variable
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot2,iPlot2);
            end
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            chartName = [char(units0),'_DummyObs-v-',estimation,'_comparison_nonvarplot'];
            IRFdisplay_asymmetry(IRFsulm_plot,nonvarplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

            case 'varendo'    %varendo only

            chartName = [char(units0),'_DummyObs-v-',estimation,'_comparison_varendo'];
            legNames = {'Symmetric Prior','Asymmetric Prior'};
            IRFdisplay_asymmetry(IRFsulm,varendo,varendolong,shockVar,legNames,'allSeries',plotOpt,chartName)

            case 'varplot1'    %varplot only

            % Varplot
            if sum(ismember(shockVar,varplot)) < numel(shockVar)
                iS = ismember(shockVar,varplot)==0;
                varplot = [varplot,shockVar(iS)];
            end
            iPlot = ismember(varendo,varplot);
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot,iPlot);
            end
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            chartName = [char(units0),'_DummyObs-v-',estimation,'_comparison_varplot'];
            legNames = {'Symmetric Prior','Asymmetric Prior'};
            IRFdisplay_asymmetry(IRFsulm_plot,varplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

            case 'varplot2'    %varplot only

            % Varplot
            if sum(ismember(shockVar,varplot)) < numel(shockVar)
                iS = ismember(shockVar,varplot)==0;
                varplot = [varplot,shockVar(iS)];
            end
            iPlot = ismember(varendo,varplot);
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot,iPlot);
            end
            varplot = varendo(iPlot);
            longnames = varendolong(iPlot);
            chartName = [char(units0),'_DummyObs-v-',estimation,'_comparison_varplot'];
            legNames = {'Symmetric Prior','Asymmetric Prior'};
            IRFdisplay_asymmetry(IRFsulm_plot,varplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = ismember(varendo,shockVar);
            iPlot2(posShock) = 1;        %reincluding policy variable
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot2,iPlot2);
            end
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            chartName = [char(units0),'_DummyObs-v-',estimation,'_comparison_nonvarplot'];
            IRFdisplay_asymmetry(IRFsulm_plot,nonvarplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

        end