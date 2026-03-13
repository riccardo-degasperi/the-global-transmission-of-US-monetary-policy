function [] = displayComparisonGroups(IRFsulm,dataOpt,plotOpt)
        
        %-------------------------------------------------------------------------%
        disp('--> Plotting comparison charts: symmetric vs. asymmetric priors...')

        % Unpacking
        shockVar        = dataOpt.shockVar;
        varendo         = dataOpt.varendo;
        varendolong     = dataOpt.varendolong;
        varplot         = dataOpt.varplot;
        units           = dataOpt.units;

        whichPlot       = plotOpt.whichPlot;
        estimation      = plotOpt.estimation;
        legNames        = plotOpt.legNames;

        % Patch for chartNames when units are more than 1
        if numel(units) > 1
            units0 = join(string(units),'_');
        else
            units0 = units;
        end
        
        % Get struc names
        names = fieldnames(IRFsulm);

        switch whichPlot

            case 'all'        %plot all versions

            % Varendo
            chartName = [estimation,'_comparison_varendo'];
            IRFdisplay_asymmetry(IRFsulm,varendo,varendolong,shockVar,legNames,'allSeries',plotOpt,chartName);

            % Varplot
            if ~strcmp(varplot,shockVar)
                varplot = [varplot,shockVar];
            end
            iPlot = ismember(varendo,varplot);
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot,iPlot);
            end
            longnames = varendolong(iPlot);
            chartName = [estimation,'_comparison_varplot'];
            IRFdisplay_asymmetry(IRFsulm_plot,varplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = find(strcmp(varendo,shockVar));
            iPlot2(posShock) = 1;        %reincluding policy variable
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot2,iPlot2);
            end
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            chartName = [estimation,'_comparison_nonvarplot'];
            IRFdisplay_asymmetry(IRFsulm_plot,nonvarplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

            case 'varendo'    %varendo only

            chartName = [estimation,'_comparison_varendo'];
            IRFdisplay_asymmetry(IRFsulm,varendo,varendolong,shockVar,legNames,'allSeries',plotOpt,chartName)

            case 'varplot'    %varplot only

            % Varplot
            if ~strcmp(varplot,shockVar)
                varplot = [varplot,shockVar];
            end
            iPlot = ismember(varendo,varplot);
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot,iPlot);
            end
            longnames = varendolong(iPlot);
            chartName = [estimation,'_comparison_varplot'];
            IRFdisplay_asymmetry(IRFsulm_plot,varplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

            % Non varplot variables
            iPlot2 = ~iPlot;             %policy variable missing
            posShock = find(strcmp(varendo,shockVar));
            iPlot2(posShock) = 1;        %reincluding policy variable
            for i = 1:numel(names)
                IRFsulm_plot.(names{i}) = IRFsulm.(names{i})(iPlot2,iPlot2);
            end
            nonvarplot = varendo(iPlot2);
            longnames = varendolong(iPlot2);
            chartName = [estimation,'_comparison_nonvarplot'];
            IRFdisplay_asymmetry(IRFsulm_plot,nonvarplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

        end