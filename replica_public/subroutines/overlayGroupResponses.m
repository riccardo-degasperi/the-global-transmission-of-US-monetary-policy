function [IRFulm, IRFulm_plot] = overlayGroupResponses(in,makeGroupCase,plotPrefix,modelOpt,dataOpt,plotOpt)
% This code aggregates the IRFs for each member of group and plots them.

% INPUT:
% - in            : structure; elements are memeber of group and each
%                   element is either another structure (if makeGroupCase 
%                   is "asymmetries" or "channels") or a cell containing IRFs.
% - makeGroupCase : default is empty, or "asymmetries" to plot the overlaid
%                   asymmetry plot, or "channels" to plot all channels plots
%                   plus the ovelaid one.
% - plotPrefix    : string that ends up in the chart names.


% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%


% Unpacking

varendo        = dataOpt.varendo;
varendolong    = dataOpt.varendolong;
varendo_global = dataOpt.varendo_global;
varplot        = dataOpt.varplot;
varplot0       = dataOpt.varplot;
%varChannel     = dataOpt.varChannel;
group          = dataOpt.group;
shockVar       = dataOpt.shockVar;
aggregation    = modelOpt.aggregation;

% System dimensions
if strcmp(makeGroupCase,'asymmetries') || strcmp(makeGroupCase,'channels')
    
    f = fieldnames(in.(group{1}));   %names of channels
    C = numel(f);                    %number of channels
    U = numel(group);                %number of countries
    N = numel(varendo);              %number of endogenous variables

elseif isempty(makeGroupCase)
    
    U = numel(group);            %number of countries
    N = numel(varendo);          %number of endogenous variables
    C = 1;                       %because there are no channels
    varChannel = {};             %to avoid inclusion of zeroed-out variables
    
else
    error('ERROR: makeGroupCase can either be "asymmetries" or "channels". Default is empty.')
end

for kk = 1:C
    
    % Get data for channel kk
    for i = 1:U
        if C == 1
        ch.(group{i}) = in.(group{i});
        else
        ch.(group{i}) = in.(group{i}).(f{kk});
        end
    end
    
    %---------------------------------------------------------------------%
    % Aggregation

    % Check that policy variable is in varendo_global
    if ~strcmp(varendo_global,shockVar)
        varendo_global = [varendo_global,shockVar];
    end

    % Find variables for which to skip 2-step aggregation
    for i = 1:numel(varendo_global)
    iGlobal0(i,:) = strcmp(varendo,varendo_global{i}); %2-step aggregation only if iGlobal = 1
    end
    iGlobal = sum(iGlobal0,1) == 1;

    % Check that the policy variable is in varplot
    if ~strcmp(varplot,shockVar)
        varplot = [varplot,shockVar];
        varplot0 = [varplot0,shockVar];
    end

    % Check if the zeroed out variables are among the variables to plot, otherwise add them
    varplot1 = varplot;
    for j = 1:numel(varChannel)
        varChannel0 = varChannel{j};
        if iscell(varChannel0)

            idx0 = [];
            for i = 1:numel(varChannel0)
                idx0(i,:) = strcmp(varplot1,varChannel0{i});
            end
            idx = sum(idx0,1);
            pos = sum(idx0,2) == 1;

            if sum(idx) == numel(varChannel0)
                %do nothing 
            else
                varplot1 = [varplot1,varChannel0(~pos)]; %here they get added disorderly
            end

        else

            idx = strcmp(varplot1,varChannel0);

            if sum(idx) == 1
                %do nothing 
            else
                varplot1 = [varplot1,varChannel0];       %here they get added disorderly
            end

        end

    end

    for i = 1:numel(varplot1)
        iPlot0(i,:) = strcmp(varendo,varplot1{i});
    end
    iPlot = sum(iPlot0,1) == 1;                          %here they are back in the right order
    varplot = varendo(iPlot);
    
    % This is in case you want to use the original varplot
    for i = 1:numel(varplot0)
        iPlot00(i,:) = strcmp(varendo,varplot0{i});
    end
    iPlot0 = sum(iPlot00,1) == 1;   
    
    %---------------------------------------------------------------------%
    switch  aggregation

        case 'all' 

        % Prepare IRFs
        for ii = 1:N
            for jj = 1:N
                tmp3 = [];                     %initialise container
                for i = 1:U
                    tmp1 = ch.(group{i});      %open IRF of unit i
                    tmp2 = tmp1{ii,jj};        %get responses of variable ii on jj
                    tmp3 = [tmp3 tmp2];        %stack it with responses for previous units 
                end
                IRFall{ii,jj} = tmp3;
            end
        end

        % Get dimensions
        N = size(IRFall,2);
        h = size(IRFall{1,1},1)-1;

        % Get IRFall median, upper and lower bounds
        IRFulm = IRFbands(IRFall,N,h);
        IRFulm_plot = IRFulm(iPlot,iPlot);
        IRFulm_plot0 = IRFulm(iPlot0,iPlot0);
        if C > 1
        chIRFulmStruc.(f{kk}) = IRFulm;
        chIRFulmPlot.(f{kk}) = IRFulm_plot;
        chIRFulmPlot0.(f{kk}) = IRFulm_plot0;
        end

        case 'median-median'

        % Prepare IRFs
        for ii = 1:N
            for jj = 1:N
                tmp3 = [];
                
                if iGlobal(ii) == 0
                
                for i = 1:U
                    tmp1 = ch.(group{i});      %open IRF of unit i
                    tmp2 = tmp1{ii,jj};        %get responses of variable ii on jj
                    tmp3 = cat(3,tmp3,tmp2);   %stack it with responses for previous units 
                end 
                tmp4 = median(tmp3,3);
                IRFall{ii,jj} = tmp4;
                
                elseif iGlobal(ii) == 1
                    
                for i = 1:U
                    tmp1 = ch.(group{i});      %open IRF of unit i
                    tmp2 = tmp1{ii,jj};        %get responses of variable ii on jj
                    tmp3 = [tmp3 tmp2];        %stack it with responses for previous units 
                end 
                IRFall{ii,jj} = tmp3;
                
                else
                    error('ERROR: aggregation median-median failed.')
                end

            end
        end

        % Get dimensions
        N = size(IRFall,2);
        h = size(IRFall{1,1},1)-1;

        % Get IRFall median, upper and lower bounds
        IRFulm = IRFbands(IRFall,N,h);
        IRFulm_plot = IRFulm(iPlot,iPlot);
        IRFulm_plot0 = IRFulm(iPlot0,iPlot0);
        if C > 1
        chIRFulmStruc.(f{kk}) = IRFulm;
        chIRFulmPlot.(f{kk}) = IRFulm_plot;
        chIRFulmPlot0.(f{kk}) = IRFulm_plot0;
        end
        
        case 'mean-median'

        % Prepare IRFs
        for ii = 1:N
            for jj = 1:N
                tmp3 = [];
                
                if iGlobal(ii) == 0
                
                for i = 1:U
                    tmp1 = ch.(group{i});      %open IRF of unit i
                    tmp2 = tmp1{ii,jj};        %get responses of variable ii on jj
                    tmp3 = cat(3,tmp3,tmp2);   %stack it with responses for previous units 
                end 
                tmp4 = mean(tmp3,3);
                IRFall{ii,jj} = tmp4;
                
                elseif iGlobal(ii) == 1
                    
                for i = 1:U
                    tmp1 = ch.(group{i});      %open IRF of unit i
                    tmp2 = tmp1{ii,jj};        %get responses of variable ii on jj
                    tmp3 = [tmp3 tmp2];        %stack it with responses for previous units 
                end 
                IRFall{ii,jj} = tmp3;
                
                else
                    error('ERROR: aggregation mean-median failed.')
                end

            end
        end

        % Get dimensions
        N = size(IRFall,2);
        h = size(IRFall{1,1},1)-1;

        % Get IRFall median, upper and lower bounds
        IRFulm = IRFbands(IRFall,N,h);
        IRFulm_plot = IRFulm(iPlot,iPlot);
        IRFulm_plot0 = IRFulm(iPlot0,iPlot0);
        if C > 1
        chIRFulmStruc.(f{kk}) = IRFulm;
        chIRFulmPlot.(f{kk}) = IRFulm_plot;
        chIRFulmPlot0.(f{kk}) = IRFulm_plot0;
        end
        
    end

end

%-------------------------------------------------------------------------%
if strcmp(makeGroupCase,'channels')
    
    % Overlaid channels plot
    legNames = [];   %legend
    chartName = ['AA_',plotPrefix,'_',aggregation,'_overlaid'];
    style = {'k-','k--','k-.','k:','ko','r-'};
    IRFdisplay_structure(chIRFulmStruc,varendo,varendolong,shockVar,legNames,style,'allSeries',plotOpt,chartName)

end

%-------------------------------------------------------------------------%
if strcmp(makeGroupCase,'asymmetries')
    
    % Overlaid asymmetries plot
    legNames = {'Tightening','Loosening'};   %legend
    chartName = ['AA_',plotPrefix,'_',aggregation,'_asymmetry'];
    IRFdisplay_asymmetry(chIRFulmStruc,varendo,varendolong,shockVar,legNames,'allSeries',plotOpt,chartName)
    
    chartName = ['AA_',plotPrefix,'_',aggregation,'_asymmetry_varplot'];
    longnames = varendolong(iPlot0);
    IRFdisplay_asymmetry(chIRFulmPlot0,varplot0,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

    
end



