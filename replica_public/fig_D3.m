%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figure D.3                                                   %
%                                                                         %
% URL to the file:                                                        %
% Fig.D.3 :  ./plots/fig_D3/Global_BVAR_multishocks_varplot_1.pdf         %
%                                                                         %
% last updated: 03/03/2026                                                %
%-------------------------------------------------------------------------%

clear all
close all
clc

warning('off','all')

addpath(genpath(pwd))  %path to child functions


shocks = {'mar_extended';'mar'};
UU = size(shocks,1);

stoBVAR = [];
stoBVARchannels = [];

% Loop over shocks
for z = 1:UU
    
    clearvars -except z UU shocks stoBVAR stoBVARchannels
    close all
    clc
    
    
    %-------------------------------------------------------------------------%
    % Type of model
    modelOpt.fullGibbs      = 1;
    modelOpt.dummies        = 0;                    %1:on; 0:off
    modelOpt.triangular     = 1;                    %1:on; 0:off
    
    % Model details
    modelOpt.p              = 12;                   %number of lags for endogenous variables
    modelOpt.pexo           = 0;                    %lags for exogenous variables (must be < p)
    modelOpt.h              = 24;                   %horizons
    modelOpt.constant       = 1;                    %0:no constant; 1:constant; 2:constant+trend
    modelOpt.sim            = 1200;                 %number of simulations
    modelOpt.burnin         = 200;                  %burnin
    modelOpt.jump           = 1;                    %stores a draw every jump draws
    
    % Date cutoffs
    dataOpt.lT              = '1990m1';             %lower date cutoff (either 'full' or a date in text format)
    dataOpt.uT              = 'full';               %upper date cutoff (either 'full' or a date in text format)
    dataOpt.frequency       = 'monthly';            %frequency of data, can be 'daily', 'weekly', 'monthly', 'quarterly', 'yearly' (also displayed in the charts)
    
    
    % Series selection
    dataOpt.dataset         = 'data_public.xlsx';   %dataset name including extension
    dataOpt.shockVar        = {'GS1'};              %shock of interest (position in the impact matrix)
    dataOpt.shockUnit       = {};                   %unit originating the shock (leave empty if policy variable is in sheet "Global")
    dataOpt.shockSign       = 'positive';           %either 'positive' or 'negative'
    dataOpt.shockSize       = 1;                    %magnitude of the shocks. Either '1sd' or a scalar (also negative)
    dataOpt.varendo         = {'OECDIP','OECDCPI','OECDSPRDUeom','GCBCFCI','GCBCRA','GCBCFIHUSD','GCBCEHUSD','EMINFLOW','EMOUTFLOW','DIRADV','EURUSD','GBPUSD','JPYUSD','GEAI'};
    dataOpt.varendo_global  = {'USIP','USCPI','USSPRDUeom','USXtoM','USTVOL','USNER','USGBY','USCBCFCI','USCBCRA','USCBCFIHUSD','USCBCEHUSD','USINFLOW','USOUTFLOW','GZEBP','VIX','RCRBPI','ROILP','GS1'};
    dataOpt.varplot         = {'OECDIP','OECDCPI','OECDSPRDUeom','GCBCFCI','GCBCRA','GCBCFIHUSD','GCBCEHUSD','EMINFLOW','EMOUTFLOW','DIRADV','EURUSD','GBPUSD','JPYUSD','RCRBPI','ROILP','GEAI'};   %selected variables to plot
    dataOpt.varexo          = {};                   %exogenous variables
    dataOpt.units           = {'Global'};           %units (can be more than one)
    
    % Identification
    modelOpt.identification = 'iv';                 %'cholesky', 'iv', 'sign', 'jarocinski-karadi', 'iv+sign'
    dataOpt.variv           = shocks(z,1);
    dataOpt.shockIV         = [];                   %position of shocks in shockVar to be identified with proxies
    
    % Options for dogmatic priors (works with triangular)
    modelOpt.shutVar      = {'OECDIP','OECDCPI','OECDSPRDUeom','GCBCFCI','GCBCRA','GCBCFIHUSD','GCBCEHUSD','EMINFLOW','EMOUTFLOW','DIRADV','EURUSD','GBPUSD','JPYUSD'};   %parameters whose prior variance is 0
    modelOpt.shutEq       = dataOpt.varendo_global; %equations on which to apply the extra shrinkage
    miscOpt.caseTriangular  = 'Chan';               %either: SW:system-wide gibbs sampler; CCM:homoscesadastic Carriero, Clark, Marcellino (2019); Chan:homoscedastic Chan (2019); SV:stochastic volatility Carriero, Clark, Marcellino (2019)
    
    % Options for BVAR
    hp.GLP                  = 1;                    %optimal prior selection
    hp.hyperpriors          = 1;                    %(if GLP=1): 1:uses hyperpriors; 0:ML only
    hp.niw                  = 1;                    %Normal Inverse Wishart prior
    hp.soc                  = 0;                    %Sum of Coefficients prior
    hp.cop                  = 0;                    %Copersistence prior
    hp.lag                  = 0;                    %optimise lag-decay prior
    hp.std                  = 0;                    %optimise diagonal elements of the scale matrix of the IW prior on the residual variance
    hp.cross                = 0;                    %optimise extra shrinkage parameter on cross-variable coefficients
    hp.covid                = 0;                    %pandemic priors (Cascaldi-Garcia, 2024)
    
    % Setting for pandemic priors time dummies
    dataOpt.lT_covid        = '2020m3';             %first time dummy for pandemic priors
    dataOpt.uT_covid        = '2020m8';             %last time dummy for pandemic priors
    
    % Initial values for priors
    hp.l1                   = 0.2;                  %Normal-Inverse-Wishart
    hp.l2                   = 1;                    %lag decay (=1 is off)
    hp.l3                   = 1;                    %Sum-of-coefficients (=1 is off)
    hp.l4                   = 1;                    %Co-persistence (=1 is off)
    hp.l5                   = 1000;                 %constant and exogenous
    hp.l6                   = 1;                    %extra shrinkage on cross-variable parameters (=1 is off)
    hp.l7                   = 0.05;                 %informativeness of prior for pandemic time dummies
    hp.eps                  = 1e-10;                %dogmatic shrinkage parameter
    
    % Options for Structural Scenario Analysis
    modelOpt.SSAfullIRF     = 0;
    modelOpt.SSAfullTS      = 0;
    miscOpt.SSAirf          = 0;                    %1:on; 2:off
    miscOpt.SSAts           = 0;                    %1:on; 2:off
    ssaOpt.SSAlT            = '';                   %(for SSAts) sets lower Xlim for scenarios plot
    ssaOpt.SSAscript        = '';                   %name of function defining SSA restrictions (customise it before launch)
    ssaOpt.KLtresh          = 0.6;                  %threshold for KL divergence
    ssaOpt.impulse          = {};                   %generate IRFs to shocks in position 'impulse' (you can set this to dataOpt.shockVar)
    ssaOpt.SSAy             = {};                   %variables whose path is restricted (can define many in this format SSAy = {{'var1','var2'};{'var1'};{'var1','var2','var3'}};)
    ssaOpt.SSAe             = {};                   %variables whose innovations are used to deliver the path
    ssaOpt.customLabels     = {};
    
    % Charts options
    plotOpt.saveCharts      = 1;                    %1:save graphs to folder ./plots
    plotOpt.saveFig         = 0;
    plotOpt.plotTimeSeries  = 0;                    %1:plot time series; 0:don't
    plotOpt.plotType        = 'singleShock';        %either 'singleShock' or 'allShocks'
    plotOpt.whichPlot       = 'varplot2';           %either 'all','varendo','varplot','colour'
    plotOpt.cb90            = 1;                    %1:plots also 90% confidence bands; 0:plots only 68% confidence bands
    plotOpt.folder          = ['fig_D3/',shocks{z}]; %either a folder name or []
    
    plotOpt.allSeries.plotRows = 6;                 %max number of rows per plot
    plotOpt.allSeries.plotCols = 6;
    plotOpt.varplot.plotRows   = 4;                 %max number of rows per plot
    plotOpt.varplot.plotCols   = 4;
    plotOpt.font               = 'Times New Roman'; %font for charts
    
    % Misc options
    miscOpt.dropExplosive   = 1;                    %drop explosive draws
    miscOpt.dateCheck       = 0;                    %check sample length when option 'full' is on. [Requires user input]
    miscOpt.FEVD            = 0;                    %Forecast Error Variance Decomposition. 0:off; 1:on
    miscOpt.HD              = 0;                    %Historical Decomposition. 0:off; 1:on
    miscOpt.interpolate     = 1;                    %1:interpolate NaNs; 0:don't
    miscOpt.interpolCase    = 'spline';             %'MA':centered MA; 'spline':cubic spline
    miscOpt.ws              = 1;                    %window size (for 'MA' option)
    
    
    %-------------------------------------------------------------------------%
    % Run analysis
    main

    %-------------------------------------------------------------------------%
    % Get varplot indices
    if sum(ismember(dataOpt.shockVar,dataOpt.varplot)) < numel(dataOpt.shockVar)
        iS = ismember(dataOpt.shockVar,dataOpt.varplot)==0;
        varplot1 = [dataOpt.varplot,dataOpt.shockVar(iS)];
    end
    iGLO = ismember(dataOpt.varendo,varplot1);
    iUS = ~ismember(dataOpt.varendo,dataOpt.varplot);
    
    % Selector for variables in US block that have an equivalent in Global one
    iUSred = logical([1 1 1 0 0 0 0 1 1 1 1 1 1 0 0 1]);
    
    % Selector for variables in Global block that have no US equivalent
    iskip = logical([1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1]);
    
    % Pack IRF for plot
    stoIRF.GLO = out.IRFsulm(iGLO,iGLO);
    tmp_US = out.IRFsulm(iUS,iUS);
    stoIRF.US = tmp_US(iUSred,iUSred);
    
    % Plot figure
    legNames = {'Global','US'};
    chartName = 'overlay_global';
    longnames = dataOpt.varendolong(iGLO);
    plotOpt.ylab = dataOpt.ylab(iGLO);
    IRFdisplay_overlay(stoIRF,varplot1,iskip,longnames,dataOpt.shockVar,legNames,'varplot',plotOpt,chartName)

    %-------------------------------------------------------------------------%
    % Save results
    stoBVAR.(shocks{z}) = out.IRFsulm;

end

%-------------------------------------------------------------------------%
%% Stack and plot

% Unpacking
plotType       = plotOpt.plotType;
whichPlot      = plotOpt.whichPlot;
varendo        = dataOpt.varendo;
varendolong    = dataOpt.varendolong;
varplot        = dataOpt.varplot;
shockVar       = dataOpt.shockVar;
plotOpt.ylab   = dataOpt.ylab;
ylab0          = dataOpt.ylab;
plotOpt.folder = 'fig_D3/';

% Varplot
if sum(ismember(shockVar,varplot)) < numel(shockVar)
    iS = ismember(shockVar,varplot)==0;
    varplot = [varplot,shockVar(iS)];
end
iPlot = ismember(varendo,varplot);
for z = 1:numel(shocks)
stoBVAR_varplot.(shocks{z}) = stoBVAR.(shocks{z})(iPlot,iPlot);
end
varplot = varendo(iPlot);
longnames = varendolong(iPlot);
plotOpt.ylab = ylab0(iPlot);
legNames = {'MAR extended','baseline'};
chartName = 'Global_BVAR_multishocks_varplot';
IRFdisplay_stack(stoBVAR_varplot,varplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)

% Non varplot variables
iPlot2 = ~iPlot;             %policy variable missing
posShock = ismember(varendo,shockVar);
iPlot2(posShock) = 1;        %reincluding policy variable
for z = 1:numel(shocks)
stoBVAR_nonvarplot.(shocks{z}) = stoBVAR.(shocks{z})(iPlot2,iPlot2);
end
nonvarplot = varendo(iPlot2);
longnames = varendolong(iPlot2);
plotOpt.ylab = ylab0(iPlot2);
legNames = {'MAR extended','baseline'};
chartName = 'Global_BVAR_multishocks_nonvarplot';
IRFdisplay_stack(stoBVAR_nonvarplot,nonvarplot,longnames,shockVar,legNames,'varplot',plotOpt,chartName)
        