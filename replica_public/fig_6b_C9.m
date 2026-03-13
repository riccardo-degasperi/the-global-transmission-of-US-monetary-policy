%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figures 6.b, C.9                                             %
%                                                                         %
% URL to the file:                                                        %
% Fig.6.b :  ./plots/fig_6b_C9/AA__median-median_overlaid_custom_1.pdf    %
% Fig.C.9 :  ./plots/fig_6b_C9/AA__median-median_overlaid_1.pdf           %
%                                                                         %
% last updated: 02/03/2026                                                %
%-------------------------------------------------------------------------%

clear all
close all
clc

warning('off','all')

addpath(genpath(pwd))  %path to child functions

%-------------------------------------------------------------------------%
% Select countries

% Advanced countries:
group = {
'AUSTRALIA';
'AUSTRIA';
'BELGIUM';
'CANADA';
'DENMARK';
'FINLAND';
'FRANCE';
'GERMANY';
'ITALY';
'JAPAN';
'NETHERLANDS';
'NORWAY';
'SPAIN';
'SWEDEN';
'UK'};

UU = size(group,1);

% Loop over countries
for z = 1:UU
    
    clearvars -except z UU group
    close all

    disp(['Bilateral VAR for ',group{z,1}])
    disp('------------------------------------------')
    %---------------------------------------------------------------------%
    % Type of model
    modelOpt.fullGibbs      = 0;
    modelOpt.dummies        = 0;                    %1:on; 0:off
    modelOpt.triangular     = 1;                    %1:on; 0:off

    % Model details
    modelOpt.p              = 12;                   %number of lags for endogenous variables
    modelOpt.pexo           = 0;                    %lags for exogenous variables (must be < p)
    modelOpt.h              = 24;                   %horizons
    modelOpt.constant       = 1;                    %0:no constant; 1:constant; 2:constant+trend
    modelOpt.sim            = 600;                  %number of simulations
    modelOpt.burnin         = 200;                  %burnin
    modelOpt.jump           = 1;                    %stores a draw every jump draws

    % Date cutoffs
    dataOpt.lT              = 'full';               %lower date cutoff (either 'full' or a date in text format)
    dataOpt.uT              = 'full';               %upper date cutoff (either 'full' or a date in text format)
    dataOpt.frequency       = 'monthly';            %frequency of data, can be 'daily', 'weekly', 'monthly', 'quarterly', 'yearly' (also displayed in the charts)

    % Series selection
    dataOpt.dataset         = 'data_public.xlsx';   %dataset name including extension
    dataOpt.shockVar        = {'GS1'};              %policy variable
    dataOpt.shockUnit       = {};                   %unit originating the shock (leave empty if policy variable is in sheet "Global")
    dataOpt.shockSign       = 'positive';           %either 'positive' or 'negative'
    dataOpt.shockSize       = 1;                    %magnitude of the shocks. Either '1sd' or a scalar (also negative)
    dataOpt.varendo         = {'IP','CPI','CCPI','SPRDUeom','XtoM','TVOL','NER','POL','IR','GBY','CBCFCI','CBCRA','CBCXFI','CBCFIHUSD','CBCEHUSD'};
    dataOpt.varendo_global  = {'USIP','USCPI','USCCPI','USSPRDUeom','USXtoM','USTVOL','USNER','USGBY','USCBCFCI','USCBCRA','USCBCXFI','USCBCFIHUSD','USCBCEHUSD','GZEBP','VIX','ROILP','RCRBPI','GEAI','GS1'};
    dataOpt.varplot         = dataOpt.varendo;
    dataOpt.varexo          = {};                   %exogenous variables
    dataOpt.units           = group(z,1);           %unit for this specific iteration
    dataOpt.group           = group;                %all units in the analysis

    % Identification
    modelOpt.identification = 'iv';                 %'cholesky', 'iv', 'sign', 'jarocinski-karadi', 'iv+sign'
    dataOpt.variv           = {'mar'};
    dataOpt.shockIV         = [1];                  %position of shocks to identify with porxies


    % Options for dogmatic priors (works with modelOpt.triangular)
    modelOpt.shutVar        = dataOpt.varendo;        %parameters whose prior variance is 0
    modelOpt.shutEq         = dataOpt.varendo_global; %equations on which to apply the extra shrinkage
    miscOpt.caseTriangular  = 'Chan';                 %either: SW:system-wide gibbs sampler; CCM:homoscesadastic Carriero, Clark, Marcellino (2019); Chan:homoscedastic Chan (2019); SV:stochastic volatility Carriero, Clark, Marcellino (2019)

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
    hp.l5                   = 100;                  %constant and exogenous
    hp.l6                   = 1;                    %extra shrinkage on cross-variable parameters (=1 is off)
    hp.l7                   = 0.05;                 %informativeness of prior for pandemic time dummies
    hp.eps                  = 1e-10;                %dogmatic shrinkage parameter

    % Options for Structural Scenario Analysis
    modelOpt.SSAfullIRF     = 1;
    modelOpt.SSAfullTS      = 0;
    miscOpt.SSAirf          = 0;                    %1:on; 2:off
    miscOpt.SSAts           = 0;                    %1:on; 2:off
    ssaOpt.SSAlT            = '';                   %(for SSAts) sets lower Xlim for scenarios plot
    ssaOpt.SSAscript        = 'SSArestrictions';    %name of function defining SSA restrictions
    ssaOpt.KLtresh          = 0.7;                  %threshold for KL divergence 
    ssaOpt.impulse          = dataOpt.shockVar;     %generate IRFs to shocks in position 'impulse' (you can set this to dataOpt.shockVar)
    ssaOpt.SSAy             = {{'ROILP','RCRBPI'};
                               {'NER','USNER'};
                               {'SPRDUeom','CBCFCI','CBCRA','USSPRDUeom','USCBCFCI','USCBCRA','GZEBP','VIX'};
                               {'POL'}};   %variables whose path is restricted
    ssaOpt.SSAe             = {{'ROILP','RCRBPI'};
                               {'NER','USNER'};
                               {'SPRDUeom','CBCFCI','CBCRA','USSPRDUeom','USCBCFCI','USCBCRA','GZEBP','VIX'};
                               {'POL'}};   %variables whose innovations are used to deliver the path
    ssaOpt.customLabels    = {'No commodity prices';'No exchange rates';'No financial variables';'No policy action';'Baseline'};

    % Charts options
    plotOpt.saveCharts      = 1;                    %1:save graphs to folder ./plots
    plotOpt.saveFig         = 0;
    plotOpt.plotTimeSeries  = 0;                    %1:plot time series; 0:don't
    plotOpt.homogeneity     = 0;                    %(for median group exercise) overlays unit median to group responses
    plotOpt.plotType        = 'singleShock';        %either 'singleShock' or 'allShocks'
    plotOpt.whichPlot       = 'none';               %either 'all','varendo','varplot','colour'
    plotOpt.cb90            = 1;                    %1:plots also 90% confidence bands; 0:plots only 68% confidence bands
    plotOpt.folder          = 'fig_6b_C9';          %either a folder name or []

    plotOpt.allSeries.plotRows    = 6;              %max number of rows per plot
    plotOpt.allSeries.plotCols    = 6;
    plotOpt.varplot.plotRows      = 3;              %max number of rows per plot
    plotOpt.varplot.plotCols      = 5;
    plotOpt.font               = 'Times New Roman'; 

    % Misc options
    miscOpt.dropExplosive   = 1;                    %drop explosive draws
    miscOpt.dateCheck       = 0;                    %check sample length when option 'full' is on. [Requires user input]
    miscOpt.FEVD            = 0;                    %Forecast Error Variance Decomposition. 0:off; 1:on
    miscOpt.HD              = 0;                    %Historical Decomposition. 0:off; 1:on
    miscOpt.interpolate     = 1;                    %1:interpolate NaNs; 0:don't
    miscOpt.interpolCase    = 'spline';             %'MA':centered MA; 'spline':cubic spline
    miscOpt.ws              = 1;                    %window size (for 'MA' option)
    miscOpt.saveSpace       = 0;                    %if =1, deletes the folder containing all IRFs at the end of the script
    miscOpt.plotAcross      = 1;                    %plot responses of specific variables across units (1:on/0:off)
    modelOpt.aggregation    = 'median-median';      %aggregation method (all, median-median, or mean-median)
    % all           : stacks all IRFs of all units and takes median and bands
    %                 out of that.
    % median-median : for each iteration and horizon, takes median across
    %                 units, then computes median and bands of this "median unit"
    % mean-median   : for each iteration and horizon, takes mean across
    %                 units, then computes median and bands of this "mean unit"

    %---------------------------------------------------------------------%
    % Run analysis for unit z
    main
    
    %---------------------------------------------------------------------%
    % Save IRFs for group analysis
    stoFolder = 'sto_AEs_CF_alt';                      %folder where to save IRF mat files
    if z==1; mkdir(stoFolder); end
    storeIRFs_fullGibbs
    
end


%-------------------------------------------------------------------------%
%% Generate group responses
clearvars -except stoFolder
close all
clc

% Load options
modelOpt = load([stoFolder,'/modelOpt']);
dataOpt  = load([stoFolder,'/dataOpt']);
plotOpt  = load([stoFolder,'/plotOpt']);
miscOpt  = load([stoFolder,'/miscOpt']);
ssaOpt   = load([stoFolder,'/ssaOpt']);

%plotOpt.ylab = dataOpt.ylab;

% Produce charts
plotOpt.whichPlot    = 'none';
plotOpt.customLabels = ssaOpt.customLabels;
plotOpt.font         = 'Times New Roman'; 

disp('--> Plotting group scenarios')
IRFssa = load([stoFolder,'/stoSSA']);
plotPrefix = '';
dataOpt.varChannel = ssaOpt.SSAy;
makeGroupResponses(IRFssa,'channels',plotPrefix,modelOpt,dataOpt,plotOpt);
clear IRFssa

%-------------------------------------------------------------------------%
% Delete folder storing the IRFs
if miscOpt.saveSpace
rmdir(stoFolder,'s');
end



