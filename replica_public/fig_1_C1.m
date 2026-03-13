%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figures 1, C.1                                               %
%                                                                         %
% URL to the file:                                                        %
% Fig.1   :  ./plots/fig_1_C1/overlay_global_1.pdf                        %
% Fig.C.1 :  ./plots/fig_1_C1/Global_FullGibbs_nonvarplot_1.pdf           %
%                                                                         %
% last updated: 02/03/2026                                                %
%-------------------------------------------------------------------------%

clear all
close all
clc

warning('off','all')

addpath(genpath(pwd))  %path to child functions

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
dataOpt.uT              = '2018m12';            %upper date cutoff (either 'full' or a date in text format)
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
modelOpt.identification = 'iv';                 %'cholesky' or 'iv'
dataOpt.variv           = {'mar'};
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
dataOpt.lT_covid        = '';                   %first time dummy for pandemic priors
dataOpt.uT_covid        = '';                   %last time dummy for pandemic priors

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
plotOpt.folder          = 'fig_1_C1';           %either a folder name or []

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

%% Store output
sto.out.IRFsulm = out.IRFsulm;
sto.out.shock   = out.shock;
sto.plotOpt  = plotOpt;
sto.dataOpt  = dataOpt;
sto.modelOpt = modelOpt;
sto.miscOpt  = miscOpt;
sto.ssaOpt   = ssaOpt;
sto.dims     = dims;
%save([pwd '/plots/',plotOpt.folder,'/output'],'-struct','sto')

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

longnames(1:10) = {'OECD & US Production','OECD & US CPI','OECD ex. NA & US Stock Price',...
      'Global & US Financial Cond. Index','Global & US Risk Appetite','Global & US Fixed Income Holdings',...
      'Global & US Equity Holdings','EMs & US Inflow','EMs & US Outflow','Interest Rate Differential'};

plotOpt.ylab = dataOpt.ylab(iGLO);
IRFdisplay_overlay_Global(stoIRF,varplot1,iskip,longnames,dataOpt.shockVar,legNames,'varplot',plotOpt,chartName)

