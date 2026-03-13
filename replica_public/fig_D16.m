%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figure D.16                                                  %
%                                                                         %
% URL to the file:                                                        %
% Fig.D.16:  ./plots/fig_D16/Global_Triangular_Chan_SSA_GS1_custom_bands_1.pdf
%                                                                         %
% last updated: 03/03/2026                                                %
%-------------------------------------------------------------------------%

clear all
close all
clc

addpath(genpath(pwd))  %path to child functions

%-------------------------------------------------------------------------%
% Type of model
modelOpt.fullGibbs      = 0;
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
dataOpt.lT              = 'full';               %lower date cutoff (either 'full' or a date in text format)
dataOpt.uT              = 'full';               %upper date cutoff (either 'full' or a date in text format)
dataOpt.frequency       = 'monthly';            %frequency of data, can be 'daily', 'weekly', 'monthly', 'quarterly', 'yearly' (also displayed in the charts)

% Series selection
dataOpt.dataset         = 'data_public.xlsx';   %dataset name including extension
dataOpt.shockVar        = {'GS1','ROILP'};      %shock of interest (position in the impact matrix)
dataOpt.shockUnit       = {};                   %unit originating the shock (leave empty if policy variable is in sheet "Global")
dataOpt.shockSign       = 'positive';           %either 'positive' or 'negative'
dataOpt.shockSize       = 1;                    %magnitude of the shocks. Either '1sd' or a scalar (also negative)
dataOpt.varendo         = {'OECDIP','OECDCPI','OECDSPRDUeom','EIA1955','OECDSTOCKS2','DIRADV','GEAI'};
dataOpt.varendo_global  = {'USIP','USCPI','USSPRDUeom','USXtoM','USTVOL','USNER','USGBY','GZEBP','VIX','ROILP','GS1'};
dataOpt.varplot         = {'OECDIP','OECDCPI','USIP','USCPI','EIA1955','OECDSTOCKS2','ROILP','GS1'};
dataOpt.varexo          = {};                   %exogenous variables
dataOpt.units           = {'Global'};           %units (can be more than one)

% Identification
modelOpt.identification = 'iv';                 %'cholesky', 'iv', 'sign', 'jarocinski-karadi', 'iv+sign'
dataOpt.variv           = {'mar_extended','bh19'};
dataOpt.shockIV         = [1 2];                %position of shocks to identify with porxies

% Options for dogmatic priors (works with triangular)
modelOpt.shutVar        = {'OECDIP','OECDCPI','OECDSPRDUeom','DIRADV'};   %parameters whose prior variance is 0
modelOpt.shutEq         = dataOpt.varendo_global; %equations on which to apply the extra shrinkage
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
modelOpt.SSAfullIRF     = 1;
modelOpt.SSAfullTS      = 0;
miscOpt.SSAirf          = 0;                    %1:on; 2:off
miscOpt.SSAts           = 0;                    %1:on; 2:off
ssaOpt.SSAlT            = '';                   %(for SSAts) sets lower Xlim for scenarios plot
ssaOpt.SSAscript        = 'SSArestrictions_fixImpulse';  %name of function defining SSA restrictions (customise it before launch)
ssaOpt.KLtresh          = 0.6;                  %threshold for KL divergence
ssaOpt.impulse          = {'GS1'};              %generate IRFs to shocks in position 'impulse' (you can set this to dataOpt.shockVar)
ssaOpt.SSAy             = {{'ROILP','GS1'}};    %variables whose path is restricted (can define many in this format SSAy = {{'var1','var2'};{'var1'};{'var1','var2','var3'}};)
ssaOpt.SSAe             = {{'ROILP'}};          %variables whose innovations are used to deliver the path
ssaOpt.customLabels     = {'No oil price';'Baseline'};

% Charts options
plotOpt.saveCharts      = 1;                    %1:save graphs to folder ./plots
plotOpt.saveFig         = 0;
plotOpt.plotTimeSeries  = 0;                    %1:plot time series; 0:don't
plotOpt.plotType        = 'singleShock';        %either 'singleShock' or 'allShocks'
plotOpt.whichPlot       = 'all';                %either 'all','varendo','varplot','colour'
plotOpt.cb90            = 1;                    %1:plots also 90% confidence bands; 0:plots only 68% confidence bands
plotOpt.folder          = 'fig_D16';            %either a folder name or []

plotOpt.allSeries.plotRows = 4;                 %max number of rows per plot
plotOpt.allSeries.plotCols = 5;
plotOpt.varplot.plotRows   = 2;                 %max number of rows per plot
plotOpt.varplot.plotCols   = 4;
plotOpt.font               = 'Times New Roman'; %font for charts

% Misc options
miscOpt.asymmetries     = 0;                    %explore asymmetric responses (only for IV identification)
miscOpt.dropExplosive   = 1;                    %drop explosive draws
miscOpt.dateCheck       = 0;                    %check sample length when option 'full' is on. [Requires user input]
miscOpt.firstDifference = 0;                    %takes first differences according to column FD in sheet 'Trans'
miscOpt.BVARcase        = 'BGR';                %either 'GLP' or 'BGR' (GLP is slower with larger systems)
miscOpt.FEVD            = 0;                    %Forecast Error Variance Decomposition. 0:off; 1:on
miscOpt.HD              = 0;                    %Historical Decomposition. 0:off; 1:on
miscOpt.interpolate     = 1;                    %1:interpolate NaNs; 0:don't
miscOpt.interpolCase    = 'spline';             %'MA':centered MA; 'spline':cubic spline
miscOpt.ws              = 1;                    %window size (for 'MA' option)


%-------------------------------------------------------------------------%
% Run analysis
main


% Display priors and save options to a text file
filename = 'SSAfullGibbs_triangular';
out_ = SSAout.SSA.GS1.No_oil_price;
printDetails(out_,filename,data,hp,modelOpt,dataOpt,miscOpt,plotOpt)


%% Custom plot
impulse = ssaOpt.impulse{:};
if modelOpt.triangular
    estimation = ['Triangular_',miscOpt.caseTriangular];
elseif modelOpt.dummies
    estimation = 'DummyObs';
end
customLabels0 = ssaOpt.customLabels;
customLabels = customLabels0([2;1]);
chIRFsulm0 = SSAout.SSAgroup.(impulse).chIRFsulm;
chIRFsulm = orderfields(chIRFsulm0,[2;1]);
chartName0 = [char(dataOpt.units),'_',estimation,'_SSA_',impulse,'_all'];
style = {'k-','k--','k-.','k:','ko','r-'};
IRFdisplay_SSA(chIRFsulm,dataOpt.varendo,dataOpt.varendolong,impulse,customLabels,'allSeries',plotOpt,chartName0)

% Custom plot
iPlot = ismember(dataOpt.varendo,dataOpt.varplot);
rest = fields(chIRFsulm);
for j = 1:numel(rest)
    chIRFsulm2.(rest{j}) = chIRFsulm.(rest{j})(iPlot,iPlot);
end
chartName0 = [char(dataOpt.units),'_',estimation,'_SSA_',impulse,'_custom'];
plotOpt0 = plotOpt;
plotOpt0.allSeries.plotRows = 2;
plotOpt0.allSeries.plotCols = 4;
IRFdisplay_SSA(chIRFsulm2,dataOpt.varendo(iPlot),dataOpt.varendolong(iPlot),impulse,customLabels,'allSeries',plotOpt0,chartName0)

chartName0 = [char(dataOpt.units),'_',estimation,'_SSA_',impulse,'_custom_bands'];
IRFdisplay_SSA_bands(chIRFsulm2,dataOpt.varendo(iPlot),dataOpt.varendolong(iPlot),impulse,customLabels,'allSeries',plotOpt0,chartName0)

