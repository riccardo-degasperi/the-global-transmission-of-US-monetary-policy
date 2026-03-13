%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figure D.6.b                                                 %
%                                                                         %
% URL to the file:                                                        %
% Fig.D.6.b :  ./plots/fig_D6/overlay/EM_JK_comparison.pdf                %
%                                                                         %
% last updated: 03/03/2026                                                %
%-------------------------------------------------------------------------%

clear
close all
clc

addpath(genpath(pwd))  %path to child functions

% Get results for emerging economies
run('runEM_JK.m')

%% Load options
AEmodelOpt = load('sto_EMs/modelOpt');
AEdataOpt  = load('sto_EMs/dataOpt');
AEplotOpt  = load('sto_EMs/plotOpt');
AEmiscOpt  = load('sto_EMs/miscOpt');

EMmodelOpt = load('sto_EMs_JK/modelOpt');
EMdataOpt  = load('sto_EMs_JK/dataOpt');
EMplotOpt  = load('sto_EMs_JK/plotOpt');
EMmiscOpt  = load('sto_EMs_JK/miscOpt');

% Load IRFs
load('sto_EMs/sto_groupBaseline_BVAR_t');
AEIRF = IRFmedian;
load('sto_EMs_JK/sto_groupBaseline_BVAR_t');
EMIRF = IRFmedian;

% Get varplot indices
if sum(ismember(AEdataOpt.shockVar,AEdataOpt.varplot)) < numel(AEdataOpt.shockVar)
    iS = ismember(AEdataOpt.shockVar,AEdataOpt.varplot)==0;
    AEdataOpt.varplot = [AEdataOpt.varplot,AEdataOpt.shockVar(iS)];
end
iAE = ismember(AEdataOpt.varendo,AEdataOpt.varplot);
if sum(ismember(EMdataOpt.shockVar,EMdataOpt.varplot)) < numel(EMdataOpt.shockVar)
    iS = ismember(EMdataOpt.shockVar,EMdataOpt.varplot)==0;
    EMdataOpt.varplot = [EMdataOpt.varplot,EMdataOpt.shockVar(iS)];
end
iEM = ismember(EMdataOpt.varendo,EMdataOpt.varplot);

% Index useful to skip plot of variables available only for AE
iskip = ismember(AEdataOpt.varplot,EMdataOpt.varplot);

% Pack IRF for plot
stoIRF.EMIRF = EMIRF(iAE,iAE);
stoIRF.AEIRF = AEIRF(iAE,iAE);

% Plot figure
AEplotOpt.folder = 'fig_D6/overlay/';
legNames = {'Jarocinski-Karadi','Baseline'};
chartName = 'EM_JK_comparison';
longnames = AEdataOpt.varendolong(iAE);
plotOpt.ylab = AEdataOpt.ylab(iAE);
IRFdisplay_overlay(stoIRF,AEdataOpt.varplot,iskip,longnames,AEdataOpt.shockVar,legNames,'varplot',AEplotOpt,chartName)

