%-------------------------------------------------------------------------%
% Checks and patches
%
% last modified: 26/03/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Check presence of exogenous for SSA
if (modelOpt.SSAfullIRF || modelOpt.SSAfullTS || miscOpt.SSAirf || miscOpt.SSAts) && ( modelOpt.constant ~= 1 || ~isempty(dataOpt.varexo) )
    error('ERROR: SSA algorithm works only with constant and covid priors. No other exogenous variables is implemented.')
end

% Check selected SSA routines
if modelOpt.SSAfullIRF && modelOpt.SSAfullTS
    error('ERROR: Run one exercise at the time between in-sample and out-of-sample structural scenario analysis.')
end

if ~modelOpt.SSAfullIRF && ~modelOpt.SSAfullTS && miscOpt.SSAirf && miscOpt.SSAts
    error('ERROR: Run one exercise at the time between in-sample and out-of-sample structural scenario analysis.')
end

% Check that dates are in the correct format
if strcmp(dataOpt.frequency,'monthly')
    if max(strlength(data.dates)) ~= 7 || min(strlength(data.dates)) ~= 6
        error('ERROR: monthly dates have to be defined as 2023m1, 1975m12, or 2023M1, 1975M12 both in the Excel source file and in the options.')
    end
elseif strcmp(dataOpt.frequency,'quarterly')
    if any(strlength(data.dates) ~= 6)
        error('ERROR: quarterly dates have to be defined as 2023q1, 1975q4, or 2023Q1, 1975Q4 both in the Excel source file and in the options.')
    end
elseif strcmp(dataOpt.frequency,'yearly')
    if any(strlength(data.dates) ~= 4)
        error('ERROR: yearly dates have to be defined as 2023 or 1975 both in the Excel source file and in the options.')
    end
elseif strcmp(dataOpt.frequency,'daily') || strcmp(dataOpt.frequency,'weekly')
    if any(strlength(data.dates) ~= 10)
        error('ERROR: daily or weekly dates have to be defined as 24/12/1990 or 03/01/2023 both in the Excel source file and in the options.')
    end
else
    error('ERROR: the selected frequency is not supported')
end

% Add frequency to plotOpt
plotOpt.frequency = dataOpt.frequency;

% Patch for shockVar as cell
if ischar(dataOpt.shockVar)
    dataOpt.shockVar = cellstr(dataOpt.shockVar);
end

% Check that sim-burnin is divisible by jump
if mod(modelOpt.sim-modelOpt.burnin,modelOpt.jump) ~= 0
    error('ERROR: sim-burnin has to be divisible by jump. Change it in the settings.')
end

% Drop explosives if you selected HD or FEVD
if (miscOpt.FEVD == 1 || miscOpt.HD == 1) && miscOpt.dropExplosive == 0
    warning('WARNING: dropExplosive has been set to 1 as you asked for FEVD or HD.')
    miscOpt.dropExplosive = 1;
end

% SSA only works with constant = 1 and without exogenous variables
if (miscOpt.SSAirf || miscOpt.SSAts || modelOpt.SSAfullIRF || modelOpt.SSAfullTS) && (modelOpt.constant ~= 1 || ~isempty(dataOpt.varexo))
    error('ERROR: Structural Scenario Analysis is only implemented for modelOpt.constant = 1.')
end