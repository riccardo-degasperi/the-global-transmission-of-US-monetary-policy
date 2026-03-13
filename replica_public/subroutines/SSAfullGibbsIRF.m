function out = SSAfullGibbsIRF(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt,ssaOpt)
% SSAfullGibbsIRF sets up the loops over shocks and structural scenarios
% and produces the overlaid chart of the different scenarios. The main
% engine is SSAsampler_irf.


% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%
disp('------------------------------------------')
disp('Full Gibbs sampler for Structural Scenario Analysis')
disp(' ')

% Initialise options for loop
ssaOpt_ = ssaOpt;

% Default chart labels
if isfield(ssaOpt,'customLabels')
    if isempty(ssaOpt.customLabels)
    customLabels = cellstr([repmat('restriction',numel(ssaOpt.SSAy),1) num2str(vec(1:numel(ssaOpt.SSAy)))]);
    else
    customLabels = ssaOpt.customLabels;
    end
else
    customLabels = cellstr([repmat('restriction',numel(ssaOpt.SSAy),1) num2str(vec(1:numel(ssaOpt.SSAy)))]);
end

% Remove whitespace
customLabels_ = cellfun(@(x) strrep(x, ' ', '_'), customLabels, 'UniformOutput', false);
customLabels_ = cellfun(@(x) strrep(x, '.', ''), customLabels_, 'UniformOutput', false);
customLabels_ = cellfun(@(x) strrep(x, ',', ''), customLabels_, 'UniformOutput', false);
customLabels_ = cellfun(@(x) strrep(x, '-', '_'), customLabels_, 'UniformOutput', false);

% Loop over shocks
for k = 1:numel(ssaOpt.impulse)

    % Select shock
    impulse = ssaOpt.impulse(k);

    % Loop over scenarios
    for j = 1:numel(ssaOpt.SSAy)

        % Display current program
        disp(['    Shock in position: ',impulse{:},'; Scenario: ',customLabels{j}])

        % Define common chartname
        if modelOpt.triangular
            estimation = ['Triangular_',miscOpt.caseTriangular];
        elseif modelOpt.dummies
            estimation = 'DummyObs';
        end
        chartName0 = [char(dataOpt.units),'_',estimation,'_',customLabels_{j},'_SSA_',impulse{:}];

        % Gibbs sampler
        ssaOpt_.SSAy = ssaOpt.SSAy{j};
        ssaOpt_.SSAe = ssaOpt.SSAe{j};
        ssaOpt_.impulse = impulse;
        SSAout = SSAsampler_irf(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt,ssaOpt_,chartName0);

        % Store output
        SSA.(impulse{:}).(customLabels_{j}) = SSAout;

        % Prepare structure for overlaid plot
        chIRFsulm.(customLabels_{j}) = SSAout.IRFulm_cf;
        chIRFs.(customLabels_{j}) = SSAout.IRF_cf;
        if j == numel(ssaOpt.SSAy)
            chIRFsulm.baseline = SSAout.IRFulm_bl;
            chIRFs.baseline = SSAout.IRF_bl;
        end

    end

    % Overlaid scenarios plot
    if numel(ssaOpt.SSAy) > 1
        if numel(customLabels) == numel(ssaOpt.SSAy)
            customLabels = [customLabels;'baseline'];   %add label for baseline if not present
        end
        chartName0 = [char(dataOpt.units),'_',estimation,'_SSA_',impulse{:},'_all'];
        style = {'k-','k--','k-.','k:','ko','r-'};
        IRFdisplay_structure(chIRFsulm,dataOpt.varendo,dataOpt.varendolong,impulse,customLabels,style,'allSeries',plotOpt,chartName0)
    end

    % Pack output for group analysis
    SSAgroup.(impulse{:}).chIRFs = chIRFs;
    SSAgroup.(impulse{:}).chIRFsulm = chIRFsulm;

end

% Pack all
out.SSAgroup = SSAgroup;
out.SSA = SSA;

% Plot counteracting shocks

B0_ = SSAout.B0nnr;

% Generate structural impulse responses
IRFs = IRFsbuild(SSAout.betas,[],B0_,dims,dataOpt.shockSign);

% Get IRFs median, upper and lower bounds
IRFsulm = IRFbands(IRFs,dims.N,dims.h);

% Plot IRFs
plotOpt.whichPlot = 'varendo';
chartName = [char(dataOpt.units),'_',estimation];
BVARplots(IRFsulm,chartName,dataOpt,plotOpt);