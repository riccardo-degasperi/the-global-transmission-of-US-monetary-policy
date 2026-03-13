function out = SSArestrictions(impulse,dims,dataOpt,ssaOpt)
% SSArestrictions defines the restrictions on the path of observables and
% on the innovations/structural shocks of the VAR that deliver the chosen
% path for the observables.
% Adapted from Antolin-Diaz, Petrella, Rubio-Ramirez (2021)
%
% Please customise the 'PART TO CUSTOMISE'
%
% Last updated: 28/02/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpack
h        = dims.h+1;           %+1 for consistency with IRFbuild function
N        = dims.N;
varendo  = dataOpt.varendo;
SSAy     = ssaOpt.SSAy;
SSAe     = ssaOpt.SSAe;

% Find positions
i_SSAy   = ismember(varendo,SSAy);
i_SSAe   = ismember(varendo,SSAe);
if ~isempty(impulse)
iImpulse = ismember(varendo,impulse);
end

%-------------------------------------------------------------------------%
% PART TO CUSTOMISE

choose = 'default_irf';

switch choose

    case 'default_irf'

    % Conditions on observables
    yCondition = nan(h,N);                                 %nan means no restriction
    yCondition(:,i_SSAy) = zeros(h,sum(i_SSAy));           %set variables SSAy to a specific path (zero is steady state)
    Omega_f = zeros(sum(sum(isfinite(yCondition))));       %"hard" restrictions
    
    % Conditions on innovations/shocks
    eCondition = zeros(h,N);                               %zero means a shock is set to its unconditional distribution
    eCondition(:,i_SSAe) = nan;                            %SSAe shocks can deviate from their unconditional distribution (you can specify for which horizon)
    eCondition(1,iImpulse) = 1;                            %shockVar shocks set to one on impact
    %Omega_g = 1*diag(ones(sum(sum(isfinite(eCondition))),1)); %Non-driving shocks restricted to their unconditional variance
    Omega_g = zeros(sum(sum(isfinite(eCondition))));       %Non-driving shocks restricted to zero ("hard" restrictions)

end
%-------------------------------------------------------------------------%

% Pack output
out.yCondition = yCondition;
out.Omega_f    = Omega_f;
out.eCondition = eCondition;
out.Omega_g    = Omega_g;
out.chooseSS   = choose;