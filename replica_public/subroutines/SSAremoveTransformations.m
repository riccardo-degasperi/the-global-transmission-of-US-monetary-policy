function out = SSAremoveTransformations(in,data)
% SSAremoveTransformations removes the transformations applied to the time
% series in the estimation phase.

% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpack
iLogs  = logical(data.iLogs);
iFD    = logical(data.iFD);

% Initialise
out = in;

% Transform
out(:,iLogs,:) = exp(in(:,iLogs,:)./100);


% EXTEND TO FIRST DIFFERENCES AND EXOGENOUS