function [outes,outiv,lTes,lTiv,uTes,uTiv] = ivMatch(p,dates,ivDates,es,iv)
% ivMatch matches the length of the external instrument and the reduced
% form VAR residuals. E.g. if the data covers the period 1975m1-2014m12 and
% the IV covers the period 1990m1-2016m1, then ivMatch trims the
% residuals (es) prior to 1990m1, while dataLoader already trimmed the IV
% after 2014m12.

% INPUT:
% - p       : number of lags in the system
% - dates   : dates (string) covered by the data, including initial p obs
% - ivDates : dates (string) covered by the external instrument
% - es      : matrix of residuals of the reduced-form VAR (T-p x N x sim)
% - iv      : vector of the external instrument (can also be more than one)

% OUTPUT:
% - outes   : 
% - outiv   :
% - lTes    :
% - lTiv    : 
% - uTes    :
% - uTiv    :

% r.degasperi@warwick.ac.uk

%-------------------------------------------------------------------------%
% Find initial date if IV starts after data
lTes = find(strcmp(ivDates(1),dates(p+1:end)),1);

% Find initial date if IV starts before data 
lTiv = find(strcmp(dates(p+1),ivDates(:)),1);

% If both lTes and lTiv are empty, there is an error
if isempty(lTes) && isempty(lTiv)
   error('ERROR: Sample matching failed. Probably IV and residuals do not overlap in time: the series you are using are too short.') 
end


%-------------------------------------------------------------------------%
% Find final date if IV ends before data
uTes = find(strcmp(ivDates(end),dates(p+1:end)),1);

% Find final date if IV ends after data
uTiv = find(strcmp(dates(end),ivDates(:)),1);

% If uTes is empty, there is an error
if isempty(uTes) && isempty(uTiv)
   error('ERROR: Sample matching failed. Probably IV and residuals do not overlap in time: the series you are using are too short.') 
end


%-------------------------------------------------------------------------%
% Match samples

if      isempty(lTiv) && isempty(uTiv)      %IV starts after and ends before
        outes = es(lTes:uTes,:,:);
        outiv = iv;
       
elseif isempty(lTes) && isempty(uTes)      %IV starts before and ends after
        outes = es;
        outiv = iv(lTiv:uTiv,:);
       
elseif isempty(lTiv) && isempty(uTes)      %IV starts after and ends after
        outes = es(lTes:end,:,:);
        outiv = iv(1:uTiv,:);
       
elseif isempty(lTes) && isempty(uTiv)      %IV starts before and ends before
        outes = es(1:uTes,:,:);
        outiv = iv(lTiv:end,:);
       
elseif ~isempty(lTiv) && ~isempty(lTes) && ~isempty(uTiv) && ~isempty(uTes)
       
        if (lTiv == lTes) && (uTes == uTiv)
            outes = es;
            outiv = iv;
        else
            error('ERROR: Sample matching failed. Check ivMatch.m')
        end
   
elseif ~isempty(lTiv) && ~isempty(lTes) && (~isempty(uTiv) || ~isempty(uTes))
    
        if lTiv == lTes == 1                   %IV and data have the same start
           if isempty(uTes)
           outes = es;
           outiv = iv(1:uTiv,:);
           elseif isempty(uTiv)
           outes = es(1:uTes,:,:);
           outiv = iv;
           end
        else
            error('ERROR: Sample matching failed. Check ivMatch.m')
        end
    
elseif (~isempty(lTiv) || ~isempty(lTes)) && ~isempty(uTiv) && ~isempty(uTes)
    
        if strcmp(ivDates(uTiv),dates(p+uTes)) %IV and data have the same end
           if isempty(lTes)
           outes = es;
           outiv = iv(lTiv:end,:);
           elseif isempty(lTiv)
           outes = es(lTes:end,:,:);
           outiv = iv;
           end
        else
            error('ERROR: Sample matching failed. Check ivMatch.m')
        end

else
    error('ERROR: Sample matching failed. Check ivMatch.m')
end



% OLD: matches only when IV is shorter or equal to data on both ends.
% % Match IV and data samples
% lT = find(strcmp(ivDates(1),dates(p+1:end)));        %find initial date
% if isempty(lT)
%    error('**iv series matching failed**') 
% end
% 
% uT = find(strcmp(ivDates(end),dates(p+1:end)));      %find final date
% if isempty(uT)
%    error('**iv series matching failed**') 
% end
% 
% % Match samples
% res = es(lT:uT,:,:);





