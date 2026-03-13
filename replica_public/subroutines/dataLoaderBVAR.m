function [out,dataOpt] = dataLoaderBVAR(modelOpt,dataOpt,miscOpt)
% dataLoaderBVAR loads multi-country data from Excel source and applies
% transformations.
%
% Last modified: 26/03/2024
% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%
disp('    Load data and apply transformations')

% Specify number of rows containing text in unit sheets
sp = 6;        

% Unpacking
identification = modelOpt.identification;
lT             = dataOpt.lT;
uT             = dataOpt.uT;
dataset        = dataOpt.dataset;
shockVar       = dataOpt.shockVar;
varendo        = dataOpt.varendo;
varendo_global = dataOpt.varendo_global;
varexo         = dataOpt.varexo;
units          = dataOpt.units;
frequency      = dataOpt.frequency;
dateCheck      = miscOpt.dateCheck;
interpolate    = miscOpt.interpolate;
interpolCase   = miscOpt.interpolCase;
ws             = miscOpt.ws;

%-------------------------------------------------------------------------%
% Identify sheets
[~,raw.sheets] = xlsfinfo(dataset);          %names of excel sheets

sGlo   = find(strcmp('Global',raw.sheets));  %position of 'Global' sheet
if isempty(sGlo)
    warning('Sheet "Global" is missing: no exogenous variables')
end
sIV    = find(strcmp('IV',raw.sheets));      %position of 'IV' sheet
if isempty(sIV)
    warning('Sheet "IV" is missing: no instrumental variables')
end
sTrans = find(strcmp('Trans',raw.sheets));   %position of 'Trans' sheet
if isempty(sTrans)
    error('Sheet "Trans" is missing: check Excel data file')
end
sUnit  = NaN(1,size(units,2));               %container for countries' sheets position
for i = 1:size(units,2)
    tmp = strcmp(units{i},raw.sheets);
    if ~ismember(1,tmp)
        error(['ERROR: unit ' units{i} ' cannot be found'])
    end
    sUnit(i) = find(tmp);                    %position of countries' sheets
end

sAll = [sGlo,sIV,sTrans,sUnit];              %combine all sheets positions


%-------------------------------------------------------------------------%
% Initialise containers
sample = NaN(size(units,2)+2,2);             %plus 2 for global endo and exo
endo0  = cell(1,size(units,2));

%-------------------------------------------------------------------------%
% Create matrices of endogenous and exogenous variables
for i = sAll  %loop over sheets
    
    % Retrieve data from Excel
    [raw.data.(raw.sheets{i}),raw.labels.(raw.sheets{i}),raw.raw.(raw.sheets{i})]=xlsread(dataset,i); 
    
    % NOTE: we store also the raw information as a cell (raw.raw) to solve
    % a problem related to xlsread trimming leading and trailing NaNs
    % automatically. In the raw.raw case, this does not happen.
    

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
    % Extract endogenous variables for each unit
    if ismember(i,sUnit)   %excluding 'Global', 'IV' and 'Trans' sheets
        
        % Find position of endogenous vars for unit i
        sEndo = NaN(1,size(varendo,2));   %container
        for ii = 1:size(varendo,2)
            tmp = strcmp(varendo{ii},raw.labels.(raw.sheets{i})(1,2:end));
            if ~ismember(1,tmp)
                error(['ERROR: endogenous variable ' varendo{ii} ' cannot be found'])
            end
            sEndo(ii) = find(tmp);        %position
        end
        
        % Find name of unit i
        ii = find(strcmp(raw.sheets(i),units));
        
        % Find longest possible sample
        tmp   = raw.raw.(raw.sheets{i})(7:end,sEndo+1); %data
        iNaN  = cellfun(@isnan,tmp);                    %find NaNs
        iNaN2 = (sum(iNaN,2)==0);                       %find rows without NaNs
        if sum(iNaN2)==0
           error(['ERROR: empty variable in unit ',units{ii}]) 
        end
        
        % Lower cutoff
        if strcmp(lT,'full')
            lTp = find(iNaN2,1,'first');  
        else
            tmp = strcmp(lT,raw.labels.(raw.sheets{i})(7:end,1));
            lTp = find(tmp);
            if isempty(lTp)
                error('ERROR: date cutoffs outside data range')
            elseif iNaN2(lTp) == 0
                warning(['date lower cutoffs outside data range for unit ',units{ii},'. Selecting longest available sample...'])
                lTp = find(iNaN2,1,'first');
                if dateCheck
                    lTdisp = raw.labels.(raw.sheets{i})(lTp+sp,1);
                    prompt = ['Starting date: ',char(lTdisp),'. Is that ok? (1:yes/0:no)  '];
                    in = input(prompt);
                    if in == 1
                        %do nothing
                    elseif in == 0
                        error('ABORTED BY USER. Pick a different starting date.')
                    else
                        error('ERROR: answer 1 for yes or 0 for no.')
                    end
                end
            end
        end
        
        % Upper cutoff
        if strcmp(uT,'full')
            uTp = find(iNaN2,1,'last');
        else
            tmp = strcmp(uT,raw.labels.(raw.sheets{i})(7:end,1));
            uTp = find(tmp);
            if isempty(uTp)
                error('ERROR: date cutoffs outside data range')
            elseif iNaN2(uTp) == 0
                warning(['date upper cutoffs outside data range for unit ',units{ii},'. Selecting longest available sample...'])
                uTp = find(iNaN2,1,'last');
                if dateCheck
                    uTdisp = raw.labels.(raw.sheets{i})(uTp+sp,1);
                    prompt = ['Ending date: ',char(uTdisp),'. Is that ok? (1:yes/0:no)  '];
                    in = input(prompt);
                    if in == 1
                        %do nothing
                    elseif in == 0
                        error('ABORTED BY USER. Pick a different ending date.')
                    else
                        error('ERROR: answer 1 for yes or 0 for no.')
                    end
                end
            end
        end

        % Store matrix of endogenous vars for unit i
        endo0{ii} = raw.raw.(raw.sheets{i})(:,sEndo+1);
        
        % Store date range
        sample(ii,:) = [lTp,uTp];
        
    end
    

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
    % Extract global (exogenous and endogenous) variables
    if i == sGlo
        
        %  Find position of exogenous vars
        sExo = NaN(1,size(varexo,2));   %container
        for ii = 1:size(varexo,2)
            tmp = strcmp(varexo{ii},raw.labels.(raw.sheets{i})(1,2:end));
            if ~ismember(1,tmp)
                error(['ERROR: exogenous variable ' varexo{ii} ' cannot be found'])
            end
            sExo(ii) = find(tmp);       %position
        end
        
        % Retrieve matrix of exogenous vars
        if isempty(sExo)
            exo = [];
        else
            
            % Find longest possible sample
            tmp   = raw.raw.(raw.sheets{i})(7:end,sExo+1);  %data
            iNaN  = cellfun(@isnan,tmp);                    %find NaNs
            iNaN2 = (sum(iNaN,2)==0);                       %find rows without NaNs
            if sum(iNaN2)==0
                error('ERROR: empty variable in sheet Global') 
            end
            
            % Lower cutoff
            if strcmp(lT,'full')
                lTp = find(iNaN2,1,'first');
            else
                tmp = strcmp(lT,raw.labels.(raw.sheets{i})(7:end,1));
                lTp = find(tmp);
                if isempty(lTp)
                    error('ERROR: date cutoffs outside data range')
                elseif iNaN2(lTp) == 0
                    warning('date lower cutoffs outside data range for sheet Global. Selecting longest available sample...')
                    lTp = find(iNaN2,1,'first');
                    if dateCheck
                        lTdisp = raw.labels.(raw.sheets{i})(lTp+sp,1);
                        prompt = ['Starting date: ',char(lTdisp),'. Is that ok? (1:yes/0:no)  '];
                        in = input(prompt);
                        if in == 1
                            %do nothing
                        elseif in == 0
                            error('ABORTED BY USER. Pick a different starting date.')
                        else
                            error('ERROR: answer 1 for yes or 0 for no.')
                        end
                    end
                end
            end

            % Upper cutoff
            if strcmp(uT,'full')
                uTp = find(iNaN2,1,'last');
            else
                tmp = strcmp(uT,raw.labels.(raw.sheets{i})(7:end,1));
                uTp = find(tmp);
                if isempty(uTp)
                    error('ERROR: date cutoffs outside data range')
                elseif iNaN2(uTp) == 0
                    warning('date upper cutoffs outside data range for sheet Global. Selecting longest available sample...')
                    uTp = find(iNaN2,1,'last');
                    if dateCheck
                        uTdisp = raw.labels.(raw.sheets{i})(uTp+sp,1);
                        prompt = ['Ending date: ',char(uTdisp),'. Is that ok? (1:yes/0:no)  '];
                        in = input(prompt);
                        if in == 1
                            %do nothing
                        elseif in == 0
                            error('ABORTED BY USER. Pick a different ending date.')
                        else
                            error('ERROR: answer 1 for yes or 0 for no.')
                        end
                    end
                end
            end
            
            % Store matrix of exogenous variables
            exo0 = raw.raw.(raw.sheets{i})(:,sExo+1);

            % Store date range
            sample(end-1,:) = [lTp,uTp];

        end
        
        % Find position of endogenous vars in global
        sEndo2 = NaN(1,size(varendo_global,2));   %container
        for ii = 1:size(varendo_global,2)
            tmp = strcmp(varendo_global{ii},raw.labels.(raw.sheets{i})(1,2:end));
            if ~ismember(1,tmp)
                error(['ERROR: global endogenous variable ' varendo_global{ii} ' cannot be found'])
            end
            sEndo2(ii) = find(tmp);       %position
        end
        
        % Retrieve matrix of endogenous global vars
        if isempty(sEndo2)
            endo_global = [];
        else
            
            
            %===============================%
            
            % Find longest possible sample
            tmp   = raw.raw.(raw.sheets{i})(7:end,sEndo2+1);  %data
            iNaN  = cellfun(@isnan,tmp);                      %find NaNs
            iNaN2 = (sum(iNaN,2)==0);                         %find rows without NaNs
            if sum(iNaN2)==0
                error(['ERROR: empty variable in unit ',units{ii}]) 
            end
            
            % Lower cutoff
            if strcmp(lT,'full')
                lTp = find(iNaN2,1,'first');
            else
                tmp = strcmp(lT,raw.labels.(raw.sheets{i})(7:end,1));
                lTp = find(tmp);
                if isempty(lTp)
                    error('ERROR: date cutoffs outside data range')
                elseif iNaN2(lTp) == 0
                    warning('date lower cutoffs outside data range for sheet Global. Selecting longest available sample...')
                    lTp = find(iNaN2,1,'first');
                    if dateCheck
                        lTdisp = raw.labels.(raw.sheets{i})(lTp+sp,1);
                        prompt = ['Starting date: ',char(lTdisp),'. Is that ok? (1:yes/0:no)  '];
                        in = input(prompt);
                        if in == 1
                            %do nothing
                        elseif in == 0
                            error('ABORTED BY USER. Pick a different starting date.')
                        else
                            error('ERROR: answer 1 for yes or 0 for no.')
                        end
                    end
                end
            end

            % Upper cutoff
            if strcmp(uT,'full')
                uTp = find(iNaN2,1,'last');
            else
                tmp = strcmp(uT,raw.labels.(raw.sheets{i})(7:end,1));
                uTp = find(tmp);
                if isempty(uTp)
                    error('ERROR: date cutoffs outside data range')
                elseif iNaN2(uTp) == 0
                    warning('date upper cutoffs outside data range for sheet Global. Selecting longest available sample...')
                    uTp = find(iNaN2,1,'last');
                    if dateCheck
                        uTdisp = raw.labels.(raw.sheets{i})(uTp+sp,1);
                        prompt = ['Ending date: ',char(uTdisp),'. Is that ok? (1:yes/0:no)  '];
                        in = input(prompt);
                        if in == 1
                            %do nothing
                        elseif in == 0
                            error('ABORTED BY USER. Pick a different ending date.')
                        else
                            error('ERROR: answer 1 for yes or 0 for no.')
                        end
                    end
                end
            end
            
            % Store matrix of global endogenous variables
            endo_global0 = raw.raw.(raw.sheets{i})(:,sEndo2+1);
            
            % Store date range
            sample(end,:) = [lTp,uTp];

        end
    end
    

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
    % Extract IV variables
    
    if i == sIV

        if strcmp(identification,'iv')

        % Unpack
        variv = dataOpt.variv;
        
        %  Find position of IV vars
        iIV = NaN(1,size(variv,2));   %container
        for ii = 1:size(variv,2)
            tmp = strcmp(variv{ii},raw.labels.(raw.sheets{i})(1,2:end));
            if ~ismember(1,tmp)
                error(['ERROR: IV variable ' variv{ii} ' cannot be found'])
            end
            tmp2 = find(tmp);       %position
            if numel(tmp2) == 1
                iIV(ii) = tmp2;
            else
                error('ERROR: there are at least 2 IVs with the same name in sheet "IV".')
            end
        end
        
        % Retrieve matrix of IV vars
        iv = raw.data.(raw.sheets{i})(:,iIV);
        
        % Eliminate NaNs
        iNaNiv = iv ~= 123456789;              %find NaNs
        iNaNiv = sum(iNaNiv,2) == size(iv,2);  %select common sample
        iv = iv(iNaNiv,:);                     %drop NaNs
        
        % Retrieve dates
        iNaNiv = logical([0;iNaNiv]);                   %remove label
        ivDates = raw.labels.(raw.sheets{i})(iNaNiv,1);

        end
    end
    

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
    % Extract transformation indices and long names for endogenous vars
    if i == sTrans
        
        % For endogenous variables
        iTrEn = NaN(1,size(varendo,2));     %container
        for ii = 1:size(varendo,2)
            tmp = strcmp(varendo{ii},raw.labels.(raw.sheets{i})(2:end,1));
            if ~ismember(1,tmp)
                error(['ERROR: variable ' varendo{ii} ' is not in sheet "Trans"'])
            end
            iTrEn(ii) = find(tmp,1,'first');  %position
        end
        
        % Retrieve transformation arrays
        iLogsEn = raw.data.(raw.sheets{i})(iTrEn,1);
        posLogsEn = find(raw.data.(raw.sheets{i})(iTrEn,1))';
        
        % RW index (used only in bayesian analysis)
        iRW_endo = raw.data.(raw.sheets{i})(iTrEn,2);    %appicable only to endogenous vars
               
        % ylabel (for charts y-axis labels)
        ylab_endo = raw.labels.(raw.sheets{i})(iTrEn+1,6);
        
        % Long names for endogenous variables
        varendolong = raw.labels.(raw.sheets{i})(iTrEn+1,2)';
        
        % For exogenous variables
        iTrEx = NaN(1,size(varexo,2));   %container
        for ii = 1:size(varexo,2)
            tmp = strcmp(varexo{ii},raw.labels.(raw.sheets{i})(2:end,1));
            if ~ismember(1,tmp)
                error(['ERROR: variable ' varexo{ii} ' is not in sheet "Trans"'])
            end
            iTrEx(ii) = find(tmp);       %position
        end
        
        % Retrieve transformation arrays
        iLogsEx = raw.data.(raw.sheets{i})(iTrEx,1);
        posLogsEx = find(raw.data.(raw.sheets{i})(iTrEx,1))';
        
        % ylabel (for charts y-axis labels)
        ylab_exo = raw.labels.(raw.sheets{i})(iTrEx+1,6);
        
        % Long names for exogenous variables
        varexolong = raw.labels.(raw.sheets{i})(iTrEx+1,2)';
        
        % For endogenous global variables
        iTrEn2 = NaN(1,size(varendo_global,2));   %container
        for ii = 1:size(varendo_global,2)
            tmp = strcmp(varendo_global{ii},raw.labels.(raw.sheets{i})(2:end,1));
            if ~ismember(1,tmp)
                error(['ERROR: variable ' varendo_global{ii} ' is not in sheet "Trans"'])
            end
            iTrEn2(ii) = find(tmp);       %position
        end
        
        % Retrieve transformation arrays
        iLogsEn_global = raw.data.(raw.sheets{i})(iTrEn2,1);
        posLogsEn_global = find(raw.data.(raw.sheets{i})(iTrEn2,1))';

        % RW index (used only in bayesian analysis)
        iRW_global = raw.data.(raw.sheets{i})(iTrEn2,2);    %appicable only to endogenous vars
        
        % ylabel (for charts y-axis labels)
        ylab_global = raw.labels.(raw.sheets{i})(iTrEn2+1,6);
        
        % Long names for endogenous global variables
        varendolong_global = raw.labels.(raw.sheets{i})(iTrEn2+1,2)';
        
    end
    
end

%-------------------------------------------------------------------------%
% Assemble data matrices
lTf = max(sample(:,1));  %minimum common starting date
uTf = min(sample(:,2));  %maximum common ending date

for j = 1:size(units,2)
    endo(:,:,j) = cell2mat(endo0{j}(lTf+sp:uTf+sp,:));
end

if ~isempty(sGlo)
    if ~isempty(sExo)
        exo = cell2mat(exo0(lTf+sp:uTf+sp,:));
    end

    if ~isempty(sEndo2)
        endo_global = cell2mat(endo_global0(lTf+sp:uTf+sp,:));
    end
end

%-------------------------------------------------------------------------%
% Retrieve dates
if strcmp(frequency,'yearly')
    dates = raw.data.(raw.sheets{i})(lTf:uTf,1);  %data dates
    dates = string(dates);
else
    dates = raw.labels.(raw.sheets{i})(lTf+sp:uTf+sp,1);  %data dates
end
if dateCheck
    prompt = ['Common sample from ',char(dates(1)),' to ',char(dates(end)),'. Is that ok? (1:yes/0:no):  '];
    in = input(prompt);
    if in == 1
        %do nothing
    elseif in == 0
        error('ABORTED BY USER.')
    else
        error('ERROR: answer 1 for yes or 0 for no.')
    end
                        
end

%-------------------------------------------------------------------------%
% Check for NaNs & Interpolation
for i = 1:size(units,2)

    [~,tmp] = find(isnan(endo(:,:,i)),1);
    if ~isempty(tmp) && ~interpolate
        error(['ERROR: there are NaNs in the series ', varendo{tmp},' for unit ', units{i}])
    elseif ~isempty(tmp) && interpolate
        warning(['there are NaNs in the series ', varendo{tmp},' for unit ', units{i},'. Interpolating...'])
        endo(:,:,i) = interpolation(endo(:,:,i),interpolCase,ws);
    end
end

if ~isempty(sGlo)
    [~,tmp] = find(isnan(exo),1);
    if ~isempty(tmp) && ~interpolate
       error(['ERROR: there are NaNs in the series ', varexo{tmp}])
    elseif ~isempty(tmp) && interpolate
       warning(['there are NaNs in the series ', varexo{tmp},'. Interpolating...'])
       exo = interpolation(exo,interpolCase,ws); 
    end

    [~,tmp] = find(isnan(endo_global),1);
    if ~isempty(tmp) && ~interpolate
       error(['ERROR: there are NaNs in the series ', varendo_global{tmp}])
    elseif ~isempty(tmp) && interpolate
       warning(['there are NaNs in the series ', varendo_global{tmp},'. Interpolating...'])
       endo_global = interpolation(endo_global,interpolCase,ws); 
    end

end

%-------------------------------------------------------------------------%
% Check that shockVar is in varendo or varendo_global:
for i = 1:numel(shockVar)
    if isempty(find(strcmp(shockVar{i},varendo),1)) && isempty(find(strcmp(shockVar{i},varendo_global),1))
        error(['ERROR: variable ',shockVar{i},' not included in the list of endogenous variables.'])
    end
end


%-------------------------------------------------------------------------%
% Log transformations
endo(:,posLogsEn,:) = 100*reallog(endo(:,posLogsEn,:));  %transform endogenous

if isempty(sGlo)
    exo = [];
else
    exo(:,posLogsEx) = 100*reallog(exo(:,posLogsEx));    %transform exogenous
end

if isempty(sGlo)
    endo_global = [];
else
    endo_global(:,posLogsEn_global) = 100*reallog(endo_global(:,posLogsEn_global));  %transform exogenous
end

%-------------------------------------------------------------------------%
% Pack output
out.dates                 = dates;
out.endo                  = endo;
out.endo_global           = endo_global;
out.exo                   = exo;
out.iRW_endo              = iRW_endo;
out.iRW_global            = iRW_global;
out.iLogsEn               = iLogsEn;
out.iLogsEn_global        = iLogsEn_global;
out.iLogsEx               = iLogsEx;
if strcmp(modelOpt.identification,'iv')
    out.iv                    = iv;
    out.ivDates               = ivDates;
end
dataOpt.varendolong           = varendolong;
dataOpt.varendolong_global    = varendolong_global;
dataOpt.varexolong            = varexolong;
dataOpt.ylab_endo             = ylab_endo;
dataOpt.ylab_exo              = ylab_exo;
dataOpt.ylab_global           = ylab_global;


