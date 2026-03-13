function [out,dateout] = SSAgetSample(in,datein0,frequency,lT,dims)
% SSAgetSample selects the sample to show in SSAplotScenarios and returns
% a compatible array of dates in datenum format.


% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpack
p = dims.p;
h = dims.h+1;

% Remove initial conditions
datein = datein0(p+1:end);

% Find cutoff
pos = strcmp(datein,lT);
if sum(pos) == 0
    error(['ERROR: there is no date that matches the selected cutoff, ',lT,'.'])
elseif sum(pos) > 1
    error(['ERROR: there is more than one date that matches the selected cutoff, ',lT,'.'])
end

% Get cutoff observation
first0 = datein(pos);
first = first0{:};

% Get last observation
last0 = datein0(end);
last = last0{:};

if strcmp(frequency,'monthly') || strcmp(frequency,'quarterly')

    % Get years and months/quarters
    yy0 = str2double(first(1:4));
    mm0 = str2double(first(6:end));
    yy1 = str2double(last(1:4));
    mm1 = str2double(last(6:end));
    
    % Determine date format
    id0 = first(5);
    id1 = last(5);
    if strcmp(id0,id1) == 0
        error('ERROR: the selected cutoff is of a different format with respect to the imputed dates.')
    end

    if strcmp(frequency,'monthly')

        dat0 = datetime(yy0,mm0,1);  %cutoff date
        dat1 = datetime(yy1,mm1,1);  %date of last data point
        dat2 = dat1+calmonths(h);    %final date

        % Get dates
        dateout = (dat0:calmonths(1):dat2)';

    elseif strcmp(frequency,'quarterly')

        dat0 = datetime(yy0,mm0+2*(mm0-1),1);  %cutoff date
        dat1 = datetime(yy1,mm1+2*(mm1-1),1);  %date of last data point
        dat2 = dat1+calmonths(3*h);  %final date

        % Get dates
        dateout = (dat0:calmonths(3):dat2)';

    end

elseif strcmp(frequency,'daily') || strcmp(frequency,'weekly')

    % Get years, months, and days
    yy0 = str2double(first(7:10));
    mm0 = str2double(first(4:5));
    dd0 = str2double(first(1:2));
    yy1 = str2double(last(7:10));
    mm1 = str2double(last(4:5));
    dd1 = str2double(last(1:2));

    if strcmp(frequency,'daily')

        dat0 = datetime(yy0,mm0,dd0);  %cutoff date
        dat1 = datetime(yy1,mm1,dd1);  %date of last data point
        dat2 = dat1+caldays(h);        %final date
    
        % Get dates
        dateout = (dat0:caldays(1):dat2)';

    elseif strcmp(frequency,'weekly')

        dat0 = datetime(yy0,mm0,dd0);  %cutoff date
        dat1 = datetime(yy1,mm1,dd1);  %date of last data point
        dat2 = dat1+calweeks(h);       %final date

        % Get dates
        dateout = (dat0:calweeks(1):dat2)';

    end


elseif strcmp(frequency,'yearly')

    % Get years
    yy0 = str2double(first);
    yy1 = str2double(last);

    dat0 = datetime(yy0,1,1);  %cutoff date
    dat1 = datetime(yy1,1,1);  %date of last data point
    dat2 = dat1+calyears(h);   %final date

    % Get dates
    dateout = (dat0:calyears(1):dat2)';
    
else
    error('ERROR: dates have to be defined as follows: 2023m1, 2023q2, 2023, 25/12/2023 both in the Excel source file and in the options.')
end

% Match data
cut = size(in,1) - numel(dateout);
out = in(cut+1:end,:,:);

