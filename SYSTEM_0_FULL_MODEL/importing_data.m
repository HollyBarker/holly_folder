diffvec=[0.1,0.2,0.4,0.6,0.8,1.0];
stepvec=[0,5,26,58,110,190,300,364,450,700];

baseName='/home/mbcx3hb2/Documents/First_Year_Report/Graphs/Data/CONCENTRATION_TIMESTEPS/diff_0._';
stepName='step'
fileType='.dat';

Solution_Concentration=zeros(5000,6,10);

for stepnum=1:1:10
    step=stepvec(stepnum); 

    %import diff=1
    filename =  [baseName,'10',stepName,num2str(step),fileType];
    delimiter = ' ';
    formatSpec = '%s%s%[^\n\r]';
    fileID = fopen(filename,'r');
    
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
    fclose(fileID);
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    
    for col=[1,2]
        % Converts strings in the input cell array to numbers. Replaced non-numeric
        % strings with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1);
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric strings to numbers.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end
    
    I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
    raw(I,:) = [];
    
    x= cell2mat(raw(:, 1));
    Solution_Concentration(:,6,stepnum) = cell2mat(raw(:, 2));
    
    
    
    
    
    
    
    
    
    
    
    %import the other diff values
    for dnum=1:1:5
        diff=diffvec(dnum);
        diffFile=diff*10;
        
        filename =  [baseName,num2str(diffFile),stepName,num2str(step),fileType];
        delimiter = ' ';
        formatSpec = '%s%s%[^\n\r]';
        fileID = fopen(filename,'r');
        
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
        fclose(fileID);
        raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
        for col=1:length(dataArray)-1
            raw(1:length(dataArray{col}),col) = dataArray{col};
        end
        numericData = NaN(size(dataArray{1},1),size(dataArray,2));
        
        for col=[1,2]
            % Converts strings in the input cell array to numbers. Replaced non-numeric
            % strings with NaN.
            rawData = dataArray{col};
            for row=1:size(rawData, 1);
                % Create a regular expression to detect and remove non-numeric prefixes and
                % suffixes.
                regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                try
                    result = regexp(rawData{row}, regexstr, 'names');
                    numbers = result.numbers;
                    
                    % Detected commas in non-thousand locations.
                    invalidThousandsSeparator = false;
                    if any(numbers==',');
                        thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                        if isempty(regexp(numbers, thousandsRegExp, 'once'));
                            numbers = NaN;
                            invalidThousandsSeparator = true;
                        end
                    end
                    % Convert numeric strings to numbers.
                    if ~invalidThousandsSeparator;
                        numbers = textscan(strrep(numbers, ',', ''), '%f');
                        numericData(row, col) = numbers{1};
                        raw{row, col} = numbers{1};
                    end
                catch me
                end
            end
        end
        
        I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
        raw(I,:) = [];
        
        x= cell2mat(raw(:, 1));
        Solution_Concentration(:,dnum,stepnum) = cell2mat(raw(:, 2));
        
        
    end
end


clearvars baseName fileType dnum stepnum diffvec stepvec diff step filename delimiter formatSpec fileID;
clearvars diffFile stepName ans row raw rawData regexstr result dataArray I invalidThousandsSeparator me numbers numericData;
clearvars col;