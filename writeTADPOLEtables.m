function [ submission_table ] = writeTADPOLEtables(rid,ADAS13_forecast,Ventricles_ICV_forecast,CLIN_STAT_forecast,ci,fn,sheet)

% display(sprintf('Constructing the output spreadsheet %s...', outputFile))
N_D2 = length(rid);
startDate = datenum('01-Jan-2018');
nForecasts = 60;
submission_table =  cell2table(cell(N_D2*nForecasts,12), ...
  'VariableNames', {'RID', 'ForecastMonth', 'ForecastDate', ...
  'CNRelativeProbability', 'MCIRelativeProbability', 'ADRelativeProbability', ...
  'ADAS13', 'ADAS1350_CILower', 'ADAS1350_CIUpper', ...
  'Ventricles_ICV', 'Ventricles_ICV50_CILower', 'Ventricles_ICV50_CIUpper' });
%* Repeated matrices - compare with submission template
submission_table.RID = reshape(repmat(rid, [1, nForecasts])', N_D2*nForecasts, 1);
submission_table.ForecastMonth = repmat((1:nForecasts)', [N_D2, 1]);
%* First subject's submission dates
for m=1:nForecasts
  submission_table.ForecastDate{m} = datestr(addtodate(startDate, m-1, 'month'), 'yyyy-mm');
end
%* Repeated matrices for submission dates - compare with submission template
submission_table.ForecastDate = repmat(submission_table.ForecastDate(1:nForecasts), [N_D2, 1]);

%* Pre-fill forecast data, encoding missing data as NaN
nanColumn = nan(size(submission_table.CNRelativeProbability));
submission_table.CNRelativeProbability = nanColumn;
submission_table.MCIRelativeProbability = nanColumn;
submission_table.ADRelativeProbability = nanColumn;
submission_table.ADAS13 = nanColumn;
submission_table.ADAS1350_CILower = nanColumn;
submission_table.ADAS1350_CIUpper = nanColumn;
submission_table.Ventricles_ICV = nanColumn;
submission_table.Ventricles_ICV50_CILower = nanColumn;
submission_table.Ventricles_ICV50_CIUpper = nanColumn;

col = 4;
t = CLIN_STAT_forecast(:,:,1)';
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);

 %*  b) MCI probabilities
col = 5;
t = CLIN_STAT_forecast(:,:,2)';
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);
  %*  c) AD probabilities
col = 6;
t = CLIN_STAT_forecast(:,:,3)';
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);
%* 2. ADAS13 score
col = 7;
t = ADAS13_forecast(:,:)';
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);
  %*  a) Lower and upper bounds (50% confidence intervals)
col = 8;
t = ADAS13_forecast(:,:)' - ci(1);
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);
col = 9;
t = ADAS13_forecast(:,:)' + ci(1);
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);
%* 3. Ventricles volume (normalised by intracranial volume)
col = 10;
t = Ventricles_ICV_forecast(:,:)';
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);
  %*  a) Lower and upper bounds (50% confidence intervals)
col = 11;
t = Ventricles_ICV_forecast(:,:)' - ci(2);
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);
col = 12;
t = Ventricles_ICV_forecast(:,:)' + ci(2);
col = submission_table.Properties.VariableNames(col);
submission_table{:,col} = t(:);

%* Convert all numbers to strings
hdr = submission_table.Properties.VariableNames;
for k=1:length(hdr)
  if ~iscell(submission_table.(hdr{k}))
    %submission_table{1:10,hdr{k}} = varfun(@num2str,submission_table{1:10,hdr{k}},'OutPutFormat','cell');
    submission_table.(hdr{k}) = strrep(cellstr(num2str(submission_table{:,hdr{k}},'%0.5f')),' ','');
  end
end

%* Use column names that match the submission template
columnNames = {'RID', 'Forecast Month', 'Forecast Date',...
'CN relative probability', 'MCI relative probability', 'AD relative probability',	...
'ADAS13',	'ADAS13 50% CI lower', 'ADAS13 50% CI upper', 'Ventricles_ICV', ...
'Ventricles_ICV 50% CI lower',	'Ventricles_ICV 50% CI upper'};
%* Convert table to cell array to write to file, line by line
%  This is necessary because of spaces in the column names: writetable()
%  doesn't handle this.
tableCell = table2cell(submission_table);
tableCell = [columnNames;tableCell];
xlswrite(fn,tableCell,sheet);

end

