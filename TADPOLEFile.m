classdef TADPOLEFile
%Wrapper class for accessing TADPOLE datasets
    
    properties(Access=public) 
        MODALITY_PATH= strcat('c:\Users\justoh\Repos\Tadpole2',filesep,'modalities',filesep);
        FILE_NAME = 'TADPOLE_D1_D2.csv';      
        PTID_COLUMN = 'PTID';
        VISCODE_COLUMN = 'VISCODE';
        MONTH_FACTOR_COLUMN = 'MONTH';
        MONTH_COLUMN = 'M';
        D1_COLUMN = 'D1';
        D2_COLUMN = 'D2';
        modalities = containers.Map
        data_file = [];
        d_values;
        EXAMDATE_COLUMN = 9; 
    end
  
    methods
        

        function obj = TADPOLEFile( d_values )
%       Constructor for the wrapper.
%       Attributes:
%           d_values: Specifies the combination of d1 and d2 rows to be
%           include. 1x2 cell array of chars. Where the first cell indicates d1
%           value and second cell d2. OPTIONAL (Default{'1', '0'}).

            if nargin < 1
                d_values = {'1','0'};         
            end

            obj.d_values = d_values;
            
            obj = obj.setDataset( d_values )
            
            %Replace bl with more flexible m00
            obj.data_file{  ismember(table2array(obj.data_file(:,3)),'bl') , 3 } = {'m00'};            
           
            %Load the columns for each modality from separate txt files.
            %Key values of the modality map define names for each modality.
            obj.modalities('cognitive') = obj.read_columns('adni_cognitive_1.txt');
            obj.modalities('cognitive_bl') = obj.read_columns('adni_cognitive_1_bl.txt');
            obj.modalities('moca_ecog') = obj.read_columns('adni_cognitive_moca_egoc.txt');
            obj.modalities('moca_ecog_bl') = obj.read_columns('adni_cognitive_moca_egoc_bl.txt');   
            
            obj.modalities('ucsf') = obj.read_columns('adni_ucsf.txt');              
           % obj.modalities('ucsf_bl') = obj.read_columns('adni_ucsf_bl.txt'); 
            
            obj.modalities('ucsffsl') = obj.read_columns('ucsffsl.txt');
            obj.modalities('ucsffsx') = obj.read_columns('ucsffsx.txt');
            
            obj.modalities('AV45_misc_jt') = obj.read_columns('AV45_misc_jt.txt');
            obj.modalities('AV145') = obj.read_columns('AV145.txt');   
            
            obj.modalities('dti') = obj.read_columns('dti.txt');   
            
            obj.modalities('csf') = obj.read_columns('csf.txt');             
            
            obj.modalities('patient') = obj.read_columns('patient_details_limited.txt');            

        end
        
        function columns = read_columns(obj, file_name)
%       Imports the columns. Saves the trouble of strcat for each modality.

             columns = importdata(strcat( obj.MODALITY_PATH,file_name));
        end
        
        function obj = setDataset(obj, d_values )
%       Sets the current dataset type.
%       Attributes:
%           d_values: 1x2 cell array defining the combination of d1 and d2
%           values used.        
            obj.data_file = readtable( strcat( obj.MODALITY_PATH,obj.FILE_NAME));            
            obj.d_values = d_values;
            % obj.d_values_vec = obj.data_file(:,obj.D2_COLUMN); 
            %dataset_indexes = ismember( table2cell(obj.data_file(:,obj.D1_COLUMN)),  obj.d_values{1}) &  ismember( table2cell(obj.data_file(:,obj.D2_COLUMN)),  obj.d_values{2});         
            % obj.data_file = obj.data_file(dataset_indexes,:);  
        end
        
        
        function [measurements, time_points, time_diff] = getMeasurements(obj,ptid, modality, numeric_values)
%       Returns the measurements and corresponding months from baseline for 
%       a specific patient from the fiven modality.
%       Attributes:
%           ptid: PTID
%           modality: Name of the modality
%           numeric_values: Will the values be casted to double. OPTIONAL (Default: True)
%           
%       Output:
%           measurements: Timepoints x Features array of measurements.
%           time_points: Timepoints x 1 array of months from baseline.
%           time_diff: difference to beginning of 2018 in months
            if nargin < 4
               numeric_values = true; 
            end
            column_names = obj.modalities(modality);
            
            column_indexes = find(ismember(obj.data_file.Properties.VariableNames, column_names));
            measurements = obj.data_file(:,column_indexes);
            patient_indexes = ismember( table2cell(obj.data_file(:,obj.PTID_COLUMN)),  ptid  );
            measurements = measurements(patient_indexes,:);
            time_points = obj.data_file(patient_indexes,'M');
            %Convert time points to double
            time_points = cellfun( @(x) str2double(x),table2array(time_points));
            
            % compute the difference to time zero (1/1/2018)
            exdate = table2array(obj.data_file(patient_indexes,obj.EXAMDATE_COLUMN));
            exyear = cellfun(@(x) str2double(x(1:4)),exdate);
            exmon  = cellfun(@(x) str2double(x(6:7)),exdate);
            time_diff = (2018 - exyear - 1)*12 + (12 - exmon);
            
            %Convert numeric features to double
            if numeric_values
                numeric_measurements = cellfun( @(x) str2double(x),table2array(measurements));
                measurements = numeric_measurements;
               %Filter empty rows
               non_nan_indexes = ~all( isnan(measurements),2);
               measurements = measurements( non_nan_indexes,:);
               time_points = time_points(non_nan_indexes,:);
               time_diff =  time_diff(non_nan_indexes,:); 
            else
               %Filter empty rows
               measurements = table2array(measurements);
               non_empty_indexes = ~any(~cellfun(@(x) length(x) > 0, measurements),2);
               measurements = measurements( non_empty_indexes,:);
               time_points = time_points(non_empty_indexes,:);  
               time_diff =  time_diff(non_empty_indexes,:); 
            end
            
            %Sort measurements by time
            [time_points I ] = sort(time_points);
            measurements = measurements(I,:);
            time_diff = time_diff(I,:);
        end
        
        function [ids,D2,idsnum] = getIds(obj)
%       Returns the number of unique PTIDs and information whether the
%       subject is in D2
            [ids,ia] = unique(obj.data_file(:, obj.PTID_COLUMN));
            ids = table2array(ids);
            idsnum = zeros(length(ids),1);
            for i = 1:length(ids)
                idsnum(i) = str2num(ids{i}(7:end));
            end
            D2 = obj.data_file(:,obj.D2_COLUMN);
            D2 = cellfun( @(x) str2double(x),table2array(D2));
            D2 = D2(ia);
        end
        
        function n_features = getFeatureN(obj, modality)
%       Returns the number of features in the modality.
%       Attributes:
%           modality: Name of the modality.
%       Output:
%           n_features: The number of features in the modality.
            column_names = obj.modalities(modality);
            n_features = size(column_names,1);
        end       
        
        function valid_patients = getValidPatients(obj, modality, threshold, numeric)
%       Returns the PTIDs of patients who have measurements from the
%       modality exceeding the given threshold.
%       Attributes:
%           modality: Name of the modality.
%           threshold: minimun number of time points required OPTIONAL (Default: 0)
%           numeric: Is the data expected to be numeric. OPTIONAL (Default: True)
%       Output:
%           valid_patients: PTIDs with valid measurements.
            if nargin < 3
               threshold = 0; 
            end

            if nargin < 4
               numeric = true; 
            end            
            
            ids = obj.getIds();
            n_ids = size(ids,1);
            valid_ids = zeros(n_ids,1);
            column_names = obj.modalities(modality);    
            
            column_indexes = find(ismember(obj.data_file.Properties.VariableNames, column_names));
            measurements = table2array(obj.data_file(:,column_indexes));
            
            for i = 1:n_ids
                
                patient_indexes = ismember( table2cell(obj.data_file(:,obj.PTID_COLUMN)), ids(i,1)  );
                patient_measurements = measurements(patient_indexes,:);

                if numeric
                    numeric_measurements = cellfun( @(x) str2double(x),patient_measurements);    
                    non_nan_indexes = ~any( isnan(numeric_measurements),2);
                else             
                   non_nan_indexes = ~any( ~cellfun(@(x) length(x) > 0, patient_measurements),2);
                end
               
               if ~isempty(find(non_nan_indexes)) && size(find(non_nan_indexes),1) > threshold
                   valid_ids(i,1) = 1;
               end
                
            end
          
            valid_patients = ids(find(valid_ids(:,1)),:);

        end
   
        function values = getFirstOrLast(obj, modality, ptids, last, numeric )
%       Returns the first or last measurement of the series from  a given
%       modality for each PTID listed.
%       Attributes:
%           modality: Name of the modality.
%           ptids: List of PTIDs. OPTIONAL (Default: All ptids with more than zero measurements from the modality.)
%           last: Will the last measurement be retrieved. OPTIONAL (Default: True)
%       Optional:
%           values: Patients x Features array.
            if isempty(ptids)
                ptids = obj.getValidPatients( modality, 0, numeric);
            end
            
            valid_patients = ptids;
            
            if nargin < 3
               last = true;
            end 
            if nargin < 4
               numeric = true;
            end
           
            n_features = getFeatureN(obj, modality);   
            

            n_patients = size(valid_patients, 1);
            
            values = cell(n_patients,n_features);
            
            for i = 1:n_patients
                [measurements time_points] = getMeasurements(obj,valid_patients(i), modality, numeric);
                          
                if last

                    values(i,:) = num2cell(measurements(end,:));
                else
                    values(i,:) = measurements(1,:);
                end
            end
            
            if numeric
                values = cell2mat(values);
            end
            
           
        end        
        
    end

end