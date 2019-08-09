function [deltas fitted ] = calculate_temporal_changes( tadpole_file, modality, fit_function)
%Calculates the changes over time for a selected variable of a specific
%modality.

%Inputs:
%   id_data: Table presenting the data related to each individual in the
%            study. This has to contain at least PTID in column 2 and VISCODE in
%            column 3
%   modality: A table of measurement data for a given modality.
%   fit_function: Function to be used to fit a line through the
%           measurements. Must take unevenly spaced time points and the
%           measurements values as an input and return a vector of coefficients.
%
%Outputs:
%   deltas: Changes in the values between visits
%   fitted: Line fitted through the selected variable

    if isempty(fit_function)
       fit_function = @fit_line; 
    end
    
    unique_ids = tadpole_file.getValidPatients(modality,1,true); 
    
    n_patients = size(unique_ids,1); 
    n_features = tadpole_file.getFeatureN(modality);
    
    deltas = cell(n_patients, n_features);
    fitted = cell(n_patients, n_features);
    
    for id = 1:size(unique_ids,1)
        
        if mod(id,100) == 0
            disp(sprintf('ID %i',id));    
        end

        [ measurements time_points ] = tadpole_file.getMeasurements(unique_ids(id),modality);
        
        for f = 1:n_features    
            
            [p]= fit_function(time_points,measurements(:,f));
            fitted{id,f}  = p; 
            
            deltas{id,f} = diff(measurements(:,f));
            
        end
    end  
    
    function slope = fit_line(x,y)
        [p S] = polyfit(x,y,1);
        slope =p(1);  
    end


end