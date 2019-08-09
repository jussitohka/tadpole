% run prepareTADPOLEdata first
% Now just cognition, MRI, and PET

clear measurements time_points time_diff
modality_indexes = [1 3 6 10];
ADAS13mod = 1;
ADAS13idx = 1;  

ventricleidx = 1;
ventriclemod = 2;
ICVidx = 7;
ventricleICVidx = 8;

diagnosisidx = 1;
diagnosismod = 4;
% read the data 
for i = 1:nsubj
    for j = 1:length(modality_indexes)
        [measurements{i,j}, time_points{i,j},time_diff{i,j}] = getMeasurements(obj,ids{i}, mods{modality_indexes(j)});         
    end
end
% compute ventricles divided by ICV
for i = 1:nsubj
    measurements{i,ventriclemod}(:,ventricleICVidx) =   measurements{i,ventriclemod}(:,ICVidx)./ measurements{i,ventriclemod}(:,ventricleidx); 
end    


% for k = 1:length(target_predictions)
for j = 1:length(modality_indexes)
    td{j} = zeros(nsubj,1);
    longi{j} = zeros(nsubj,1);
    interval_len{j} = zeros(nsubj,1);
    for i = 1:nsubj
        longi{j}(i) = length(time_points{i,j});
        if length(time_points{i,j}) > 1
            td{j}(i) = time_diff{i,j}(end -1) - time_diff{i,j}(end);            
            if length(time_points{i,j}) > 2
                interval_len{j}(i) = time_diff{i,j}(end -2) - time_diff{i,j}(end - 1);
            end
         end 
    end
end

if 0
for j = 1:5
    figure
    hist(td{j})
    figure
    hist(longi{j})
end
end


% prepare training data


lidx = [ADAS13idx,ventricleICVidx,diagnosisidx];
lmod = [ADAS13mod,ventriclemod,diagnosismod];

target_predictions = [12:6:48];
interval = 3;
ulim = target_predictions + interval;
dlim = target_predictions - interval;

for j = 1:length(modality_indexes)
    trainidx{j} = zeros(nsubj,length(target_predictions));
    testlabel{j} = zeros(nsubj,length(target_predictions),3);
    testlabelch{j} = zeros(nsubj,length(target_predictions),3);
    for i = 1:nsubj
        if length(time_points{i,j}) > 1
            ti{i,j} = time_diff{i,j}(1:(end - 1)) - time_diff{i,j}(end); 
            for k = 1:length(target_predictions)
                iii = find(ti{i,j} < ulim(k) & ti{i,j} > dlim(k));
                if ~isempty(iii)
                    trainidx{j}(i,k) = iii(end); %  the subject i is training material for the modality j and time interval k 
                                                 % the last index is suitable for training is trainidx{j}(i,k)
                end
                if trainidx{j}(i,k) %interval_len{j}(i) < ulim(k) & interval_len{j}(i) > dlim(k)                    
                    ttpp = time_points{i,j}(end); % this is the final time point for the subject i 
                    for s = 1:3
                        idx = find(time_points{i,lmod(s)} == ttpp);
                        if ~isempty(idx) 
                            if ~isnan(measurements{i,lmod(s)}(idx,lidx(s)))
                                testlabel{j}(i,k,s) = measurements{i,lmod(s)}(idx,lidx(s));
                                
                                % change this for rate of change (per month)
                                % if a regression problem
                                if s < 3 % s == 3 classification problem
                                    if idx > 1
                                       testlabelch{j}(i,k,s) = measurements{i,lmod(s)}(idx - 1,lidx(s)) - testlabel{j}(i,k,s);
                                       % in months
                                       testlabelch{j}(i,k,s) = testlabelch{j}(i,k,s)/(time_points{i,lmod(s)}(idx - 1) - time_points{i,lmod(s)}(idx));
                                    end
                                end
                            end
                        end
                    end                                                   
                end
            end              
        end
    end
end

for i = 1:10
    foldid{i} = balanced_crossval([1:nsubj],10,[],0,1);
end

%%%%%%%%%%%%%%%%%%
% Cross-sectional training
%%%%%%%%%%%%%%%%%%



for j = 1:length(modality_indexes)
    for k = 1:length(target_predictions)
        for s = 1:2
           complete_idx = find(trainidx{j}(:,k) & testlabel{j}(:,k,s));
           lls(j,k,s) =  length(complete_idx);
           if lls(j,k,s) > 20
           X = zeros(length(complete_idx),size(measurements{complete_idx(1),j},2));
           Y = zeros(length(complete_idx),1);
           for i = 1:length(complete_idx)
               X(i,:) = measurements{complete_idx(i),j}(trainidx{j}(complete_idx(i),k),:);
               Y(i) = testlabel{j}(complete_idx(i),k,s);
           end
           [mae(j,k,s),cr(j,k,s)] = cv_regression(X,Y,foldid{1}(complete_idx));
           if s == 1 & j == 1
              priormae(j,k,s) = nanmean(abs(Y - X(:,ADAS13idx)));
           elseif s == 2 & j == 2
               priormae(j,k,s) = nanmean(abs(Y - X(:,ventricleICVidx)));
           else
               priormae(j,k,s) = 0;
           end    
           else     
              mae(j,k,s) = 0;
              cr(j,k,s) = 0;
           end
        end
    end
end

           
           
            
              
              


