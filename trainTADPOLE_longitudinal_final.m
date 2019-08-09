% run prepareTADPOLEdata first
% Now just cognition, MRI, and PET
% Training cross-sectional tree-bagger models as have noted
% that longitudinal info contributes essentially nothing. 
% Integrating longitudinal one variable models into these models

clear measurements time_points time_diff td
modality_indexes = [1 3 6 10];
ADAS13mod = 1;
ADAS13idx = 1;  

ventricleidx = 1;
ventriclemod = 2;
ICVidx = 7;
ventricleICVidx = 8;
% longitudinal = 1;

diagnosisidx = 1;
diagnosismod = 4;
compactdiagnosisidx = 4;
% read the data 
for i = 1:nsubj
    for j = 1:length(modality_indexes)
        [measurements{i,j}, time_points{i,j},time_diff{i,j}] = getMeasurements(obj,ids{i}, mods{modality_indexes(j)});         
    end
end
% compute ventricles divided by ICV
% and compact diagnosis 
 % 1=Stable:NL to NL, 
        % 2=Stable:MCI to MCI, 
        % 3=Stable:AD to AD, 
        % 4=Conv:NL to MCI, 
        % 5=Conv:MCI to AD, 
        % 6=Conv:NL to AD, 
        % 7=Rev:MCI to NL, 
        % 8=Rev:AD to MCI, 
        % 9=Rev:AD to NL, -1=Not available
for i = 1:nsubj
    measurements{i,ventriclemod}(:,ventricleICVidx) =   measurements{i,ventriclemod}(:,ventricleidx)./ measurements{i,ventriclemod}(:,ICVidx); 
    xxx = measurements{i,diagnosismod}(:,diagnosisidx);
    measurements{i,diagnosismod}(:,compactdiagnosisidx) = ((xxx == 1) | (xxx == 7) | (xxx == 9))*1 +... 
        2*((xxx == 2) | (xxx == 4) | (xxx == 8)) + 3*((xxx == 3) | (xxx == 5) | (xxx == 6));
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

% delay_time = [24 24 24 24];

for j = 1:length(modality_indexes)
    trainidx{j} = zeros(nsubj,length(target_predictions));
    trainidx2{j} = zeros(nsubj,length(target_predictions));
    testlabel{j} = zeros(nsubj,length(target_predictions),3);
    testlabelch{j} = zeros(nsubj,length(target_predictions),3);
    trainlabel{j} = zeros(nsubj,length(target_predictions),3);
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
                    s = 3;
                    trainlabel{j}(i,k,s) = measurements{i,lmod(s)}(trainidx{j}(i,k),lidx(s)); % this to simplify the classification 
                    ttpp = time_points{i,j}(end); % this is the final time point for the subject i 
                    idx = find(time_points{i,j} == (time_points{i,j}(trainidx{j}(i,k)) - delay_time(j)));
                    if length(idx) == 1
                        trainidx2{j}(i,k) = idx;
                    end
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

if exist('foldid','file')
    load('foldid')
else
    for i = 1:10
        foldid{i} = balanced_crossval([1:nsubj],10,[],0,1);        
    end
    save('foldid','foldid');
end

% compute estimates using longitudinal data for each subject
indOI = [ADAS13idx ventricleICVidx];
fmeas = zeros(nsubj,2);
fR = zeros(nsubj,2);
fError = zeros(nsubj,2);
for j = 1:2
   for i = 1:nsubj
        if length(time_points{i,j}) > 2 
            [fmeas(i,j),fR(i,j)] = longitudinal_estimate(measurements{i,j}(1:(end - 1),indOI(j)),time_points{i,j}(1:(end - 1)),time_points{i,j}(end));
             fError(i,j) = abs(fmeas(i,j) -  measurements{i,j}(end,indOI(j)));
        elseif length(time_points{i,j}) == 2  
            fmeas(i,j) = measurements{i,j}(1,indOI(j));
            fR(i,j) = 0.25;
            fError(i,j) = abs(fmeas(i,j) -  measurements{i,j}(end,indOI(j)));
        else
            fError(i,j) = NaN;          
        end
   end
end

%%%%%%%%%%%%%%%%%%
% Cross-sectional training
%%%%%%%%%%%%%%%%%%


for k = 1:length(target_predictions)
   for s = 1:4
       Yhat{s,k} = zeros(nsubj,length(modality_indexes));
       Ytrue{s,k} = zeros(nsubj,length(modality_indexes));
   end
    
   
    for s = 1:2
        for j = 1:length(modality_indexes)
           complete_idx = find(trainidx{j}(:,k) & testlabel{j}(:,k,s) & trainidx2{j}(:,k));
           lls(j,k,s) =  length(complete_idx);
           if lls(j,k,s) > 20
               sz = size(measurements{complete_idx(1),j},2);
               if longitudinal 
                   X = zeros(length(complete_idx),2*sz);
               else
                   X = zeros(length(complete_idx),sz);
               end
               
               Y = zeros(length(complete_idx),1);
               for i = 1:length(complete_idx)
                   X(i,1:sz) = measurements{complete_idx(i),j}(trainidx{j}(complete_idx(i),k),:);
                   if longitudinal 
                       tt2 = trainidx2{j}(complete_idx(i),k);
                       tt = trainidx{j}(complete_idx(i),k);
                       tdd = time_points{complete_idx(i),j}(tt2:(tt - 1)) - time_points{complete_idx(i),j}(tt);

                       X(i,(sz +1):(2*sz)) = nanmean((measurements{complete_idx(i),j}(tt2:(tt - 1),:) - measurements{complete_idx(i),j}(tt,:))./tdd);
                   end
                   Y(i) = testlabel{j}(complete_idx(i),k,s);
               end
               [mae(j,k,s),cr(j,k,s),bbb{j,k,s},Yhat{s,k}(complete_idx,j)] = cv_regression(X,Y,foldid{1}(complete_idx));
               
               Ytrue{s,k}(complete_idx,j) = Y;
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
         X = [Yhat{s,k}(complete_idx,:) fmeas(complete_idx,j),fR(i,j)
        [tbmae(k,s),cr(k,s),tb{k,s},tbYhat{s,k}] = cv_treebagger(X,Ytrue{s,k}(complete_idx,j),foldid{1}(complete_idx,cls);
       %  [tbmae(k,s),cr(k,s),tb{k,s},tbYhat{s,k}] = cv_treebagger(Yhat{s,k}(complete_idx,:),Y,foldid,cls);
    end
    s = 3;
        % 1=Stable:NL to NL, 
        % 2=Stable:MCI to MCI, 
        % 3=Stable:AD to AD, 
        % 4=Conv:NL to MCI, 
        % 5=Conv:MCI to AD, 
        % 6=Conv:NL to AD, 
        % 7=Rev:MCI to NL, 
        % 8=Rev:AD to MCI, 
        % 9=Rev:AD to NL, -1=Not available
        %*****************************
        % For NC -> MCI
        % *****************************
    for j = 1:length(modality_indexes)   
        diaglog1 = (testlabel{j}(:,k,s) == 2) | (testlabel{j}(:,k,s) == 4) | (testlabel{j}(:,k,s) == 1); 
        diaglog2 = (trainlabel{j}(:,k,s) == 1) | (trainlabel{j}(:,k,s) == 7);
        complete_idx = find(trainidx{j}(:,k) &  trainidx2{j}(:,k) & diaglog1 & diaglog2);
        lls(j,k,s) =  length(complete_idx);
        if lls(j,k,s) > 20
           sz = size(measurements{complete_idx(1),j},2);
           if longitudinal 
               X = zeros(length(complete_idx),2*sz);
           else
               X = zeros(length(complete_idx),sz);
           end
           
           Y = zeros(length(complete_idx),1);
           for i = 1:length(complete_idx)
               X(i,1:sz) = measurements{complete_idx(i),j}(trainidx{j}(complete_idx(i),k),:);
               if longitudinal 
                   tt2 = trainidx2{j}(complete_idx(i),k);
                   tt = trainidx{j}(complete_idx(i),k);
                   tdd = time_points{complete_idx(i),j}(tt2:(tt - 1)) - time_points{complete_idx(i),j}(tt);

                   X(i,(sz +1):(2*sz)) = nanmean((measurements{complete_idx(i),j}(tt2:(tt - 1),:) - measurements{complete_idx(i),j}(tt,:))./tdd); 
               end
               Y(i) = testlabel{j}(complete_idx(i),k,s) == 1; % label 1 NC, label 0 MCI
           end
           [mae(j,k,s),cr(j,k,s),bbb{j,k,s},Yhat{s,k}(complete_idx,j)] = cv_regression(X,Y,foldid{1}(complete_idx),1);
           Ytrue{s,k}(complete_idx,j) = Y;
        else     
            mae(j,k,s) = 0;
            cr(j,k,s) = 0;
        end
    end
   
        %*****************************
        % For MCI -> AD
        % *****************************
        diaglog1 = (testlabel{j}(:,k,s) == 3) | (testlabel{j}(:,k,s) == 5) | (testlabel{j}(:,k,s) == 2); 
        diaglog2 = (trainlabel{j}(:,k,s) == 2) | (trainlabel{j}(:,k,s) == 4);
        complete_idx = find(trainidx{j}(:,k) &  trainidx2{j}(:,k) & diaglog1 & diaglog2);
        lls(j,k,s + 1) =  length(complete_idx);
        if lls(j,k,s + 1) > 20
           sz = size(measurements{complete_idx(1),j},2);
           if longitudinal 
               X = zeros(length(complete_idx),2*sz);
           else
               X = zeros(length(complete_idx),sz);
           end
         
           Y = zeros(length(complete_idx),1);
           for i = 1:length(complete_idx)
               X(i,1:sz) = measurements{complete_idx(i),j}(trainidx{j}(complete_idx(i),k),:);
               if longitudinal 
                   tt2 = trainidx2{j}(complete_idx(i),k);
                   tt = trainidx{j}(complete_idx(i),k);
                   tdd = time_points{complete_idx(i),j}(tt2:(tt - 1)) - time_points{complete_idx(i),j}(tt);

                   X(i,(sz +1):(2*sz)) = nanmean((measurements{complete_idx(i),j}(tt2:(tt - 1),:) - measurements{complete_idx(i),j}(tt,:))./tdd); 
               end
               Y(i) = testlabel{j}(complete_idx(i),k,s) == 2; % label 1 MCI, label 0 AD
           end
           [mae(j,k,s + 1),cr(j,k,s + 1),bbb{j,k,s + 1},Yhat{s + 1,k}(complete_idx,j)] = cv_regression(X,Y,foldid{1}(complete_idx),1);
           Ytrue{s + 1,k}(complete_idx,j) = Y;
        else     
            mae(j,k,s + 1) = 0;
            cr(j,k,s + 1) = 0;
        end
        
        
        
    end
end

           
           
            
              
              


