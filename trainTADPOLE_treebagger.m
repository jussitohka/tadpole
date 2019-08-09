% train TreeBagger

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

delay_time = [24 24 24 24];

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

load('foldid')
% stage124 = load('stage1results_longitudinal24');
% stage112 = load('stage1results_longitudinal12');
stage10 = load('stage1results_longitudinal0');

if 1
% compute estimates using longitudinal data for each subject
indOI = [ADAS13idx ventricleICVidx];
fmeas = zeros(nsubj,2,2);
fR = zeros(nsubj,2);
fError = zeros(nsubj,2,3);
flag = 0;
tpthr = [7 4];
con = [0 0.5];

for i = 1:3
    mean_values{1,i} = zeros(1,6);
    mean_values{2,i} = zeros(1,8);
    nnn{1,i} = 0;
    nnn{2,i} = 0;
end

for j = 1:2
   for i = 1:nsubj
        if length(time_points{i,j}) > 1
            X = measurements{i,j}((end - 1),:);
            
            t0 = time_points{i,2}(end);
            inddiag = find(time_points{i,4} == t0);
            comdiag = measurements{i,diagnosismod}(inddiag,compactdiagnosisidx);
            if comdiag > 0 & ~any(isnan(X))
            
            mean_values{j,comdiag} = mean_values{j,comdiag} + X;
            nnn{j,comdiag} = nnn{j,comdiag} + 1;
            end
            X = [X ones(size(X,1),1)]; 
            ttt = time_points{i,j}(end) - time_points{i,j}(end - 1);
            ttt = round((ttt - con(j))/12);
            if ttt > 4
                ttt = 4;
            end
            if ttt == 0
                ttt = 1;
                flag = 1;
            else
                flag = 0;
            end
            ttt = 2*ttt -1; % don't use 6 month models
            fmeas(i,2,j) = X*stage10.bbb{j,ttt,j}{foldid{1}(i)};
          %  fmeas(i,2,j) = max(fmeas(i,2,j),0);
            fError(i,j,3) = abs(fmeas(i,2,j) -  measurements{i,j}(end,indOI(j)));
        end
        if length(time_points{i,j}) > 2 
            [fmeas(i,1,j),fR(i,j)] = longitudinal_estimate(measurements{i,j}(1:(end - 1),indOI(j)),time_points{i,j}(1:(end - 1)),time_points{i,j}(end));
           %  fmeas(i,1,j) = max(fmeas(i,1,j),0);
            fError(i,j,1) = abs(fmeas(i,1,j) -  measurements{i,j}(end,indOI(j)));
             tll = (4/length(time_points{i,j}));
             if ~flag
                 if length(time_points{i,j}) > tpthr(j)
                     alpha = [fR(i,j)^tll 1 - fR(i,j)^tll];
                 else
                     alpha = [0.1 0.9];
                 end
                 fError(i,j,2) = abs(sum(alpha.*fmeas(i,:,j)) -  measurements{i,j}(end,indOI(j)));
                 
             else
                 alpha = [1 0];
                 fError(i,j,2) = abs(sum(alpha.*fmeas(i,1:2,j)) -  measurements{i,j}(end,indOI(j)));
             end
        elseif length(time_points{i,j}) == 2  
            fmeas(i,1,j) = measurements{i,j}(1,indOI(j));
            fR(i,j) = 0.25;
            fError(i,j,1) = abs(fmeas(i,1,j) -  measurements{i,j}(end,indOI(j)));
            if ~flag
                 alpha = [0 1];  
                 fError(i,j,2) = abs(sum(alpha.*fmeas(i,1:2,j)) -  measurements{i,j}(end,indOI(j)));
             else
                 alpha = [1 0];
                 fError(i,j,2) = abs(sum(alpha.*fmeas(i,1:2,j)) -  measurements{i,j}(end,indOI(j)));
             end
        else
            fError(i,j,1:3) = NaN;          
        end
   end
end

ci = nanmedian(fError);
for j = 1:2
    for i = 1:3
        mean_values{j,i} = mean_values{j,i}/nnn{j,i};
    end
end


end
%% classification problems

% Yd = [];
% Xd = [];
% for i = 1:nsubj
%     
%     % mri measurement is the key
%     t0 = time_points{i,2}(end);
%     indadas = find(time_points{i,1} == t0);
%     inddiag = find(time_points{i,4} == t0);
%     if ~isempty(inddiag) & ~isempty(indadas)
%         comdiag = measurements{i,diagnosismod}(inddiag,compactdiagnosisidx);
%         fmri = measurements{i,2}(end,indOI(2));
%         fadas = measurements{i,1}(indadas,indOI(1));
%         % first nc vs. mci
%         if comdiag == 1 | comdiag == 2 | comdiag == 3
%             if fmri > 0 & fadas > 0
%                 Yd = [Yd;comdiag];
%                 Xd = [Xd;[fadas fmri]];
%             end
%         end
%     end
% end

        





% First ADAS13 
% s = 1;
% Yhat = stage10.Yhat;
% Ytrue = stage10.Ytrue;
% for k = 1:2:length(stage10.target_predictions) % predictions at 12 months at 24 months, at 36 months at 48 months
%     % first only using np 
%     cls = 0; % regression
%     complete_idx = find(trainidx{s}(:,k) & testlabel{s}(:,k,s));
%    % X = [Yhat{s,k}(complete_idx,s) fmeas(complete_idx,s),fR(complete_idx,s)];
%     sz = 6;
%     X1 = zeros(length(complete_idx),sz);
%     for i = 1:length(complete_idx)
%         X1(i,1:sz) = measurements{complete_idx(i),s}(trainidx{s}(complete_idx(i),k),:);
%     end    
%     [tbmae(k,s,1),tbcr(k,s,1),tb{k,s,1}] = cv_treebagger(X1,Ytrue{s,k}(complete_idx,s),foldid{1}(complete_idx),cls);
%     % then using np and mri
%     complete_idx = find(trainidx{s}(:,k) & trainidx{2}(:,k) & testlabel{s}(:,k,s));
%     sz2 = 8;
%     X2 = zeros(length(complete_idx),sz2);
%     X1 = zeros(length(complete_idx),sz);
%     for i = 1:length(complete_idx)
%         X1(i,1:sz) = measurements{complete_idx(i),s}(trainidx{s}(complete_idx(i),k),:);
%         X2(i,1:sz2) = measurements{complete_idx(i),2}(trainidx{2}(complete_idx(i),k),:);
%     end 
%   %  X = [Yhat{s,k}(complete_idx,[1 2]) fmeas(complete_idx,s),fR(complete_idx,s)];
%     X = [X1 X2];
%     [tbmae(k,s,2),tbcr(k,s,2),tb{k,s,2}] = cv_treebagger(X,Ytrue{s,k}(complete_idx,s),foldid{1}(complete_idx),cls);
% end
% % For MRI
% s = 2;
% for k = 1:2:length(stage10.target_predictions) % predictions at 12 months at 24 months, at 36 months at 48 months
%     % only using MRI 
%     cls = 0; % regression
%     complete_idx = find(trainidx{s}(:,k) & testlabel{s}(:,k,s));
%    %  X = [Yhat{s,k}(complete_idx,s) fmeas(complete_idx,s),fR(complete_idx,s)];
%     X2 = zeros(length(complete_idx),sz2);
%     for i = 1:length(complete_idx)
%         X2(i,1:sz2) = measurements{complete_idx(i),2}(trainidx{2}(complete_idx(i),k),:);
%     end 
%     [tbmae(k,s,1),tbcr(k,s,1),tb{k,s,1}] = cv_treebagger(X2,Ytrue{s,k}(complete_idx,s),foldid{1}(complete_idx),cls);
% end
% save results_treebagger_regression tbmae tbcr tb  

% evaluation 

if 0 
hippind = 2;
Yd = [];
Xd = [];
compind = [];
for i = 1:nsubj
    
    % mri measurement is the key
    if length(time_points{i,2}) > 2
        t0 = time_points{i,2}(end);
        t1 = time_points{i,2}(end - 1);
        indadas = find(time_points{i,1} == t0);
        inddiag = find(time_points{i,4} == t0);
        if ~isempty(inddiag) & ~isempty(indadas)
            comdiag = measurements{i,diagnosismod}(inddiag,compactdiagnosisidx);
            fmri = (measurements{i,2}(end,indOI(2))); % -   measurements{i,2}(end - 1,indOI(2)))/(t0 - t1);
            fadas = measurements{i,1}(indadas,indOI(1));
            % first nc vs. mci
            if comdiag == 1 | comdiag == 2 | comdiag == 3
                if ~isnan(fmri) & fadas > 0
                    Yd = [Yd;comdiag];
                    Xd = [Xd;[fadas fmri measurements{i,4}(inddiag,3) (measurements{i,4}(inddiag,2) + t0/12)]];
                    compind = [compind; i];
                end
            end
        end
    end
end
clsind = find(Yd < 3);
[mae(1),cr(1),tb{1}] = cv_treebagger(Xd(clsind,:),Yd(clsind) - 1,foldid{1}(compind(clsind)),1);
meth = 'classification';
tb_complete{1} = TreeBagger(500,Xd(clsind,:),Yd(clsind) - 1,'OOBPrediction','on','method',meth); 
clsind = find(Yd > 1);
[mae(2),cr(2),tb{2}] = cv_treebagger(Xd(clsind,:),Yd(clsind) - 2,foldid{1}(compind(clsind)),1);
tb_complete{1} = TreeBagger(500,Xd(clsind,:),Yd(clsind) - 2,'OOBPrediction','on','method',meth);

save tb_complete tb_complete     
end      
      
