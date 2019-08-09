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




% prepare training data


lidx = [ADAS13idx,ventricleICVidx,diagnosisidx];
lmod = [ADAS13mod,ventriclemod,diagnosismod];

target_predictions = [12:6:48];
interval = 3;
ulim = target_predictions + interval;
dlim = target_predictions - interval;

delay_time = [24 24 24 24];



load('foldid')
% stage124 = load('stage1results_longitudinal24');
% stage112 = load('stage1results_longitudinal12');
stage10 = load('stage1results_longitudinal0');


% compute estimates using longitudinal data for each subject
indOI = [ADAS13idx ventricleICVidx];


flag = 0;
tpthr = [7 4];
con = [0 0.5];

load mean_values
srids = sort(rids(d2flag == 1));
fmeas = zeros(60,2);
fR = zeros(60,1);
for j = 1:2
    disp(j)
   for kkk = 1:length(srids)
       % reogranize the data so that rids are ordered
       i = find(rids == srids(kkk));
       X = measurements{i,j}(end,:);
       t0 = time_points{i,j}(end);
       inddiag = find(time_points{i,4} == t0);
       comdiag = measurements{i,diagnosismod}(inddiag,compactdiagnosisidx);
       X = [X ones(size(X,1),1)]; 
       isn = isnan(X) | X == 0;
       D2diagnosis = comdiag;
       if isnan(D2diagnosis) | D2diagnosis == 0
           D2diagnosis = 2;
       end
       X(isn) = mean_values{j,D2diagnosis}(isn);
       t0 = time_diff{i,j}(end);     
       for tt = 1:60
           ttt = tt + t0;
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
           fmeas(tt,2) = X*stage10.bbb{j,ttt,j}{1};
           if length(time_points{i,j}) > 1 
              [fmeas(tt,1),fR(tt)] = longitudinal_estimate(measurements{i,j}(:,indOI(j)),time_points{i,j},t0 + tt);
           
              tll = (4/length(time_points{i,j}));
              if ~flag
                 if length(time_points{i,j}) > tpthr(j)
                     alpha = [fR(tt)^tll 1 - fR(tt)^tll];
                 else
                     alpha = [0.1 0.9];
                 end   
                 D2pred(kkk,j,tt) = sum(alpha.*fmeas(tt,:));
              else
                 alpha = [1 0];
                 D2pred(kkk,j,tt) = sum(alpha.*fmeas(tt,:));
              end
           elseif length(time_points{i,j}) == 1  
               fmeas(tt,1) = measurements{i,j}(1,indOI(j));
               fR(tt) = 0.25;
               if ~flag
                 alpha = [0 1];  
                 D2pred(kkk,j,tt) = sum(alpha.*fmeas(tt,1:2));
               else
                 alpha = [1 0];
                 D2pred(kkk,j,tt) = sum(alpha.*fmeas(tt,1:2)); 
               end
           else
               fmeas(tt,1) = 0;
               alpha = [0 1];
               D2pred(kkk,j,tt) = sum(alpha.*fmeas(tt,1:2)); 
               
           end
           if isnan(D2pred(kkk,j,tt))
               if ~isnan(fmeas(tt,1))
                   D2pred(kkk,j,tt) = fmeas(tt,1);
               elseif ~isnan(fmeas(tt,2))
                   D2pred(kkk,j,tt) = fmeas(tt,2);
               else
                 D2pred(kkk,j,tt) = mean_values{j,D2diagnosis}(indOI(j));
               end
           end
        end
        
   end
end
D2pred = real(D2pred);
D2pred = max(D2pred,0);
load tb_complete
load conversion_probabilities
meth = 'classification';
ttvec = [1 12 24 36 48 60];
for kkk = 1:length(srids)
   disp(kkk)
   i = find(rids == srids(kkk));
   t0 = time_diff{i,4}(end);
   age = (measurements{i,4}(end,2));
   apo =  measurements{i,4}(end,3);
   D2diagnosis = measurements{i,4}(end,4);
   if isnan(D2diagnosis) | D2diagnosis == 0
       D2diagnosis = 2;
   end
   if D2diagnosis == 1
       
       for tt = 1:length(ttvec)
           ttt = round(ttvec(tt) + t0);
           X = [D2pred(kkk,1,tt) D2pred(kkk,2,tt) apo age + ttt/12];
           [yyy scores]  = predict(tb_complete{1},X);
           if ttt > 48
                ttt = 48;
           end
           D2prob_scores_small(kkk,1:2,tt) = scores.*[(1 - conversion_probability_nl_i(ttt,2)) conversion_probability_nl_i(ttt,2)];
           D2prob_scores_small(kkk,1:2,tt) = D2prob_scores_small(kkk,1:2,tt)/sum(D2prob_scores_small(kkk,1:2,tt));
       end
       D2prob_scores(kkk,1,1:60) = interp1(ttvec,squeeze(D2prob_scores_small(kkk,1,:)),1:60);
       D2prob_scores(kkk,2,1:60) = 1 - D2prob_scores(kkk,1,1:60); 
       D2prob_scores(kkk,3,1:60) = 0;
       
   end
   if D2diagnosis == 2
       for tt = 1:length(ttvec)
           ttt = round(ttvec(tt) + t0);
           X = [D2pred(kkk,1,tt) D2pred(kkk,2,tt) apo age + ttt/12];
           [yyy scores]  = predict(tb_complete{2},X);
           if ttt > 48
                ttt = 48;
           end
           D2prob_scores_small(kkk,2:3,tt) = scores.*[(1 - conversion_probability_ad_i(ttt,2)) conversion_probability_ad_i(ttt,2)];
           D2prob_scores_small(kkk,2:3,tt) = D2prob_scores_small(kkk,2:3,tt)/sum(D2prob_scores_small(kkk,2:3,tt));
       end
       D2prob_scores(kkk,3,1:60) = interp1(ttvec,squeeze(D2prob_scores_small(kkk,3,:)),1:60);
       D2prob_scores(kkk,2,1:60) = 1 - D2prob_scores(kkk,3,1:60); 
       D2prob_scores(kkk,1,1:60) = 0;
   end
   if D2diagnosis == 3
        for tt = 1:60
            D2prob_scores(kkk,1,tt) = 0;
            D2prob_scores(kkk,2,tt) = 0;
            D2prob_scores(kkk,3,tt) = 1;
        end
    end
end
save D2predictions D2prob* D3pred srids
%% classification problems


        







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
      
