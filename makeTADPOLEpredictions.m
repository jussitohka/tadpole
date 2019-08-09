% makeTADPOLEpredictions
% run prepareTADPOLEdata and trainTADPOLE treebagger before

D3 = readtable('./modalities/TADPOLE_D3.csv');
D3{  ismember(table2array(D3(:,2)),'bl') , 2} = {'m00'};
D3_rids = D3.RID;
D3_time_points =  D3.VISCODE;
% D3_time_points = cellfun( @(x) str2double(x),table2array(D3_time_points));
D3_time_points = cellfun( @(x) str2double(x(2:end)),D3_time_points);
D3_exdate = cellstr(table2array(D3(:,3))); 
D3_exyear = cellfun(@(x) str2double(x(1:4)),D3_exdate);
D3_exmon  = cellfun(@(x) str2double(x(6:7)),D3_exdate);
D3_time_diff = (2018 - D3_exyear - 1)*12 + (12 - D3_exmon);
nD3subj = length(D3_rids);

D3{  ismember(table2array(D3(:,4)),'NL') , 4} = {'1'};
D3{  ismember(table2array(D3(:,4)),'MCI') , 4} = {'2'};  
D3{  ismember(table2array(D3(:,4)),'MCI to Dementia') , 4} = {'3'};
D3{  ismember(table2array(D3(:,4)),'Dementia') , 4} = {'3'};
D3{  ismember(table2array(D3(:,4)),'NL to MCI') , 4} = {'2'};
D3diagnosis =  cellfun( @(x) str2double(x),D3.DX);
D3diagnosis(isnan(D3diagnosis)) = 2;
D3diagnosis(D3diagnosis == 0) = 2;
load mean_values
j = 1;
con = [0 0.5];
for i = 1:nD3subj
    X = [D3.ADAS13(i) D3.MMSE(i) mean_values{j,D3diagnosis(i)}(3:6)];
    X = [X ones(size(X,1),1)]; 
    isn = isnan(X);
    X(isn) = mean_values{j,D3diagnosis(i)}(isn);
    t0 = D3_time_diff(i);
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
            ttt = 2*ttt -1;
            D3pred(i,j,tt) = X*stage10.bbb{j,ttt,j}{1};
    end
end
j = 2;
for i = 1:nD3subj
    X = [D3.Ventricles(i) D3.Hippocampus(i) D3.WholeBrain(i) D3.Entorhinal(i) D3.Fusiform(i) D3.MidTemp(i) D3.ICV(i) D3.Ventricles(i)/D3.ICV(i)];
    X = [X ones(size(X,1),1)]; 
    isn = isnan(X);
    X(isn) = mean_values{j,D3diagnosis(i)}(isn);
    t0 = D3_time_diff(i);
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
        ttt = 2*ttt -1;
        D3pred(i,j,tt) = X*stage10.bbb{j,ttt,j}{1};
    end
end
load tb_complete
load conversion_probabilities
meth = 'classification';
ttvec = [1 12 24 36 48 60];
for i = 1:nD3subj
    disp(i)
    t0 = D3_time_diff(i);
    age = D3.AGE(i) + t0/12 + D3_time_points(i)/12;
    if D3diagnosis(i) == 1
        for tt = 1:length(ttvec)
            
            ttt = round(ttvec(tt) + t0);
            X = [D3pred(i,1,tt) D3pred(i,2,tt) 0 age + ttt/12];
            [yyy scores]  = predict(tb_complete{1},X);
            if ttt > 48
                ttt = 48;
            end
            D3prob_scores_small(i,1:2,tt) = scores.*[(1 - conversion_probability_nl_i(ttt,2)) conversion_probability_nl_i(ttt,2)];
            D3prob_scores_small(i,1:2,tt) = D3prob_scores_small(i,1:2,tt)/sum(D3prob_scores_small(i,1:2,tt));
            
        end
        D3prob_scores(i,1,1:60) = interp1(ttvec,squeeze(D3prob_scores_small(i,1,:)),1:60);
        D3prob_scores(i,2,1:60) = 1 - D3prob_scores(i,1,1:60); 
        D3prob_scores(i,3,1:60) = 0;
    end
    if D3diagnosis(i) == 2
        for tt = 1:length(ttvec)
            ttt = round(ttvec(tt) + t0);
            X = [D3pred(i,1,tt) D3pred(i,2,tt) 0 age + ttt/12];
            [yyy scores]  = predict(tb_complete{2},X);
            if ttt > 48
                ttt = 48;
            end
            D3prob_scores_small(i,2:3,tt) = scores.*[(1 - conversion_probability_ad_i(ttt,2)) conversion_probability_ad_i(ttt,2)];
            D3prob_scores_small(i,2:3,tt) = D3prob_scores_small(i,2:3,tt)/sum(D3prob_scores_small(i,2:3,tt));
            
        end
        D3prob_scores(i,3,1:60) = interp1(ttvec,squeeze(D3prob_scores_small(i,3,:)),1:60);
        D3prob_scores(i,2,1:60) = 1 - D3prob_scores(i,3,1:60); 
        D3prob_scores(i,1,1:60) = 0;
    end
    if D3diagnosis(i) == 3
        for tt = 1:60
            D3prob_scores(i,1,tt) = 0;
            D3prob_scores(i,2,tt) = 0;
            D3prob_scores(i,3,tt) = 1;
        end
    end
end
save D3predictions D3prob* D3pred D3_rids
        