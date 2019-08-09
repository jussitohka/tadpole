% preparing the data for machine learning
clear
obj = TADPOLEFile; % D1 and D2
% obj2 = TADPOLEFile({'1','0'}); % only D1



[ids,d2flag,rids] = obj.getIds;
% ids2 = obj2.getIds;
% for i = 1:length(ids)
 %   idsnum(i) = str2num(ids{i}(7:end));
% end

% for i = 1:length(ids2)
 %   idsnum2(i) = str2num(ids2{i}(7:end));
%end
%d1flag = ismember(idsnum,idsnum2);
% d2flag = ~d1flag;
nsubj = length(ids);
nd2subj =  sum(double(d2flag));
mods = {'cognitive'
        'moca_ecog'
        'ucsf'
        'ucsffsl'
        'ucsffsx'
        'AV45_misc_jt'
        'AV145'
        'dti'
        'csf'
        'patient'};
    
    
ntp = zeros(nsubj,length(mods));    
tpmax = zeros(nsubj,length(mods));    

if 0
for i = 1:nsubj
    i
    for j = 1:length(mods)
        [measurements, time_points,time_diff] = getMeasurements(obj,ids{i}, mods{j});
         if isempty(time_points)
             ntp(i,j) = 0;
             tpmax(i,j) = 0;
             td(i,j) = 0;
         else
             ntp(i,j) = length(time_points);
             tpmax(i,j) = max(time_points);
             tdmax(i,j) = max(time_diff);
             tdmin(i,j) = min(time_diff);
         end
    end
end
end

if 0
    j = 10;
    for i = 1:nsubj
        i
        [measurements, time_points,time_diff] = getMeasurements(obj,ids{i}, mods{j});
        if length(time_points) > 2
            final_diff(i) = time_diff(end -1) - time_diff(end);
        else
            final_diff(i) = 0;
        end
    end
end
