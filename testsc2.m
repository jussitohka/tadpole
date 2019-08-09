%% classification problems
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
clsind = find(Yd < 3); % NC vs. MCI
[mae(1),cr(1),tb{1}] = cv_treebagger(Xd(clsind,:),Yd(clsind) - 1,foldid{1}(compind(clsind)),1);
meth = 'classification';
tb_complete{1} = TreeBagger(200,Xd(clsind,:),Yd(clsind) - 1,'OOBPrediction','on','method',meth); 
clsind = find(Yd > 1);  % MCI vs AD
[mae(2),cr(2),tb{2}] = cv_treebagger(Xd(clsind,:),Yd(clsind) - 2,foldid{1}(compind(clsind)),1);
tb_complete{2} = TreeBagger(200,Xd(clsind,:),Yd(clsind) - 2,'OOBPrediction','on','method',meth);

save tb_complete tb_complete
% for i = 1:3
% iiii{i} = find(Yd == i);
%   mmm(i) = mean(Xd(iiii{i},2));
%   sss(i) = std(Xd(iiii{i},2));
% end
% close all
% scatter(Xd(iiii{1},1),Xd(iiii{1},2),'r')
% hold
% scatter(Xd(iiii{2},1),Xd(iiii{2},2),'b')
% scatter(Xd(iiii{3},1),Xd(iiii{3},2),'g')
% mmm
% sss
% (mmm(1) - mmm(2))/sss(2)