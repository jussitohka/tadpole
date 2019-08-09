indOI = [ADAS13idx ventricleICVidx];
fmeas = zeros(nsubj,2,2);
fR = zeros(nsubj,2);
fError = zeros(nsubj,2,3);
flag = 0;
tpthr = [7 4];
con = [0 0.5];
for j = 1:2
   for i = 1:nsubj
        if length(time_points{i,j}) > 1
            X = measurements{i,j}((end - 1),:);
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
            fError(i,j,3) = abs(fmeas(i,2,j) -  measurements{i,j}(end,indOI(j)));
        end
        if length(time_points{i,j}) > 2 
            [fmeas(i,1,j),fR(i,j)] = longitudinal_estimate(measurements{i,j}(1:(end - 1),indOI(j)),time_points{i,j}(1:(end - 1)),time_points{i,j}(end));
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