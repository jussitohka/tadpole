function [ mae,cr,tb,Yhat ] = cv_treebagger(X,Y,foldid,cls )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('cls','var')
    cls = 0;
end
if cls
    meth = 'classification';
else
    meth = 'regression';
end

Yhat = zeros(size(Y));
for fold = 1:max(foldid)
    trainsubj = ~(foldid == fold);
    testsubj = (foldid == fold);
    [tb{fold}]= TreeBagger(200,X(trainsubj,:),Y(trainsubj,1),'OOBPrediction','on','method',meth);
    if ~cls
        Yhat(testsubj) = predict(tb{fold},X(testsubj,:));
    else
        [yyy scores]  = predict(tb{fold},X(testsubj,:));
        Yhat(testsubj) = scores(:,2);
    end
end
if ~cls
    mae = mean(abs(Y - Yhat));
    cr = corr(Y,Yhat);
else
    mae = mean(abs(Y - Yhat));
    [~,~,~,cr] = perfcurve(Y,Yhat,1);
end



end

