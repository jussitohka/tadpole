function [ fmeas, fR] = longitudinal_estimate(meas,tp,ftp)
tp = [ones(length(tp),1) tp];
meas(isnan(meas)) = nanmean(meas);
b = tp\meas;
fmeas = b(1) + b(2)*ftp;
measest = tp*b;
fR = 1 - sum((meas - measest).^2)/sum((meas - mean(meas)).^2);


end

