% compute_priors
% computing the priors for all kinds of conversion intervals 
priorvec = [];
for i = 1:nsubj
    if length(time_points{i,4}) > 1
        if measurements{i,4}(1,4) == 1 % normal at the baseline
            idx = find(measurements{i,4}(:,4) == 2);
            if ~isempty(idx)
                priorvec = [priorvec; [0 time_points{i,4}(min(idx))]]; % converted at xxx months
            else
                priorvec = [priorvec;[1 time_points{i,4}(end)]]; % did not convert; data for xxx months
            end
        end
    end
end
months = [0:6:48];
conversion_probability_nl = [];
for m = 1:length(months)
    n1 = length(find(priorvec(:,1) == 0 & priorvec(:,2) < (months(m) + 1)));
    n2 = length(find(priorvec(:,1) == 0 & priorvec(:,2) > (months(m))));
    n3 = length(find(priorvec(:,1) == 1 & priorvec(:,2) > (months(m) - 1)));
    conversion_probability_nl = [conversion_probability_nl; months(m) n1/(n1 + n2 + n3)];
end
priorvec2 = priorvec;

priorvec = [];
for i = 1:nsubj
    if length(time_points{i,4}) > 1
        if measurements{i,4}(1,4) == 2 % mci at the baseline
            idx = find(measurements{i,4}(:,4) == 3);
            if ~isempty(idx)
                priorvec = [priorvec; [0 time_points{i,4}(min(idx))]]; % converted at xxx months
            else
                priorvec = [priorvec;[1 time_points{i,4}(end)]]; % did not convert; data for xxx months
            end
        end
    end
end
months = [0:6:48];
conversion_probability_ad = [];
for m = 1:length(months)
    n1 = length(find(priorvec(:,1) == 0 & priorvec(:,2) < (months(m) + 1)));
    n2 = length(find(priorvec(:,1) == 0 & priorvec(:,2) > (months(m))));
    n3 = length(find(priorvec(:,1) == 1 & priorvec(:,2) > (months(m) - 1)));
    conversion_probability_ad = [conversion_probability_ad; months(m) n1/(n1 + n2 + n3)];
end

% interpolate 
xxx = [0:1:48];
conversion_probability_ad_i = [interp1(months,conversion_probability_ad,xxx)];
conversion_probability_nl_i = [interp1(months,conversion_probability_nl,xxx)];
