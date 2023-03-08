function [X, num_features, num_samples_total] = call_data_set()

load fisheriris

X = meas.';
[num_features,num_samples_total] = size(X);

end
