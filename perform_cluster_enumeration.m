% This function estimates the number of clusters based on modularity
% maximization.
% For details, see:
%
% [1] A. Taştan, M. Muma and A. M. Zoubir, “Robust Regularized
% Locality Preserving Indexing for Fiedler Vector Estimation,”
% Signal Process. (accepted), 2023.
%
% [2] A. Clauset, M. E. J. Newman and C. Moore, “Finding community
% structure in very large networks,” Phys. Rev. E, vol. 70, pp. 066111,
% 2004.
%
% Copyright (C) 2023 Aylin Tastan. All rights reserved.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%
%Inputs:
%      y_hat             : (A numeric) estimated Fiedler vector of size nx1
%      W                 : (A numeric) affinity matrix of size n x n
%      spec_num_clusters : a vector of candidate number of clusters
%      partitioning_opt  : a string for preferred partitioning function
%                          'kmeans','kmedoids' or 'kmedoidswithTukey'
%                          (default is 'kmeans').
%
% Outputs:
%         est_num_clusters     : estimated number of clusters in the
%                                observed data set
%         partition_modularity : modularity value of partition
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [est_num_clusters,partition_modularity,c_hat] = performClusterEnumeration(y_hat,W,spec_num_clusters,partitioning_opt)

if nargin < 4 || isempty(partitioning_opt)
    partitioning_opt = 'kmeans';
end

%% Perform cluster enumeration
MAXIT = 10;

for iter = 1:MAXIT
    for i = 1:length(spec_num_clusters)

        if(isequal(partitioning_opt ,'kmeans'))

            [c_hat,~] = kmeans(y_hat,spec_num_clusters(i),'Replicates',10);

        else if(isequal(partitioning_opt ,'kmedoids'))

                [c_hat,~] = kmedoids(y_hat,spec_num_clusters(i),'Replicates',10);

        else
            [c_hat,~] = kmedoids(y_hat,spec_num_clusters(i),'Distance',@distfun_tuk, 'Start', 'sample');

        end
        end

        mod(i) = calculate_modularity(c_hat,W);
    end

    [~,indcand_K] = max(mod);
    cand_num_clusters(iter) = indcand_K;
end

est_num_clusters = mode(cand_num_clusters,'all');
if(isequal(partitioning_opt ,'kmeans'))

    [c_hat,~] = kmeans(y_hat,est_num_clusters,'Replicates',10);

else if(isequal(partitioning_opt ,'kmedoids'))

        [c_hat,~] = kmedoids(y_hat,est_num_clusters,'Replicates',10);

else

    [c_hat,~] = kmedoids(y_hat,est_num_clusters,'Distance',@distfun_tuk, 'Start', 'sample');

end
end

partition_modularity = calculate_modularity(c_hat,W);

end
