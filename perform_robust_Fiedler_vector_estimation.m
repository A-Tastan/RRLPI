% This function performs Robust Fiedler Vector Estimation based on RRLPI.
% For details, see:
%
% [1] A. Taştan, M. Muma and A. M. Zoubir, “Robust Regularized
% Locality Preserving Indexing for Fiedler Vector Estimation,”
% Signal Process. (accepted), 2023.
%
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
% Inputs:
%       X                 : (numeric) data matrix of size m x n
%                           (m: num_features, n: num_samples_total)
%       num_samples_total : a scalar contains the total number of samples
%       gamma_cand        : a vector containing the candidate penalty
%                           parameters
%       min_num_samples   : minimum number of samples in per set for delta
%                           seperated sets
%       sim_measure       : a string for preferred similarity measure,
%                           'cosine','Pearson' or 'enet'
%                           (default is 'cosine')
%       func_opt          : a scalar contains function option(default is 1)
%                           1-> overall edge weight-based regularized
%                               locality preserving indexing,
%                           2-> residual-based regularized locality
%                               preserving indexing
%       decision_rule     : a scalar contains decision rule for delta
%                           seperated sets; 1 for cutting using zero,
%                           2 for cutting using median (default is 1)
%       lambda            : (optional)for elastic net similarity penalty
%                           value (default is 0.5)
%       plotting          : a binary variable that indicates whether a
%                           figure will be generated or not (default is 0)
%
%
% Outputs:
%        y_hat            : (numeric) Fiedler vector of size n x 1
%        beta_hat         : (numeric) the regression coefficient vector of
%                           size m x 1
%        W                : (numeric) affinity matrix of size n x n
%
% Version      : March 8, 2023
% Dependencies : Normalize_Feature_Vectors,
%                Compute_Affinity_Matrix,
%                Initialize_Parameters,
%                Robust_Regularized_Locality_Preserving_Indexing,
%                Check_Delta_Seperated_Sets.
% Author       : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [y_hat,beta_hat,W] = perform_robust_Fiedler_vector_estimation(X,num_samples_total,gamma_cand,min_num_samples,sim_measure,func_opt,decision_rule,lambda,plotting)

%% Define parameters
if nargin < 9 || isempty(plotting)
    plotting = 0;
end

if nargin < 8 || isempty(lambda)
    lambda = 0.5;
end

if nargin < 7 || isempty(decision_rule)
    decision_rule = 1;
end

if nargin < 6 || isempty(func_opt)
    func_opt = 1;
end

if nargin < 5 || isempty(sim_measure)
    sim_measure = 'cosine';
end

%% Similarity Measure
X = normalize_feature_vectors(X,num_samples_total); %the feature vectors are normalized i.e. ||x||_2
W = compute_affinity_matrix(X,num_samples_total,sim_measure,lambda);

%% Initialization
[y,beta_0] = initialize_parameters(W,X); %Initial estimate of Fiedler vector y and regression coefficient vector beta_0

%% Penalty Parameter Selection
for i = 1:length(gamma_cand)
    
    %Robust Fiedler Vector Estimation
    [y_cand,~] = robust_regularized_locality_preserving_indexing(X,W,y,gamma_cand(i),func_opt,beta_0);
    
    %Delta Separated Sets
    [score_delta,gap_initial] = check_delta_separated_sets(y_cand,min_num_samples,decision_rule);
    
    %Store Initial gap and Delta scores
    Gap_val(i) = gap_initial;
    Score_delta_val(i) = score_delta;
    
end

if (any(Score_delta_val)) %If any subsets do not provide delta separation,
    [~,ind_gamma_est] = max(Score_delta_val);
else %pick the maximum initial gap.
    [~,ind_gamma_est] = max(Gap_val);
end

gamma_est = gamma_cand(ind_gamma_est);

%% Perform Fiedler Vector Estimation for Selected Penalty
[y_hat,beta_hat] = robust_regularized_locality_preserving_indexing(X,W,y,gamma_est,func_opt,beta_0);

%% Plotting
if(plotting)
    %LE
    figure(1);
    scatter(1:num_samples_total,y,'d','MarkerFaceColor',[0.7294    0.8706    0.6078],'MarkerEdgeColor',[0.6196    0.7608    0.6157])
    set(gca,'Fontname','Times','FontSize',12)
    set(gca,'FontSize',16)
    grid on
    ylabel('$\hat{\mathbf{y}}$','interpreter','latex','Fontsize',16);
    xlabel('$n$','interpreter','latex','Fontsize',16);
    title('Fiedler vector estimation based on  LE');
    
    %RRLPI
    figure(2);
    scatter(1:num_samples_total,y_hat,'p','MarkerFaceColor',[0.8392    0.4863    0.4863],'MarkerEdgeColor',[0.0784    0.3059    0.4588])
    set(gca,'Fontname','Times','FontSize',12)
    set(gca,'FontSize',16)
    grid on
    ylabel('$\hat{\mathbf{y}}$','interpreter','latex','Fontsize',16);
    xlabel('$n$','interpreter','latex','Fontsize',16);
    title('Robust Fiedler vector estimation based on  RRLPI');
end

end
