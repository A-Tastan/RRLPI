% This function performs Robust Fiedler Vector Estimation. For details:
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
%Inputs:
%      X          :(A numeric) data matrix of size m x n
%                  (m: number of features, n: number of observations)
%      W          :(A numeric) affinity matrix of size n x n
%      y          : initial Fiedler vector of size n x 1
%      gamma      : a scalar contains penalty parameter
%      func_opt   : a scalar contains function option (default is 1)
%                   1-> degree-based regularized locality preserving indexing
%                   2-> residual-based regularized locality preserving indexing
%      beta_0     :(A numeric) regression vector of size m x 1
%
% Outputs:
%      y_hat      : (A numeric) estimated Fiedler vector of size n x 1
%      beta_hat   : (A numeric) estimated regression vector of size m x 1
%
% Version      : March 8, 2023
% Dependencies : whub,
%                RankDeficientFastDecomposition (if XX is rank deficient).
%                An example MATLAB function is provided in :
%                [2]P. Courrieu, “Fast computation of Moore-Penrose inverse
%                matrices.” Online-Edition: https://arxiv.org/abs/0804.4809,
%                2008.
% Author       : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [y_hat,beta_hat] = robust_regularized_locality_preserving_indexing(X,W,y,gamma,func_opt,beta_0)

%% Define Parameters
ITERMAX = 1000;
TOL = 1.0e-5;

if isreal(y)
    const = 1.4815;
    c = 1.345;    % 95 percent ARE
    %c = 0.7317;  % 85 percent ARE
else
    const = 1.20112;
    c = 1.214;     % 95 percent ARE
    %c = 0.515;   % 85 percent ARE
end

wfun = @(x) whub(x,c);

%% Overall edge weight-based Regularized Locality Preserving Indexing
% Solve regularized least squares problem
if(func_opt == 1)
    
    %Auxiliary scale estimate based on degree
    D = sum(W); %diaognal weight matrix
    resid = abs(D - median(D)).';
    sig = const * median(resid(resid~=0)); % auxiliary scale estimate
    
    %Regularized least squares
    w = wfun(resid/sig);
    Xstar = bsxfun(@times, X.', w);
    XX = full(Xstar'*X.');
    if(gamma)
        for i=1:size(XX,1)
            XX(i,i) = XX(i,i) + gamma*sig^2;
        end
    end
    
    XX = max(XX,XX');
    B = Xstar'*y;
    if rank(XX) < size(XX,2)
        R = RankDeficientFastDecomposition(XX);
    else
        R = chol(XX);
    end
    beta_hat = R \ (R' \ B);
    
    
    %% Residual-based Regularized Locality Preserving Indexing
    % Solve iteratively reweighted regularized least squares
else
    
    %Auxiliary scale estimate
    resid = abs( y - X.' * beta_0);
    sig = const * median(resid(resid~=0)); % auxiliary scale estimate
    
    %Iteratively reweighted regularized least squares
    for iter = 1:ITERMAX
        
        resid(resid <.000001) =.000001;
        w = wfun(resid / sig);
        Xstar = bsxfun(@times, X.', w);
        XX = full(Xstar'*X.');
        if(gamma)
            for i=1:size(XX,1)
                XX(i,i) = XX(i,i) + gamma*sig^2;
            end
        end
        XX = max(XX,XX');
        B = Xstar' * y;
        if rank(XX) < size(XX,2)
            R = RankDeficientFastDecomposition(XX);
        else
            R = chol(XX);
        end
        
        R = chol(XX); %XX=R'R
        beta_hat = R \ (R' \ B);
        crit = norm(beta_hat-beta_0)/norm(beta_0);
        
        if crit < TOL
            break;
        end
        
        beta_0 = beta_hat;
        resid = abs( y - X.' * beta_0);
        
    end
    
end

y_hat = X.' * beta_hat;  %y = X.'beta

end