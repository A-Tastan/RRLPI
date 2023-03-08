# RRLPI
Robust Regularized Locality Preserving Indexing for Fiedler Vector Estimation

The Fiedler vector is the eigenvector associated with the algebraic connectivity of the graph Laplacian. It is central to graph analysis as it provides substantial information to learn the latent structure of a graph. In real-world applications, however, the data may be subject to heavy-tailed noise and outliers which deteriorate the structure of the Fiedler vector estimate and lead to a breakdown of popular methods. Thus, we propose a Robust Regularized Locality Preserving Indexing (RRLPI) Fiedler vector estimation method that approximates the nonlinear manifold structure of the Laplace Beltrami operator while minimizing the impact of outliers. To achieve this aim, an analysis of the effects of two fundamental outlier types on the eigen-decomposition of block affinity matrices is conducted. Then, an error model is formulated based on which the RRLPI method is developed. It includes an unsupervised regularization parameter selection algorithm that leverages the geometric structure of the projection space. The performance is benchmarked against existing methods in terms of detection probability, partitioning quality, image segmentation capability, robustness and computation time using a large variety of synthetic and real data experiments.

For details, see:

[1] A. Taştan, M. Muma and A. M. Zoubir, “Robust Regularized Locality Preserving Indexing for Fiedler Vector Estimation,” Signal Process. (accepted), 2023.

The codes can be freely used for non-commercial use only. Please make appropriate references to our article.

NOTE : This code requires an additional function 'whub' which is available in :

https://github.com/RobustSP/toolbox/blob/master/codes/07_AuxiliaryFunctions/whub.m

Further, it may require "enet.m" and "SoftThresh.m" functions if the preferred similarity measure is the elastic net similarity. The codes are, respectively, available in :

https://github.com/RobustSP/toolbox/blob/master/codes/02_Regression/enet.m

and

https://github.com/RobustSP/toolbox/blob/master/codes/07_AuxiliaryFunctions/SoftThresh.m

Lastly, it may require an additional function if the matrix XX in "Robust_Regularized_Locality_Preserving_Indexing.m" is rank deficient. An example function that is named "RankDeficientFastDecomposition.m" code is available in:
 
[2] P. Courrieu, “Fast computation of Moore-Penrose inverse matrices.” Online-Edition: https://arxiv.org/abs/0804.4809, 2008


