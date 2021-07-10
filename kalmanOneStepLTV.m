function [K,Ppred,Pfilt] = kalmanOneStepLTV(system,E,Pprev)
%% Description
% This function computes the one-step kalman filter gain for a time-instant
% k, i.e., K(k), subject to a sparsity constraint.
% Input:    - system: 1x4 cell whose rows contain matrices A,C,Q and R i.e.,
%               - system{1,1} = A(k)
%               - system{2,2} = C(k)
%               - system{3,3} = Q(k)
%               - system{4,4} = R(k)
%           - E: sparsity pattern
%           - Pprev: Previous predicted error covariance matrix, i.e., P(k|k-1)
% Output:   - K: filter gain K(k)
%           - Ppred: Predicted error covarinace matrix P(k+1|k)
%           - Pfilt: Predicted error covarinace matrix P(k|k)

%% Gain computation 
n = size(system{1,1},1); % Get value of n from the size of A 
% Compute gain with sparse matrix solver [1]
K = sparseEqSolver(eye(n),system{1,2}*Pprev*system{1,2}'+system{1,4},Pprev*system{1,2}',E);
% Update the covariance matrix after the filtering step, i.e., P(k|k)
Pfilt = K*system{1,4}*K'+...
    (eye(n)-K*system{1,2})*Pprev*((eye(n)-K*system{1,2})');
% Compute P(k+1|k)
Ppred = system{1,1}*Pfilt*system{1,1}'+system{1,3};
end

%[1] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497
