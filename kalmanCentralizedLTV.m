function [K,Ppred,Pfilt] = kalmanCentralizedLTV(system,Pprev)
%% Description
% This function computes the centralized kalman filter gain for a time-instant
% k, i.e., K(k), subject to a sparsity constraint.
% Input:    - system: 1x4 cell whose rows contain matrices A,C,Q and R i.e.,
%               - system{1,1} = A(k)
%               - system{2,2} = C(k)
%               - system{3,3} = Q(k)
%               - system{4,4} = R(k)
%           - Pprev: Previous predicted error covariance matrix, i.e., P(k|k-1)
% Output:   - K: filter gain K(k)
%           - Ppred: Predicted error covarinace matrix P(k+1|k)
%           - Pfilt: Predicted error covarinace matrix P(k|k)

%% Gain computation 
n = size(system{1,1},1); % Get value of n from the size of A 
% Compute gain
K = Pprev*transpose(system{1,2})/(system{1,2}*Pprev*transpose(system{1,2})+system{1,4});
% Update the covariance matrix after the filtering step, i.e., P(k|k)
Pfilt = K*system{1,4}*K'+...
    (eye(n)-K*system{1,2})*Pprev*((eye(n)-K*system{1,2})');
% Compute P(k+1|k)
Ppred = system{1,1}*Pfilt*system{1,1}'+system{1,3};
end