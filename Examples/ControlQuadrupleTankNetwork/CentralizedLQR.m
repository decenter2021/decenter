function [P,K] = CentralizedLQR(system,w)
%% Description
% This function computes the centralized LQR regulator gain for a window w
% Input:    - system: (w+1)x4 cell whose rows contain matrices A,B,Q and R 
%           for the whole window. Note that for the instant (w+1) only Q is
%           used. For instant 1 
%           - E: sparsity pattern
%           - w: window length (w gains computed)
% Output:   - K: wx1 cell of gains for the whole window
%           - P: (w+1)x1 cell of P matrices for the whole window 
% Important notes: 
%           - output gain corresponds to the control law: u(k)=-K(k)^T*x(k)

%% Gain computation
    n = size(system{1,1},1); % Get value of n from the size of A   
    P = cell(w+1,1);  % Initialize cell arrays
    K = cell(w,1);
 
    P{w+1,1} = system{w+1,3}; % terminal condition 
    % Generate gain and update covariance for centralized LQR
    for k = w:-1:1
       P_ = (system{k,1})'*P{k+1,1}*system{k,1};
       K{k,1} = (P_/system{k,1})*system{k,2}/(system{k,4}+(system{k,2})'*(P{k+1,1})*system{k,2});
       P{k,1} = system{k,3}+K{k,1}*system{k,4}*(K{k,1})'+...
           (eye(n)-K{k,1}*(system{k,2})'/((system{k,1}))')*...
           P_*(eye(n)-K{k,1}*(system{k,2})'/((system{k,1}))')';
    end

end