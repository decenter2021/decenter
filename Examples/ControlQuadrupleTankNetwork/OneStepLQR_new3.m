function [P,K] = OneStepLQR(system,E,w)
%% Description
% This function computes the One-step LQR regulator gain for a window w
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
    m = size(system{1,2},2); % Get value of n from the size of B 
    
    P = cell(w+1,1); % Initialize cell arrays
    K = cell(w,1);

    P{w+1,1} = system{w+1,3}; % terminal condition  
    for k = w:-1:1
       K{k,1} = zeros(n,m);
%        S = system{k,4}+transpose(system{k,2})*P{k+1,1}*system{k,2};
%        P_ = system{k,1}'*P{k+1,1}*system{k,1};
       for i = 1:n
            L = zeros(n);
            L (i,i) = 1; % Generate matrix L_i
            M = zeros(m);
            for j = 1:m % Gererate matrix M_i
                if E(i,j) ~= 0
                    M(j,j) = 1;
                end
            end
            % Compute the ith term of the summation 
            K{k,1} = K{k,1} + (L*((system{k,1}'*P{k+1,1}*system{k,1})/system{k,1})*system{k,2}*M)/...
                (eye(m)-M+M*(system{k,4}+transpose(system{k,2})*P{k+1,1}*system{k,2})*M);
       end
       % Update P
       P{k,1} = system{k,3}+K{k,1}*system{k,4}*transpose(K{k,1})+...
           (eye(n)-K{k,1}*system{k,2}'/(system{k,1}'))*...
           (system{k,1}'*P{k+1,1}*system{k,1})*transpose(eye(n)-K{k,1}*system{k,2}'/(system{k,1}'));
    end       
end
