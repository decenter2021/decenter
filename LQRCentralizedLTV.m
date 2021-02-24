function [K,P] = LQRCentralizedLTV(system,T,opts)
%% Description
% This function computes the one-step LQR regulator gain for a window T
% Input:    - system: (T+1)x4 cell whose rows contain matrices A,B,Q and R 
%           for the whole window, i.e.,
%               - system{i,1} = A(k+i-1), i = 1,...,T
%               - system{i,2} = B(k+i-1), i = 1,...,T
%               - system{i,3} = Q(k+i-1), i = 1,...,T+1
%               - system{i,4} = R(k+i-1), i = 1,...,T
%           Note that for the entry (T+1) only Q is used.
%           - T: window length (T gains computed)
% Output:   - K: Tx1 cell of gains for the whole window,i.e.,
%                 K{i,1}, tau = k+i-1,...,k+i-2+T
%           - P: (T+1)x1 cell of P matrices for the whole window, i.e.,
%                 P{i,1}, tau = k+i-1,...,k+T
% Important notes: 
%           - output gain corresponds to the control law: u(k)=-K(k)*x(k)

%% Argument handling
if ~exist('opts','var') 
    opts.verbose = false; % Default is not to display algorithm status messages
elseif ~isfield(opts,'verbose')
    opts.verbose = false; % Default is not to display algorithm status messages
end
if opts.verbose
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Running centralized algorithm with T = %d.\n',T);
    fprintf('----------------------------------------------------------------------------------\n');
end
%% Gain computation
persistent n
if isempty(n)
    n = size(system{1,1},1); % Get value of n from the size of A 
end
persistent m
if isempty(m)
    m = size(system{1,2},2); % Get value of n from the size of B 
end
P = cell(T+1,1); % Initialize cell arrays
K = cell(T,1);
P{T+1,1} = system{T+1,3}; % terminal condition  
for k = T:-1:1
    % Compute P
    S = system{k,4}+system{k,2}'*P{k+1,1}*system{k,2};
    K{k,1} = S\system{k,2}'*P{k+1,1}*system{k,1};
    % Update P
    P{k,1} = system{k,3}+K{k,1}'*system{k,4}*K{k,1}+...
       (system{k,1}-system{k,2}*K{k,1})'*P{k+1,1}*(system{k,1}-system{k,2}*K{k,1});       
end
end