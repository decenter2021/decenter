function [K,P] = LQROneStepLTV(system,E,T,opts)
%% Description
% This function computes the one-step LQR regulator gain for a window T
% Input:    - system: (T+1)x4 cell whose rows contain matrices A,B,Q and R 
%           for the whole window, i.e.,
%               - system{i,1} = A(k+i-1), i = 1,...,T
%               - system{i,2} = B(k+i-1), i = 1,...,T
%               - system{i,3} = Q(k+i-1), i = 1,...,T+1
%               - system{i,4} = R(k+i-1), i = 1,...,T
%           Note that for the entry (T+1) only Q is used.
%           - E: sparsity pattern
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
    fprintf('Running one-step algorithm with T = %d.\n',T);
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
persistent M L
if isempty(M) || isempty(L)
   M = cell(m,1);
   L = cell(m,1);
   for j = 1:n
        L{j,1} = zeros(n);
        L{j,1}(j,j) = 1; % Generate matrix L_i
        M{j,1} = zeros(m);
        for i = 1:m % Gererate matrix M_i
            if E(i,j) ~= 0
                M{j,1}(i,i) = 1;
            end
        end
   end
end
P = cell(T+1,1); % Initialize cell arrays
K = cell(T,1);
P{T+1,1} = system{T+1,3}; % terminal condition  
for k = T:-1:1
   K{k,1} = zeros(m,n);
   S = system{k,4}+system{k,2}'*P{k+1,1}*system{k,2};
   for j = 1:n
        % Compute the ith term of the summation   
        K{k,1} = K{k,1} + ...
         (eye(m)-M{j,1}+M{j,1}*S*M{j,1})\M{j,1}*system{k,2}'*P{k+1,1}*system{k,1}*L{j,1};
   end
   % Update P
   P{k,1} = system{k,3}+K{k,1}'*system{k,4}*K{k,1}+...
       (system{k,1}-system{k,2}*K{k,1})'*P{k+1,1}*(system{k,1}-system{k,2}*K{k,1});
end       
end