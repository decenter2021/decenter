function [K,P] = LQROneStepLTV(system,E,T,opts)
%% Description
% This function computes the one-step LQR regulator gain for a window T
% according to [1]
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
n = size(system{1,1},1); % Get value of n from the size of A 
P = cell(T+1,1); % Initialize cell arrays
K = cell(T,1);
P{T+1,1} = system{T+1,3}; % terminal condition  
for k = T:-1:1
   % Compute gain with efficient algorithm [2]
   K{k,1} = sparseEqSolver(system{k,4}+system{k,2}'*P{k+1,1}*system{k,2},eye(n),system{k,2}'*P{k+1,1}*system{k,1},E);
   % Update P
   P{k,1} = system{k,3}+K{k,1}'*system{k,4}*K{k,1}+...
       (system{k,1}-system{k,2}*K{k,1})'*P{k+1,1}*(system{k,1}-system{k,2}*K{k,1});
end       
end

%% References
% [1] Pedroso L, Batista P. Discrete?time decentralized linear quadratic 
% control for linear time?varying systems. International Journal of Robust 
% and Nonlinear Control. 2021 Sep 8. https://doi.org/10.1002/rnc.5772

% [2] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497