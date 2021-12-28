function E = sparseEqSolver(A,B,C,E)
%% Description
% This function solves the sparse matrix equation 
% [AXB-C]_ij = 0    if E_ij ~= 0 
% [X]_ij = 0        if E_ij = 0 
% where A (nxn), B(oxo), C(nxo), and E(nxo) are known.
% It is assumed that there is one solution to this equation.
% Algorithm described in [1].
% Input:    - A,B,C,E
%           - E: a matrix that defines the sparsity pattern      
% Output:   - X: the solution
n = size(E,1); % Get value of n from the size of E 
o = size(E,2); % Get value of o from the size of E 
ksi = n*o-sum(sum(E==0)); % Get number of nonzero elements of E
R = zeros(ksi,1); % Allocate vector r
S = zeros(ksi,ksi); % Allocate vector S
p = 0; % Init entry counter of nonzero entries of vec(X)
for i = 1:n     % For (i,j) in chi
    for j = 1:o
        if E(i,j)==0 continue; end
        p = p + 1;
        R(p) = C(i,j);
        q = 0; % Init entry counter of nonzero entries of vec(X)
        for k = 1:n     % For (k,l) in chi
           for l = 1:o
               if E(k,l)==0 continue; end
               q = q+1;
               S(p,q) = A(i,k)*B(l,j);
           end
        end
    end 
end
vecX = S\R; % Solve the resulting system of linear equations
p = 0;
for i = 1:n     % For (i,j) in chi
    for j = 1:o
        if E(i,j)==0 continue; end
        p = p + 1;
        E(i,j) = vecX(p);
    end
end
end

%% References
%[1] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497

