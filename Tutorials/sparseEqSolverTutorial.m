%% Tutorial of sparseEqSolver 
% Algorithm developed in [1]
%% Generate random matrices A,B,C, and sparsity parttern E 
n = 5;
o = 4;
A = rand(n,n)
B = rand(o,o)
C = rand(n,o)
E = round(rand(n,o))
%% Solve 
% [AXB-C]_ij = 0    if E_ij ~= 0 
% [X]_ij = 0        if E_ij = 0 
X = sparseEqSolver(A,B,C,E)
%% Verify solution
sum(sum(abs((A*X*B-C).*E)))

%% References
% [1] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497