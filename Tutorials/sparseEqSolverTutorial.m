%% Example of sparseEqSolver 
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