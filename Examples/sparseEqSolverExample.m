rng(1);
A = rand(5,5);
B = rand(4,4);
E = round(rand(5,4));
C = rand(5,4);
X = sparseEqSolver(A,B,C,E)
sum(sum(abs((A*X*B-C).*E)))