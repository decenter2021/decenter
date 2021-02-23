function [A,B] = Dxdot(x_eq,Cte)
%% Description 
% This function computes matrices A and B of the continuous linearized 
% model, given an equilibrium state.
% Input:    - x_eq: equilibrium state vector
%           - Cte: struct with constants of the model dynamics

%% Compute A and B

% Initialize matrices
A = zeros(4,4);
B = zeros(4,2);

% Evaluate the non null partial derivatives in equilibrium
if x_eq(1) > 0
    A(1,1) = -((Cte.a1/Cte.A1)*Cte.g)/sqrt(2*Cte.g*x_eq(1));
else
    A(1,1) = 0;
end
if x_eq(2) > 0
    A(2,2) = -((Cte.a2/Cte.A2)*Cte.g)/sqrt(2*Cte.g*x_eq(2));
else
    A(2,2) = 0;
end

if x_eq(3) > 0
    A(1,3) = ((Cte.a3/Cte.A1)*Cte.g)/sqrt(2*Cte.g*x_eq(3));
    A(3,3) = -((Cte.a3/Cte.A3)*Cte.g)/sqrt(2*Cte.g*x_eq(3));
else
    A(1,3) = 0;
    A(3,3) = 0;
end

if x_eq(4) > 0
    A(2,4) = ((Cte.a2/Cte.A2)*Cte.g)/sqrt(2*Cte.g*x_eq(4));
    A(4,4) = -((Cte.a4/Cte.A4)*Cte.g)/sqrt(2*Cte.g*x_eq(4));
else
    A(2,4) = 0;
    A(4,4) = 0;
end
B(1,1) = Cte.g1*Cte.k1/Cte.A1;
B(2,2) = Cte.g2*Cte.k2/Cte.A2;
B(3,2) = (1 - Cte.g2)*Cte.k2/Cte.A3;
B(4,1) = (1-Cte.g1)*Cte.k1/Cte.A4;

end

