function xdot = xdotContinuous(x,u,Cte)
%% Description 
% This function computes the derivative of the state vector using the non 
% linear dynamics of the model, for a given state and control vector.
% Input:    - x: state vector 
%           - u: actuation vector
% Output:   - xdot: derivative of the state vector

%% xdot computation

% Add non negative level and non negative actuation constraint
x(x<0) = 0;
u(u<0) = 0;

% Initialize xdot
xdot = zeros(4,1);

% Compute xdot acconding to the non linear model
xdot(1) = -(Cte.a1/Cte.A1)*sqrt(2*Cte.g*x(1))+...
    (Cte.a3/Cte.A1)*sqrt(2*Cte.g*x(3))+Cte.g1*Cte.k1*u(1)/Cte.A1;
xdot(2) = -(Cte.a2/Cte.A2)*sqrt(2*Cte.g*x(2))+...
    (Cte.a4/Cte.A2)*sqrt(2*Cte.g*x(4))+Cte.g2*Cte.k2*u(2)/Cte.A2;
xdot(3) = -(Cte.a3/Cte.A3)*sqrt(2*Cte.g*x(3))+(1-Cte.g2)*Cte.k2*u(2)/Cte.A3;
xdot(4) = -(Cte.a4/Cte.A4)*sqrt(2*Cte.g*x(4))+(1-Cte.g1)*Cte.k1*u(1)/Cte.A4;   

end