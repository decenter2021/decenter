function LinDynamics = getDiscreteDynamicsWL(x,Cte)
%% Description
% This function computes the linearized discrete model for a given state.
% The equilibrium state is computed around the level in tanks 1 and 2.
% Input:    - x: state vector
%           - dT: sampling period
%           - Cte: struct with constants of the model dynamics 
% Output:   - LinDynamics: 1x7 cell with matrices A,C,Q,R,B, equilibrium
% state vector, xEq, and equilibrium control vector, uEq.

%% Compute linearized model

% Initialize 1x7 cell to store the model matrices and equilibrium vectors
LinDynamics = cell(1,6);
% Supress output of the numeric solver
%options = optimset('Display','off');

% Initialize equilibrium state vector
xEq = zeros(4,1);
% Level in the first and second tanks define the equilibrium
xEq(1:2,1) = x(1:2,1);

% In equilibrium solve for level in the third and fouth tanks and the 
% control vector, defined by the level on the first and second tanks
%Xsol = real(fsolve(@(X) xdot([x(1:2);X(1:2)],X(3:4),Cte),[x(3:4);0;0], options));

xEq(3:4) = Cte.alpha*[x(1,1);x(2,1); sqrt(x(1,1)*x(2,1))];
uEq = Cte.beta*[sqrt(x(1,1)); sqrt(x(2,1))];

% Update equilibrium vectors
%xEq(3:4) = Xsol(1:2);
%uEq = Xsol(3:4);

% Compute linearized entries of matrices A,B,C,Q, and R for a coninuous
% process
[A,B] = Dxdot(xEq,Cte);

G = expm([A B; zeros(Cte.m,Cte.n) zeros(Cte.m,Cte.m)]*Cte.dT);
A = G(1:Cte.n,1:Cte.n);
B = G(1:Cte.n,Cte.n+1:end);
LinDynamics{1,1}(1:Cte.n,1:Cte.n) = A;
LinDynamics{1,1}(Cte.n+1:3*Cte.n/2,1:Cte.n) = Cte.h*A;
LinDynamics{1,1}(1:Cte.n,Cte.n+1:3*Cte.n/2) = zeros(Cte.n,Cte.n/2);
LinDynamics{1,1}(Cte.n+1:3*Cte.n/2,Cte.n+1:3*Cte.n/2) = eye(Cte.n/2);
LinDynamics{1,2}(1:Cte.n,1:Cte.m) = B;
LinDynamics{1,2}(Cte.n+1:3*Cte.n/2,1:Cte.m) = Cte.h*B;
LinDynamics{1,4} = Cte.R;
% Equilibrium states
LinDynamics{1,5} = [xEq;zeros(Cte.n/2,1)];
LinDynamics{1,6} = uEq;

%%%%%%%%%%% ANtes
LinDynamics{1,3} = Cte.H'*Cte.Qt*Cte.H;

%%%% Agora
%LinDynamics{1,3} = Cte.QtNovo;

end

