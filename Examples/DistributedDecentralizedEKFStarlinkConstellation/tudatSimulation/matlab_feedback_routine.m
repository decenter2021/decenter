%% Package: tudat-matlab-thrust-feedback
% Author: Leonardo Pedroso
%% Feedback loop implementation
% Implement the feedback-loop in this script 
% Input:   t  
%          x_t dim: 7N x 1
% Output:  u_t dim: 3 x N
% The controller parameters should be defined in matlab-app.m

%% Thrust
u_t = zeros(3*N,1);