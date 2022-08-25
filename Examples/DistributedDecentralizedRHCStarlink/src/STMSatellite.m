%% Computation of the system transition matrix
function [Ad,Bd] = STMSatellite(OENominal_t0,Ts)
    
    %% Constants 
    % Orbital dynamics
    mu = 3.986004418e14; %(m^3 s^-2)
    RE = 6378.137e3; %(m)
    J2 = 1082.6267e-6;
    gamma = (3/4)*J2*RE^2;
    
    %% Covert to notation in Mauro 2018
    % Chief elements
    a_c = OENominal_t0(1);
    u_c_t0 = OENominal_t0(2);
    ex_c = OENominal_t0(3);
    ey_c = OENominal_t0(4);
    i_c = OENominal_t0(5);
    %Omega_c = OENominal_t0(6);
    
    %% Compute parameters
    % Chief 
    n_c = sqrt(mu/(a_c^3));
    e_c = sqrt(ex_c^2+ey_c^2);
    eta_c = sqrt(1-e_c^2);
    K_c = gamma*n_c/(a_c^2*eta_c^4);
    Q_c = 5*cos(i_c)^2-1;
    P_c = 3*cos(i_c)^2-1;
    F_c = 4+3*eta_c;
    E_c = 1+eta_c;
    S_c = sin(2*i_c);
    T_c = sin(i_c)^2;
    Lambda_c = (3/2)*n_c + (7/2)*E_c*K_c*P_c;
    
    %% Computation of STM (Eq (14), Mauro et al., 2018)
    Dw = K_c*Q_c*Ts;
    Ad = eye(6);
    Ad(2,1) = -Lambda_c*Ts;
    Ad(2,5) = -K_c*F_c*S_c*Ts;
    Ad(3,3) = cos(Dw);
    Ad(3,4) = -sin(Dw);
    Ad(4,3) = sin(Dw);
    Ad(4,4) = cos(Dw);
    Ad(6,1) = (7/2)*K_c*S_c*Ts;
    Ad(6,5) = 2*K_c*T_c*Ts;
    
    %% Compute more parameters
    W_c = n_c + K_c*Q_c +eta_c*K_c*P_c;
    C = K_c*Q_c/W_c;
    u_c_t = u_c_t0 + Ts*W_c;
    delta_u = Ts*W_c;
    
    %% Compute convolution matrix (Eq (16), Mauro et al., 2018)
    % In TNW local frame (X: velocity; Y: - along radius; Z: angular momentum)
    
    Bd = zeros(6,3);
    Bd(1,1) = 2*delta_u/(n_c*a_c*W_c);
    
    Bd(2,2) = (-1)*-2*delta_u/(n_c*a_c*W_c);
    %Bd(2,2) = -(-1)*-2*delta_u/(n_c*a_c*W_c);
    Bd(2,1) = -Lambda_c*delta_u^2/(n_c*a_c*W_c^2);
    Bd(2,3) = F_c*K_c*S_c*(cos(u_c_t)-cos(u_c_t0)+sin(u_c_t0)*delta_u)/...
        (n_c*a_c*W_c^2);
    
    Bd(3,2) = (-1)*-(cos(u_c_t)-cos(u_c_t0+C*delta_u))/(n_c*a_c*(1-C)*W_c);
    %Bd(3,2) = -(-1)*-(cos(u_c_t)-cos(u_c_t0+C*delta_u))/(n_c*a_c*(1-C)*W_c);
    Bd(3,1) = 2*(sin(u_c_t)-sin(u_c_t0+C*delta_u))/(n_c*a_c*(1-C)*W_c);
    
    Bd(4,2) = (-1)*-(sin(u_c_t)-sin(u_c_t0+C*delta_u))/(n_c*a_c*(1-C)*W_c);
    %Bd(4,2) = -(-1)*-(sin(u_c_t)-sin(u_c_t0+C*delta_u))/(n_c*a_c*(1-C)*W_c);
    Bd(4,1) = -2*(cos(u_c_t)-cos(u_c_t0+C*delta_u))/(n_c*a_c*(1-C)*W_c);
    
    Bd(5,3) = (sin(u_c_t)-sin(u_c_t0))/(n_c*a_c*W_c);
    %Bd(5,3) = -(sin(u_c_t)-sin(u_c_t0))/(n_c*a_c*W_c);
    
    Bd(6,1) = (7/2)*K_c*S_c*delta_u^2/(n_c*a_c*W_c^2);

    Bd(6,3) = -(W_c+2*K_c*T_c)*(cos(u_c_t)-cos(u_c_t0))/(n_c*a_c*W_c^2)...
        - 2*K_c*T_c*sin(u_c_t0)*delta_u/(n_c*a_c*W_c^2);
    %Bd(6,3) = +(W_c+2*K_c*T_c)*(cos(u_c_t)-cos(u_c_t0))/(n_c*a_c*W_c^2)...
    %    + 2*K_c*T_c*sin(u_c_t0)*delta_u/(n_c*a_c*W_c^2);
    
end
