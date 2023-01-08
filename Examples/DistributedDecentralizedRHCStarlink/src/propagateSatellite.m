function [x_pred] = propagateSatellite(x,u,Ts)
    %% State prediction with J2 perturbation by soving ODE
    t = [0 Ts]; % Set time instants to evaluate x
    options = odeset('RelTol',1e-10); % ODE solver options
    [~,x] = ode45(@(t,x) f(t,x,u),t,x,options); % Solve ODE
    x_pred = x(end,:)'; % Output solution at next time instant   
end

function dxdt = f(t,x,u)
    %% Constants 
    mu = 3.986004418e14; %(m^3 s^-2)
    RE = 6371e3; %(m)
    J2 = 1082.6267e-6;
    wE = 7.2921150e-5; %(rad/s)
    % Atmosphere cf. Table 7-4 Vallado1997
    h0 = 500e3; %(m)
    rho0 = 6.967e-13; %(Kg/m^3)
    H = 63.822e3; %(m)
    % Satellite
    Cd = 2.2;
    Ad = 24; % (m^2)
    % Thrust
    Ct1 = 0.068; % (N)
    Isp = 1640; % (s)
    g0 = 9.81; % (ms^-2)
    fs = 1;
    %% State vector 
    rx = x(1);
    ry = x(2);
    rz = x(3);
    vx = x(4);
    vy = x(5);
    vz = x(6);
    m = x(7);

    % TNW frame 
    % (X: velocity; Y: - along radius; Z: angular momentum)
    ux = x(4:6)/norm(x(4:6));
    uz = cross(x(1:3),x(4:6));
    uz = uz/norm(uz);
    uy = cross(uz,ux);
    
    % Saturate thrust
    u(u>Ct1) = Ct1;
    u(u<-Ct1) = -Ct1;
    
    % Thrust in cartesian frame
    f = u(1)*ux + u(2)*uy + u(3)*uz;
    fdx = f(1);
    fdy = f(2);
    fdz = f(3);
    
    %% Compute spacecraft dynamics (cf. Mathematica)
    dxdt = zeros(6,1);
    % Kinematics
    dxdt(1:3) = x(4:6);
    % Dynamics
    % Central body acceleration
    accKep = zeros(3,1);
    accKep(1) = (-1)*mu*rx*(abs(rx)^2+abs(ry)^2+abs(rz)^2)^(-3/2);
    accKep(2) = (-1)*mu*ry*(abs(rx)^2+abs(ry)^2+abs(rz)^2)^(-3/2);
    accKep(3) = (-1)*mu*rz*(abs(rx)^2+abs(ry)^2+abs(rz)^2)^(-3/2);
    % J2 Acceleration
    accJ2 = zeros(3,1);
    accJ2(1) = (15/2)*J2*mu*RE^2*rx*rz^2*(rx^2+ry^2+rz^2)^(-7/2)+(-3/2)*J2*mu* ...
                RE^2*rx*(rx^2+ry^2+rz^2)^(-5/2);
    accJ2(2) = (15/2)*J2*mu*RE^2*ry*rz^2*(rx^2+ry^2+rz^2)^(-7/2)+(-3/2)*J2*mu* ...
                RE^2*ry*(rx^2+ry^2+rz^2)^(-5/2);
    accJ2(3) = (15/2)*J2*mu*RE^2*rz^3*(rx^2+ry^2+rz^2)^(-7/2)+(-9/2)*J2*mu*RE^2* ...
                rz*(rx^2+ry^2+rz^2)^(-5/2);
    % Atmospheric drag acceleration
    accD = zeros(3,1);
    accD(1) = (-1/2)*Ad*Cd*exp(1)^(H^(-1)*(h0+RE+(-1)*(abs(rx)^2+abs(ry)^2+abs( ...
                rz)^2)^(1/2)))*m^(-1)*rho0*(vx+ry*wE)*(abs(vz)^2+abs(vy+(-1)*rx* ...
                wE)^2+abs(vx+ry*wE)^2)^(1/2);
    accD(2) = (-1/2)*Ad*Cd*exp(1)^(H^(-1)*(h0+RE+(-1)*(abs(rx)^2+abs(ry)^2+abs( ...
                rz)^2)^(1/2)))*m^(-1)*rho0*(vy+(-1)*rx*wE)*(abs(vz)^2+abs(vy+(-1)* ...
                rx*wE)^2+abs(vx+ry*wE)^2)^(1/2);
    accD(3) = (-1/2)*Ad*Cd*exp(1)^(H^(-1)*(h0+RE+(-1)*(abs(rx)^2+abs(ry)^2+abs( ...
                rz)^2)^(1/2)))*m^(-1)*rho0*vz*(abs(vz)^2+abs(vy+(-1)*rx*wE)^2+abs( ...
                vx+ry*wE)^2)^(1/2);
    %Thruster 
    accT = zeros(3,1);
    accT(1) = fdx*fs/m;
    accT(2) = fdy*fs/m;
    accT(3) = fdz*fs/m;
    % Complete dynamics
    dxdt(4:6) = accD + accJ2 +  accT + accKep;
    % Mass dynamics
    dxdt(7) = -(sum(abs(u)))/(Isp*g0);
end
