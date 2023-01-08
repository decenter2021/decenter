function OE = nominalConstellationOEMean(t,i,walkerParameters,semiMajorAxis,walkerAnchor,parallel)
    % Constants
    mu = 3.986004418e14; %(m^3 s^-2)
    RE = 6378.137e3; %(m)
    J2 = 1082.6267e-6;
    % Compute parameters
    n = sqrt(mu/(semiMajorAxis)^3);
    gamma = (J2/2)*(RE/semiMajorAxis)^2;
    Omega_dot = -3*gamma*n*cos(walkerParameters(1));
    arg_perigee_dot = (3/2)*gamma*n*(5*cos(walkerParameters(1))^2-1);
    M_dot = (3/2)*gamma*n*(3*cos(walkerParameters(1))^2-1);
    % Compute nominal orbital elements at t0
    OE = nominalOEWalkerAtAnchor(i,walkerParameters,semiMajorAxis,walkerAnchor,parallel);
    % Compute nominal orbital elements at t
    du = (n+M_dot+arg_perigee_dot)*(t-walkerAnchor(1));
    % Nominal u
    OE(2,:) = OE(2,:) + du;
    % Fix angle intervals
    OE(2,OE(2,:)>2*pi) = OE(2,OE(2,:)>2*pi)-floor(OE(2,OE(2,:)>2*pi)/(2*pi))*2*pi;
    OE(2,OE(2,:)<0) = OE(2,OE(2,:)<0)+ceil(-OE(2,OE(2,:)<0)/(2*pi))*2*pi;
    % Nominal Omega
    dOmega = Omega_dot*(t-walkerAnchor(1));
    OE(6,:) = OE(6,:) + dOmega;
    % Fix angle intervals
    OE(6,OE(6,:)>2*pi) = OE(6,OE(6,:)>2*pi)-floor(OE(6,OE(6,:)>2*pi)/(2*pi))*2*pi;
    OE(6,OE(6,:)<0) = OE(6,OE(6,:)<0)+ceil(-OE(6,OE(6,:)<0)/(2*pi))*2*pi;
end
