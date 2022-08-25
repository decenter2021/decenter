function [Du,DOmega,anchor] = constellationRigidity(OE_u_Omega, walkerParameters)
    % Inputs
    T = walkerParameters(2); % (T) Number of satellites
    P = walkerParameters(3); % (P) Number of orbital planes
    F = walkerParameters(4); % (F) Phasing parameter
    % Compute optimal anchor at this instant
    [u0, Omega0] = nominalConstellationAnchor(OE_u_Omega,walkerParameters);
    Du = zeros(T,1);
    DOmega = zeros(T,1);
    for i = 1:T
        error = OE_u_Omega(1,i)- (u0 + rem(i-1,T/P)*2*pi*P/T + floor((i-1)*P/T)*2*pi*F/T);       
        if error < -pi
            error = error + (floor(abs(error-pi)/(2*pi)))*2*pi;
        elseif error > pi
            error = error - floor((error+pi)/(2*pi))*2*pi;
        end
        %Du = Du + abs(error)/T;
        Du(i) = error;
        error = OE_u_Omega(2,i)- (Omega0 + floor((i-1)*P/T)*2*pi/P);
        if error < -pi
            error = error + (floor(abs(error-pi)/(2*pi)))*2*pi;
        elseif error > pi
            error = error - floor((error+pi)/(2*pi))*2*pi;
        end
        %DOmega = DOmega + abs(error)/T;
        DOmega(i) = error;
    end
    anchor = [u0; Omega0];
end

