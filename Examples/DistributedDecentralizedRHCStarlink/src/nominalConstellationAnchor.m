function [u0,Omega0] = nominalConstellationAnchor(OE_u_Omega,walkerParameters)
    % Inputs
    T = walkerParameters(2); % (T) Number of satellites
    P = walkerParameters(3); % (P) Number of orbital planes
    F = walkerParameters(4); % (F) Phasing parameter
    u0 = 0;
    angleRange_u0 = 0;
    Omega0 = 0;
    angleRange_Omega0 = 0;
    for i = 1:T
        % Compute raw angle
        error = OE_u_Omega(1,i)- (rem(i-1,T/P)*2*pi*P/T + floor((i-1)*P/T)*2*pi*F/T);
        % The first satellite is used to select the angle range to be
        % considered
        if i == 1
            % Set range to -pi, pi
            if abs(error - floor((error+pi)/(2*pi))*2*pi) < pi/2
                angleRange_u0 = 0;              
            else % Set range to 0, 2*pi
                angleRange_u0 = 1;
            end
        end
        % Fix error range
        if ~angleRange_u0 % Range in -pi, pi
            if error < -pi
                error = error + (floor(abs(error-pi)/(2*pi)))*2*pi;
            elseif error > pi
                error = error - floor((error+pi)/(2*pi))*2*pi;
            end    
        else % Range in 0,2*pi
            if error >2*pi
                error = error-floor(error/(2*pi))*2*pi;
            elseif error<0
                error = error+ceil(-error/(2*pi))*2*pi;
            end
        end
        % Contribution of this satellite to the solution
        u0 = u0 + error/T; % u (rad)

        % Compute raw angle
        error = OE_u_Omega(2,i)- floor((i-1)*P/T)*2*pi/P;
        % The first satellite is used to select the angle range to be
        % considered
        if i == 1
            % Set range to -pi, pi
            if abs(error - floor((error+pi)/(2*pi))*2*pi) < pi/2
                angleRange_Omega0 = 0;              
            else % Set range to 0, 2*pi
                angleRange_Omega0 = 1;
            end
        end
        % Fix error range
        if ~angleRange_Omega0 % Range in -pi, pi
            if error < -pi
                error = error + (floor(abs(error-pi)/(2*pi)))*2*pi;
            elseif error > pi
                error = error - floor((error+pi)/(2*pi))*2*pi;
            end    
        else % Range in 0,2*pi
            if error >2*pi
                error = error-floor(error/(2*pi))*2*pi;
            elseif error<0
                error = error+ceil(-error/(2*pi))*2*pi;
            end
        end
        % Contribution of this satellite to the solution
        Omega0 = Omega0 + error/T;  % Omega (rad)
    end 

    % Fix angles of final solution
    if u0>2*pi
        u0 = u0-floor(u0/(2*pi))*2*pi;
    elseif u0<0
        u0 = u0+ceil(-u0/(2*pi))*2*pi;
    end
    if Omega0>2*pi
        Omega0 = Omega0-floor(Omega0/(2*pi))*2*pi;
    elseif Omega0<0
        Omega0 = Omega0+ceil(-Omega0/(2*pi))*2*pi;
    end

    % We could just program this for a generic component: a function that
    % receives as arguments the data and a function hanle to compute
    % nominal OE as a function of i - not worth it yet

end