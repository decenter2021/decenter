% OE to Relative OE
function [dalpha] = OEMean2OERel(OE_d,OE_c)    
    %% Compute dalpha 
    % Eq (1) Mauro et al. (2018)
    dalpha = zeros(6,1);
    dalpha(1) = (OE_d(1)/OE_c(1))-1; % a
    du = OE_d(2)-OE_c(2);
    if du < -pi
        du = du + 2*pi;
    elseif du > pi
        du = du - 2*pi;
    end
    dOmega = OE_d(6)-OE_c(6);
    if dOmega < -pi
        dOmega = dOmega + 2*pi;
    elseif dOmega > pi
        dOmega = dOmega - 2*pi;
    end 
    dalpha(2) = du+(dOmega)*cos(OE_c(5)); % lambda
    dalpha(3) = OE_d(3)-OE_c(3); % ex
    dalpha(4) = OE_d(4)-OE_c(4); % ey
    dalpha(5) = OE_d(5)-OE_c(5); % ix
    dalpha(6) = (OE_d(6)-OE_c(6))*sin(OE_c(5)); % iy
end