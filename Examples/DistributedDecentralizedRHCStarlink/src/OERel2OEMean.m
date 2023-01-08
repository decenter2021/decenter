% Relative OE to OE
function [OE_d] = OERel2OEMean(dalpha,OE_c)    
    %% Compute dalpha 
    % Eq (8) Mauro et al. (2018)
    OE_d = zeros(6,1);
    OE_d(6) = OE_c(6)+dalpha(6)/sin(OE_c(5)); % Omega
    OE_d(1) = (1+dalpha(1))*OE_c(1); % a
    OE_d(2) = OE_c(2) + dalpha(2) - (OE_d(6)-OE_c(6))*cos(OE_c(5)); % u
    OE_d(3) = OE_c(3) + dalpha(3);
    OE_d(4) = OE_c(4) + dalpha(4);
    OE_d(5) = OE_c(5) + dalpha(5);
    
    u = OE_d(2);
    if u>2*pi
        OE_d(2) = u-floor(u/(2*pi))*2*pi;
    elseif u<0
        OE_d(2) = u+ceil(-u/(2*pi))*2*pi;
    end
    
    Omega = OE_d(6);
    if Omega>2*pi
        OE_d(6) = Omega-floor(Omega/(2*pi))*2*pi;
    elseif Omega<0
        OE_d(6) = Omega+ceil(-Omega/(2*pi))*2*pi;
    end  
    
end