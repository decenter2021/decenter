function OE = nominalOEWalkerAtAnchor(idx,walkerParameters,sma,walkerAnchor,parallel)
%% Notation of input parameters for clarity
inclination = walkerParameters(1);
T = walkerParameters(2); % (T) Number of satellites
P = walkerParameters(3); % (P) Number of orbital planes
F = walkerParameters(4); % (F) Phasing parameter
numberOfSatellitesPerPlane = T/P;
eccentricity = 0; 
argumentOfPeriapsis = 0;

%% Compute parameters
longitudeOfAscendingNodeSpacing = 2*pi/P;
trueAnomalySpacing = 2*pi/numberOfSatellitesPerPlane;
adjacentPlaneSpacing = F*2*pi/T;

%% Compute OE
OE = zeros(6,length(idx));
OE(1,:) = sma; % a (m)
OE(3,:) = eccentricity*cos(argumentOfPeriapsis); % ex
OE(4,:) = eccentricity*sin(argumentOfPeriapsis); % ey
OE(5,:) = inclination;

if parallel
    parfor i = 1:length(idx)
        OE(2,i) = walkerAnchor(2)+(rem(idx(i)-1,numberOfSatellitesPerPlane))*trueAnomalySpacing...
        + floor((idx(i)-1)/numberOfSatellitesPerPlane)*adjacentPlaneSpacing; % u (rad)
    end
    parfor i = 1:N    
        OE(6,i) = walkerAnchor(3)+floor((idx(i)-1)/numberOfSatellitesPerPlane)*longitudeOfAscendingNodeSpacing;% Omega (rad)
    end
else
    for i = 1:length(idx)
        OE(2,i) = walkerAnchor(2)+(rem(idx(i)-1,numberOfSatellitesPerPlane))*trueAnomalySpacing...
        + floor((idx(i)-1)/numberOfSatellitesPerPlane)*adjacentPlaneSpacing; % u (rad)
        OE(6,i) = walkerAnchor(3)+floor((idx(i)-1)/numberOfSatellitesPerPlane)*longitudeOfAscendingNodeSpacing;% Omega (rad)
    end  
end
%% Fix angle intervals
OE(2,OE(2,:)>2*pi) = OE(2,OE(2,:)>2*pi)-floor(OE(2,OE(2,:)>2*pi)/(2*pi))*2*pi;
OE(2,OE(2,:)<0) = OE(2,OE(2,:)<0)+ceil(-OE(2,OE(2,:)<0)/(2*pi))*2*pi;
OE(6,OE(6,:)>2*pi) = OE(6,OE(6,:)>2*pi)-floor(OE(6,OE(6,:)>2*pi)/(2*pi))*2*pi;
OE(6,OE(6,:)<0) = OE(6,OE(6,:)<0)+ceil(-OE(6,OE(6,:)<0)/(2*pi))*2*pi;
% if OE(2)>2*pi
%     OE(2) = OE(2)-floor(OE(2)/(2*pi))*2*pi;
% elseif OE(2)<0
%     OE(2) = OE(2)+ceil(-OE(2)/(2*pi))*2*pi;
% end
% if OE(6)>2*pi
%     OE(6) = OE(6)-floor(OE(6)/(2*pi))*2*pi;
% elseif OE(6)<0
%     OE(6) = OE(6)+ceil(-OE(6)/(2*pi))*2*pi;
% end  

end