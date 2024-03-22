function [LB, UB] = MoveLimits(ini,FEM,ro)

%% Move limits
LB = ro - ini.step;    % Lower Move Limit
UB = ro + ini.step;    % Upper Move Limit
fprintf('Move limits step size: %d \n',ini.step)

%% Volume analysis
% prevent non-feasible solutions/relax the solution space
Volume = (FEM.mesh.Ve')*ro;
VolumeFraction = Volume/FEM.TotalVolume;
sigma = 2; % factor to expand limits

% If the volumen is higher than the maximum limit:
% Expand the lower limit to allow for lower density elements
if VolumeFraction > (ini.MaxVol - ini.step)
    LB = ro - sigma*ini.step; % Lower Limit
    fprintf('Lower move limit expanded \n');
end

% If the volume is lower than the minimum limit:
% Expand the upper limit to allow for higher density elements
if VolumeFraction < (ini.MinVol + ini.step)
    UB = ro + sigma*ini.step; % Upper Limit
    fprintf('Upper move limit expanded \n');
end

%% Check limits
LB(LB<1e-9) = 1e-9;     UB(UB>1) = 1;

% Protect the limits in extreme situations
UB(UB<ini.step) = ini.step;     LB(LB>1-ini.step) = 1-ini.step;

%% Set limits for forced elements
if any(ini.VoidNonDesign)
    LB(ini.Sets.VoidNonDesign) = 1e-9;
    UB(ini.Sets.VoidNonDesign) = 1e-9 + ini.step;
end
if any(ini.SolidNonDesign)
    LB(ini.Sets.SolidNonDesign) = 1 - ini.step;
    UB(ini.Sets.SolidNonDesign) = 1;
end
fprintf('Move limits defined \n\n')