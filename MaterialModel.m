function [roK,roM,roC] = MaterialModel(MatModel,ro,pK,pM,pC,limit)
% Material model interpolation

if nargin == 3
    pM = pK;
    pC = pK;
    limit = 1e-9;
elseif nargin == 4
    pC = pM;
    limit = 1e-9;
elseif nargin == 5
    limit = 1e-9;
end

% Interpolate material properties
switch MatModel
    case 'SIMP'
        roK = ro.^pK;
        roM = ro.^pM;
        roC = ro.^pC;
        
    case 'RAMP'
        roK = ro./(1 + pK*(1-ro));
        roM = ro./(1 + pM*(1-ro));
        roC = ro./(1 + pC*(1-ro));

    case 'testRAMP'
        roK = ro./(1 + pK*(1-ro));
        roM = ro./(1 + pM*(1-ro));

        idVoid = roK <= 0.1;
        roK(idVoid) = roK(idVoid)./( 10*(1e-6.^roK(idVoid)) );
        idVoid = roM <= 0.1;
        roM(idVoid) = roM(idVoid)./( 10*(1e-6.^roM(idVoid)) );
        % roK(idVoid) = ro(idVoid)./100;
        % roK(idVoid) = 0.01*( 1./roK(idVoid) );
        % roK(idVoid) = 0.002*( 1./(ro(idVoid).^2) );
        
    otherwise
        error('Invalid Material Model definition. It must be SIMP or RAMP')
end

% Set density lower limit
roK(roK<limit) = limit;
roM(roM<limit) = limit;
roC(roC<limit) = limit;