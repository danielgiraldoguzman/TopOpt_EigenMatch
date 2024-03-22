function [dStress,dStrain,eps] = dStressStrain(mesh,ini,dU,ro,E,k)

% pre-allocation
dStress = zeros(6,mesh.nel);
dStrain = zeros(6,mesh.nel);
eps = zeros(6,8*mesh.nel);

for e = 1:mesh.nel
    enodes = mesh.conect(e,:);   % Nodal idex for element i
    xyz = mesh.ncoord(enodes,:); % Nodal coordinates
    DOFs = [3*enodes - 2;
            3*enodes - 1;
            3*enodes];
    edof = DOFs(:); % Element DOFs
    du = dU(edof);
    
    % Stress-strain computation
    switch ini.ElementType
        case 8
            dStrain(:,e) = StrainBrick8(xyz,du);
        case 20
            [dStrain(:,e),eps(:,(e-1)*8+1:e*8)] = StrainBrick20(xyz,du);
    end
    
    %drho_e/drho_k
    ek = e == k;
    
    switch ini.MatModel
    case 'SIMP'
        MatFactor = ro(e)^ini.p;
        derivFact = ini.p*ro(e)^(ini.p-1);
    case 'RAMP'
        MatFactor = ro(e)/(1 + ini.p*(1-ro(e)));
        derivFact = (1+ini.p)./(1 + ini.p*( 1-ro(e) ))^2;
    end
    
    % Stress
    dStress(:,e) = derivFact*ek*ini.cE*E(e,:).' + ...
                   MatFactor*ini.cE*dStrain(:,e);
end

end