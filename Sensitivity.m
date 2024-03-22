function dObjFun =...
    Sensitivity(ini,FEM,U_stX,U_stY,U_stZ,E_stZ,...
    AntiEigF1,AntiEigF2,AntiEigF3,EigF1,EigF3,...
    AntiEigV1,AntiEigV2,AntiEigV3,EigV1,EigV3,...
    ro,drho_d,drho_h)
% Sensitivity Analysis

%% Check material model to decide which derivative factor apply
switch ini.MatModel
    case 'SIMP'
        derivFact = ini.p*ro.^(ini.p-1);
    case 'RAMP'
        derivFact = (1 + ini.p)./( 1+ini.p*(1-ro) ).^2;
    otherwise
        error('Invalid Material Model definition. It must be SIMP or RAMP')
end

%% Apply boundary conditions for BC-based approach
Dirichlet0 = [ini.Sets.BCdof; ini.Sets.SymmBCdof];   % Pre-allocation

% Modify K matrix due to Boundary Conditions
K0 = FEM.K;
K0(:,Dirichlet0) = 0;       % Put 0 in all DOFs rows restricted
K0(Dirichlet0,:) = 0;       % Put 0 in all DOFs columns restricted
diagK = diag(K0);           % Select the diagonal from K
diagK(Dirichlet0) = 1;      % Put ones for DOFs with restriction
K0 = spdiags(diagK,0,K0);   % Replace the diagonal of K with diagK

% Dirichlet Boundary Conditions for Mass Matrix
M0 = FEM.M;
M0(:,Dirichlet0) = 0;       % Put 0 in all DOFs rows restricted
M0(Dirichlet0,:) = 0;       % Put 0 in all DOFs columns restricted
diagM = diag(M0);           % Select the diagonal from M
diagM(Dirichlet0) = 1;      % Put ones for DOFs with restriction
M0 = spdiags(diagM,0,M0);   % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions for Mass Matrix
C0 = FEM.C;
C0(:,Dirichlet0) = 0;       % Put 0 in all DOFs rows restricted
C0(Dirichlet0,:) = 0;       % Put 0 in all DOFs columns restricted
diagM = diag(C0);           % Select the diagonal from M
diagM(Dirichlet0) = 1;      % Put ones for DOFs with restriction
C0 = spdiags(diagM,0,C0);   % Replace the diagonal of M with diagM

%% Apply boundary conditions for Eigen-based approach
% Dirichlet Boundary Conditions to Mass Matrix for antiresonances
Dirichlet1 = [ini.Sets.AntiRes1BCdof; ini.Sets.SymmBCdof];
M1 = FEM.M;
M1(:,Dirichlet1) = 0;       % Put 0 in all DOFs rows restricted
M1(Dirichlet1,:) = 0;       % Put 0 in all DOFs columns restricted
diagM = diag(M1);           % Select the diagonal from M
diagM(Dirichlet1) = 1;      % Put ones for DOFs with restriction
M1 = spdiags(diagM,0,M1);   % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions to Mass Matrix for antiresonances
Dirichlet2 = [ini.Sets.AntiRes2BCdof; ini.Sets.SymmBCdof];
M2 = FEM.M;
M2(:,Dirichlet2) = 0;       % Put 0 in all DOFs rows restricted
M2(Dirichlet2,:) = 0;       % Put 0 in all DOFs columns restricted
diagM = diag(M2);           % Select the diagonal from M
diagM(Dirichlet2) = 1;      % Put ones for DOFs with restriction
M2 = spdiags(diagM,0,M2);   % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions to Mass Matrix for antiresonances
Dirichlet3 = [ini.Sets.AntiRes3BCdof; ini.Sets.SymmBCdof];
M3 = FEM.M;
M3(:,Dirichlet3) = 0;       % Put 0 in all DOFs rows restricted
M3(Dirichlet3,:) = 0;       % Put 0 in all DOFs columns restricted
diagM = diag(M3);           % Select the diagonal from M
diagM(Dirichlet3) = 1;      % Put ones for DOFs with restriction
M3 = spdiags(diagM,0,M3);   % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions to Mass Matrix for resonances
Dirichlet4 = ini.Sets.SymmBCdof;
M4 = FEM.M;
M4(:,Dirichlet4) = 0;       % Put 0 in all DOFs rows restricted
M4(Dirichlet4,:) = 0;       % Put 0 in all DOFs columns restricted
diagM = diag(M4);           % Select the diagonal from M
diagM(Dirichlet4) = 1;      % Put ones for DOFs with restriction
M4 = spdiags(diagM,0,M4);   % Replace the diagonal of M with diagM

%% Preallocation
dFun = zeros(FEM.mesh.nel,1);
w1 = ini.w1; % Stresses
w2 = ini.w2;    w3 = ini.w3;    w4 = ini.w4; % Displacements
w5 = ini.w5;    w6 = ini.w6;    w7 = ini.w7; % Antiresonances
w8 = ini.w8;    w9 = ini.w9; % Eigenfrequencies

ux_dof = 3*ini.Sets.OutNodes - 2;
uy_dof = 3*ini.Sets.OutNodes - 1;
uz_dof = 3*ini.Sets.OutNodes - 0;

f_T1 = ini.Target1;
f_T2 = ini.Target2;
f_T3 = ini.Target3;
omega1 = 2*pi*f_T1;
omega2 = 2*pi*f_T2;
omega3 = 2*pi*f_T3;

% Dynamic Stiffness inverse
tic; fprintf('Inverse of matrix 1 ...');
if w2 ~= 0
    invS_dyn1 = inv(K0 - omega1^2*M0 + 1i*omega1*C0);
else
    invS_dyn1 = 0;
end
fprintf('Done.\n'); toc

tic; fprintf('Inverse of matrix 2 ...');
if w3 ~= 0
    if omega2 == omega1 && w2 ~= 0
        invS_dyn2 = invS_dyn1;
    else
        invS_dyn2 = inv(K0 - omega2^2*M0 + 1i*omega2*C0);
    end
else
    invS_dyn2 = 0;
end
fprintf('Done.\n'); toc

tic; fprintf('Inverse of matrix 3 ...');
if w4 ~= 0 || w1 ~= 0
    if omega3 == omega1 && w2 ~= 0
        invS_dyn3 = invS_dyn1;
    elseif omega3 == omega2 && w3 ~= 0
        invS_dyn3 = invS_dyn2;
    else
        invS_dyn3 = inv(K0 - omega3^2*M0 + 1i*omega3*C0);
    end
else 
    invS_dyn3 = 0;
end
fprintf('Done.\n'); toc


% Rearrange U vectors
if w2 ~= 0; U_stX = U_stX.'; U_stX = U_stX(:); else; U_stX = 0; end
if w3 ~= 0; U_stY = U_stY.'; U_stY = U_stY(:); else; U_stY = 0; end
if w4 ~= 0 || w1 ~= 0; U_stZ = U_stZ.'; U_stZ = U_stZ(:); else; U_stZ = 0; end


% For Antiresonance A1
if any(AntiEigF1)
    f_A1 = AntiEigF1;
    A = AntiEigV1.'; antiV1 = A(:);
    L_A1 = (2*pi*f_A1)^2; % Lambda for Antiresonance f_A1
else
    f_A1 = f_T1;
    antiV1 = [];
    L_A1 = 0;
end

% For Antiresonance A2
if any(AntiEigF2)
    f_A2 = AntiEigF2;
    A = AntiEigV2.'; antiV2 = A(:);
    L_A2 = (2*pi*f_A2)^2; % Lambda for Antiresonance f_A2
else
    f_A2 = f_T2;
    antiV2 = [];
    L_A2 = 0;
end

% For Antiresonance A3
if any(AntiEigF3)
    f_A3 = AntiEigF3;
    A = AntiEigV3.'; antiV3 = A(:);
    L_A3 = (2*pi*f_A3)^2; % Lambda for Antiresonance f_A2
else
    f_A3 = f_T3;
    antiV3 = [];
    L_A3 = 0;
end

% For the resonance terms
if any(EigF1)
    f_R1 = EigF1;
    modes1 = length(EigF1);
    Lambda1 = (2*pi*EigF1).^2;
    % Rearrange Eigenvector
    V1 = zeros(FEM.mesh.dof,modes1);
    for q = 1:modes1
        A = EigV1(:,:,q).';  V1(:,q) = A(:);
    end
else
    f_R1 = 0;
    modes1 = 0;
    Lambda1 = 0;
    V1 = 0;
end

% if any(EigF2)
%     f_R2 = EigF2;
%     modes2 = length(EigF2);
%     Lambda2 = (2*pi*EigF2).^2;
%     Rearrange Eigenvector
%     V2 = zeros(FEM.mesh.dof,modes2);
%     for q = 1:modes2
%         A = EigV2(:,:,q).';  V2(:,q) = A(:);
%     end
% else
%     f_R2 = 0;
%     modes2 = 0;
%     Lambda2 = 0;
%     V2 = 0;
% end

if any(EigF3)
    f_R3 = EigF3;
    modes3 = length(EigF3);
    Lambda3 = (2*pi*EigF3).^2;
    % Rearrange Eigenvector
    V3 = zeros(FEM.mesh.dof,modes3);
    for q = 1:modes3
        A = EigV3(:,:,q).';  V3(:,q) = A(:);
    end
else
    f_R3 = 0;
    modes3 = 0;
    Lambda3 = 0;
    V3 = 0;
end

DesignElements = ...
 setdiff(1:FEM.mesh.nel,[ini.Sets.SolidNonDesign ini.Sets.VoidNonDesign]);

for k = 1:FEM.mesh.nel
if ismember(k,DesignElements) % to skip computation of non-design elments
%% Derivative of matrices w.r.t to element k
% Stiffness matrix derivative
dK = derivFact(k)*...
    sparse(FEM.mesh.Ai(:,k),FEM.mesh.Aj(:,k),FEM.k(:,k),...
    FEM.mesh.dof,FEM.mesh.dof); %#ok

% Mass matrix derivative
dM = derivFact(k)*...
    sparse(FEM.mesh.Ai(:,k),FEM.mesh.Aj(:,k),FEM.m(:,k),...
    FEM.mesh.dof,FEM.mesh.dof);

% Damping matrix derivative
dC = derivFact(k)*...
    sparse(FEM.mesh.Ai(:,k),FEM.mesh.Aj(:,k),FEM.c(:,k),...
    FEM.mesh.dof,FEM.mesh.dof);

%% Apply boundary conditions to matrix derivatives dK, dM and dC
% Dirichlet Boundary Conditions for Stiffness Matrix
dK0 = dK;
dK0(:,Dirichlet0) = 0;      % Put 0 in all DOFs rows restricted
dK0(Dirichlet0,:) = 0;      % Put 0 in all DOFs columns restricted
diagK = diag(dK0);          % Select the diagonal from K
diagK(Dirichlet0) = 1;      % Put ones for DOFs with restriction
dK0 = spdiags(diagK,0,dK0); % Replace the diagonal of K with diagK

% Dirichlet Boundary Conditions for Mass Matrix
dM0 = dM;
dM0(:,Dirichlet0) = 0;      % Put 0 in all DOFs rows restricted
dM0(Dirichlet0,:) = 0;      % Put 0 in all DOFs columns restricted
diagM = diag(dM0);          % Select the diagonal from M
diagM(Dirichlet0) = 1;      % Put ones for DOFs with restriction
dM0 = spdiags(diagM,0,dM0); % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions for Damping Matrix
dC0 = dC;
dC0(:,Dirichlet0) = 0;      % Put 0 in all DOFs rows restricted
dC0(Dirichlet0,:) = 0;      % Put 0 in all DOFs columns restricted
diagM = diag(dC0);          % Select the diagonal from M
diagM(Dirichlet0) = 1;      % Put ones for DOFs with restriction
dC0 = spdiags(diagM,0,dC0); % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions for dK for antiresonance #1
dK1 = dK;
dK1(:,Dirichlet1) = 0;      % Put 0 in all DOFs rows restricted
dK1(Dirichlet1,:) = 0;      % Put 0 in all DOFs columns restricted
diagK = diag(dK1);          % Select the diagonal from K
diagK(Dirichlet1) = 1;      % Put ones for DOFs with restriction
dK1 = spdiags(diagK,0,dK1); % Replace the diagonal of K with diagK

% Dirichlet Boundary Conditions for dM for antiresonance #1
dM1 = dM;
dM1(:,Dirichlet1) = 0;      % Put 0 in all DOFs rows restricted
dM1(Dirichlet1,:) = 0;      % Put 0 in all DOFs columns restricted
diagM = diag(dM1);          % Select the diagonal from M
diagM(Dirichlet1) = 1;      % Put ones for DOFs with restriction
dM1 = spdiags(diagM,0,dM1); % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions for dK for antiresonance #2
dK2 = dK;
dK2(:,Dirichlet2) = 0;      % Put 0 in all DOFs rows restricted
dK2(Dirichlet2,:) = 0;      % Put 0 in all DOFs columns restricted
diagK = diag(dK2);          % Select the diagonal from K
diagK(Dirichlet2) = 1;      % Put ones for DOFs with restriction
dK2 = spdiags(diagK,0,dK2); % Replace the diagonal of K with diagK

% Dirichlet Boundary Conditions for dM for antiresonance #2
dM2 = dM;
dM2(:,Dirichlet2) = 0;      % Put 0 in all DOFs rows restricted
dM2(Dirichlet2,:) = 0;      % Put 0 in all DOFs columns restricted
diagM = diag(dM2);          % Select the diagonal from M
diagM(Dirichlet2) = 1;      % Put ones for DOFs with restriction
dM2 = spdiags(diagM,0,dM2); % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions for dK for antiresonance #3
dK3 = dK;
dK3(:,Dirichlet3) = 0;      % Put 0 in all DOFs rows restricted
dK3(Dirichlet3,:) = 0;      % Put 0 in all DOFs columns restricted
diagK = diag(dK3);          % Select the diagonal from K
diagK(Dirichlet3) = 1;      % Put ones for DOFs with restriction
dK3 = spdiags(diagK,0,dK3); % Replace the diagonal of K with diagK

% Dirichlet Boundary Conditions for dM for antiresonance #3
dM3 = dM;
dM3(:,Dirichlet3) = 0;      % Put 0 in all DOFs rows restricted
dM3(Dirichlet3,:) = 0;      % Put 0 in all DOFs columns restricted
diagM = diag(dM3);          % Select the diagonal from M
diagM(Dirichlet3) = 1;      % Put ones for DOFs with restriction
dM3 = spdiags(diagM,0,dM3); % Replace the diagonal of M with diagM

% Dirichlet Boundary Conditions for dK for resonances
dK4 = dK;
dK4(:,Dirichlet4) = 0;      % Put 0 in all DOFs rows restricted
dK4(Dirichlet4,:) = 0;      % Put 0 in all DOFs columns restricted
diagK = diag(dK4);          % Select the diagonal from K
diagK(Dirichlet4) = 1;      % Put ones for DOFs with restriction
dK4 = spdiags(diagK,0,dK4); % Replace the diagonal of K with diagK

% Dirichlet Boundary Conditions for dM for resonances
dM4 = dM;
dM4(:,Dirichlet4) = 0;      % Put 0 in all DOFs rows restricted
dM4(Dirichlet4,:) = 0;      % Put 0 in all DOFs columns restricted
diagM = diag(dM4);          % Select the diagonal from M
diagM(Dirichlet4) = 1;      % Put ones for DOFs with restriction
dM4 = spdiags(diagM,0,dM4); % Replace the diagonal of M with diagM

%% Derivate of displacements for a single frequency
if w2 ~= 0
    dS_dyn = dK0 - omega1^2*dM0 + 1i*omega1*dC0;
    dU_stX = -(invS_dyn1*dS_dyn)*U_stX;
    dUx = dU_stX(ux_dof); % Separate desired data
else
    dUx = 0;
end
if w3 ~= 0
    dS_dyn = dK0 - omega2^2*dM0 + 1i*omega2*dC0;
    dU_stY = -(invS_dyn2*dS_dyn)*U_stY;
    dUy = dU_stY(uy_dof); % Separate desired data
else
    dUy = 0;
end
if w4 ~= 0 && w1 ~= 0
    dS_dyn = dK0 - omega3^2*dM0 + 1i*omega3*dC0;
    dU_stZ = -(invS_dyn3*dS_dyn)*U_stZ;
    dUz = dU_stZ(uz_dof); % Separate desired data
    % Derivative of stress (and strain if needed)
    [dStress,~,~] = dStressStrain(FEM.mesh,ini,dU_stZ,ro,E_stZ,k);
    dSzz = dStress(3,ini.Sets.OutElements).'; % for selected elements
elseif w4 ~= 0 && w1 == 0
    dS_dyn = dK0 - omega3^2*dM0 + 1i*omega3*dC0;
    dU_stZ = -(invS_dyn3*dS_dyn)*U_stZ;
    dUz = dU_stZ(uz_dof); % Separate desired data
    dSzz = 0;
elseif w4 == 0 && w1 ~= 0
    dS_dyn = dK0 - omega3^2*dM0 + 1i*omega3*dC0;
    dU_stZ = -(invS_dyn3*dS_dyn)*U_stZ;
    dUz = dU_stZ(uz_dof); % Separate desired data
    % Derivative of stress (and strain if needed)
    [dStress,~,~] = dStressStrain(FEM.mesh,ini,dU_stZ,ro,E_stZ,k);
    dSzz = dStress(3,ini.Sets.OutElements).'; % for selected elements
else
    dUz = 0;
    dSzz = 0;
end

%% Sensitiviy with respect to element k
if any(FEM.Szz)
    dNormSzz = (real(FEM.Szz)'*real(dSzz) + imag(FEM.Szz)'*imag(dSzz))/...
           ( norm(FEM.Szz,2) );
else
    dNormSzz = 0;
end

if any(FEM.Ux)
    dNormUx = ( real(FEM.Ux)'*real(dUx) + imag(FEM.Ux)'*imag(dUx) )/...
          ( norm(FEM.Ux,2) );
else
    dNormUx = 0;
end

if any(FEM.Uy)
    dNormUy = ( real(FEM.Uy)'*real(dUy) + imag(FEM.Uy)'*imag(dUy) )/...
          ( norm(FEM.Uy,2) );
else
    dNormUy = 0;
end

if any(FEM.Uz)
    dNormUz = ( real(FEM.Uz)'*real(dUz) + imag(FEM.Uz)'*imag(dUz) )/...
          ( norm(FEM.Uz,2) );
else
    dNormUz = 0;
end

%% Derivate of eigenvalues w.r.t element k
% For antiresonance f_A1
if any(AntiEigF1)
    dL_A1 = ( antiV1'*(dK1 - L_A1*dM1)*antiV1)/( antiV1'*M1*antiV1 );
else
    dL_A1 = 0;
end

% For antiresonance f_A2
if any(AntiEigF2)
    dL_A2 = ( antiV2'*(dK2 - L_A2*dM2)*antiV2)/( antiV2'*M2*antiV2 );
else
    dL_A2 = 0;
end

% For antiresonance f_A3
if any(AntiEigF3)
    dL_A3 = ( antiV3'*(dK3 - L_A3*dM3)*antiV3)/( antiV3'*M3*antiV3 );
else
    dL_A3 = 0;
end

% for resonances R1
if any(EigF1)
    dL_R1 = zeros(modes1,1);
    for q = 1:modes1
        dL_R1(q) = (V1(:,q)'*(dK4 - Lambda1(q)*dM4)*V1(:,q))/...
                     (V1(:,q)'*M4*V1(:,q)); %#ok
    end
else
    dL_R1 = 0;
end

% for resonances R3
if any(EigF3)
    dL_R3 = zeros(modes3,1);
    for q = 1:modes3
        dL_R3(q) = (V3(:,q)'*(dK4 - Lambda3(q)*dM4)*V3(:,q))/...
                     (V3(:,q)'*M4*V3(:,q)); %#ok
    end
else
    dL_R3 = 0;
end

%% Sensitiviy w.r.t element k
if any(AntiEigF1)
    TermA1 =  ( (f_A1-f_T1)/(4*pi^2*f_A1*f_T1^2) )*dL_A1;
else
    TermA1 = 0;
end

if any(AntiEigF2)
    TermA2 =  ( (f_A2-f_T2)/(4*pi^2*f_A2*f_T2^2) )*dL_A2;
else
    TermA2 = 0;
end

if any(AntiEigF3)
    TermA3 =  ( (f_A3-f_T3)/(4*pi^2*f_A3*f_T3^2) )*dL_A3;
else
    TermA3 = 0;
end

if any(EigF1)
    TermRes1 = ((f_A1)./(8*pi^2.*(f_R1-f_A1).^2)).*...
                sqrt(((f_R1-f_A1)/f_A1).^2).*...
                ( (1/f_A1)*dL_A1 -...
                (f_A1.*dL_R1)./(f_R1.*(f_R1-f_A1)) + ...
                dL_A1./(f_R1-f_A1) );
else
    TermRes1 = 0;
end

if any(EigF3)
    TermRes3 = ((f_A3)./(8*pi^2.*(f_R3-f_A3).^2)).*...
                sqrt(((f_R3-f_A3)/f_A3).^2).*...
                ( (1/f_A3)*dL_A3 -...
                (f_A3.*dL_R3)./(f_R3.*(f_R3-f_A3)) +...
                dL_A3./(f_R3-f_A3) );
else
    TermRes3 = 0;
end

dFun(k) = w1*dNormSzz + ...
          w2*dNormUx + w3*dNormUy + w4*dNormUz + ...
          w5*TermA1 + w6*TermA2 + w7*TermA3 +...
          w8*sum(TermRes1) + w9*sum(TermRes3);
% fprintf('Sensitivity for element %d completed \n', k)

else
    dFun(k) = 0;
    fprintf('Sensitivity for non-design element %d set to 0 \n', k)
end

end % End of sensitivity parfor loop

% Chain rule derivative
dObjFun = sum(dFun.*drho_h.*drho_d)';

end % End of sensitivity function