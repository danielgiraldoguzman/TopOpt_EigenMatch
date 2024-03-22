% Write updated matrices to ABAQUS
fprintf('Writing updated FE Matrices...  ')
[FEM, MatFactorK] = WriteUpdatedMTX(ini,rho_h(:,it),FEM);
fprintf('Done\n\n')

% The input files must be modified as needed before running an optimization.
% Make sure to verify the matrix information transferred from and to
% ABAQUS, especially if the mesh has been changed.

%% BCx simulation
if ini.w2 ~= 0
Inp_BC = ['BC_Ux_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_BC ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('Simulation for BC done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
[U_stX,~,~] = ... % values subject to X-load
    ReadFILResults(Inp_BC,ini.recordKeysBC,ini.ElementType,FEM.mesh.nel);
fprintf('Results file processed.\n\n')
% Extract data
FEM.Ux = U_stX(ini.Sets.OutNodes,1);
else
    U_stX = 0;
    FEM.Ux = 0;
end

%% BCy simulation
if ini.w3 ~= 0
Inp_BC = ['BC_Uy_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_BC ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('Simulation for BC done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
[U_stY,~,~] = ... % values subject to Y-load
    ReadFILResults(Inp_BC,ini.recordKeysBC,ini.ElementType,FEM.mesh.nel);
fprintf('Results file processed.\n\n')
% Extract data
FEM.Uy = U_stY(ini.Sets.OutNodes,2);
else
    U_stY = 0;
    FEM.Uy = 0;
end

%% BCZ simulation
if ini.w4 ~= 0 && ini.w1 ~= 0  % Both BC-zz and Antires z-axis active
Inp_BC = ['BC_Uz_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_BC ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('Simulation for BC done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
[U_stZ,S_stZ,E_stZ] = ... % values subject to Z-load
    ReadFILResults(Inp_BC,ini.recordKeysBC,ini.ElementType,FEM.mesh.nel);
% Ee_stZ = Avg8IntPoints(E_stZ,ini);  % Average integration points
Ss_stZ = Avg8IntPoints(S_stZ,ini);  % Average integration points
Ss_stZ = MatFactorK.*Ss_stZ;        % Apply material-model penalization
fprintf('Results file processed.\n\n')
% Extract data
FEM.Uz = U_stZ(ini.Sets.OutNodes,3);
FEM.Szz = Ss_stZ(ini.Sets.OutElements,3);

elseif ini.w4 ~= 0 && ini.w1 == 0  % Only AntiRes z-axis active
Inp_BC = ['BC_Uz_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_BC ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('Simulation for BC done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
U_stZ = ... % values subject to Z-load
    ReadFILResults(Inp_BC,ini.recordKeysBC,ini.ElementType,FEM.mesh.nel);
fprintf('Results file processed.\n\n')
% Extract data
FEM.Uz = U_stZ(ini.Sets.OutNodes,3);
E_stZ = 0;
FEM.Szz = 0;

elseif ini.w4 == 0 && ini.w1 ~= 0 % Only BC-zz active
Inp_BC = ['BC_Uz_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_BC ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('Simulation for BC done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
[U_stZ,S_stZ,E_stZ] = ... % values subject to Z-load
    ReadFILResults(Inp_BC,ini.recordKeysBC,ini.ElementType,FEM.mesh.nel);
% Ee_stZ = Avg8IntPoints(E_stZ,ini);  % Average integration points
Ss_stZ = Avg8IntPoints(S_stZ,ini);  % Average integration points
Ss_stZ = MatFactorK.*Ss_stZ;        % Apply material-model penalization
fprintf('Results file processed.\n\n')
% Extract data
FEM.Uz = U_stZ(ini.Sets.OutNodes,3);
FEM.Szz = Ss_stZ(ini.Sets.OutElements,3);

else
    U_stZ = 0;
    FEM.Uz = 0;
    E_stZ = 0;
    FEM.Szz = 0;
end

%% Harmonic FEA for x-axis force
if ini.w5 ~= 0
Inp_Harm1 = ['FreqResX_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_Harm1 ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('FRF for x-axis forces done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Extract frequency vector
if isempty(Freq)
    Freq = ReadFreq(Inp_Harm1);
    ObjFun.FRF1 = zeros(length(Freq),ini.MaxIter);
    fprintf('ODB file processed.\n');
end
% Solution
fprintf('Reading results file...  ')
tic; HarmDisp1 = ReadFILResults(Inp_Harm1,ini.recordKeysEig);
fprintf('Done\n'); toc

% Extract FRF on the Aniresonance nodes only
FRF1 = abs( HarmDisp1(ini.AntiResBCnodes,:,:) );
Uabs1 = reshape( FRF1(:,3-ini.DOFkey1,:),length(ini.AntiResBCnodes),[] );
Uavg1 = mean(Uabs1,1)';
ObjFun.FRF1(1:length(Uavg1),it) = Uavg1;
fprintf('FRF for x-axis force processed.\n\n')
end

%% Harmonic FEA for y-axis force
if ini.w6 ~= 0
Inp_Harm2 = ['FreqResY_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_Harm2 ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('FRF for y-axis forces done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Extract frequency vector
if isempty(Freq)
    Freq = ReadFreq(Inp_Harm2);
    ObjFun.FRF2 = zeros(length(Freq),ini.MaxIter);
    fprintf('ODB file processed.\n');
end

% Solution
fprintf('Reading results file...  ')
tic; HarmDisp2 = ReadFILResults(Inp_Harm2,ini.recordKeysEig);
fprintf('Done\n'); toc

% Extract FRF on the Aniresonance nodes only
FRF2 = abs( HarmDisp2(ini.AntiResBCnodes,:,:) );
Uabs2 = reshape( FRF2(:,3-ini.DOFkey2,:),length(ini.AntiResBCnodes),[] );
Uavg2 = mean(Uabs2,1)';
ObjFun.FRF2(1:length(Uavg2),it) = Uavg2;
fprintf('FRF for y-axis force processed.\n\n')
end

%% Harmonic FEA for z-axis force
if ini.w7 ~= 0
Inp_Harm3 = ['FreqResZ_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_Harm3 ini.cpus ' interactive ask_delete=OFF']);
if SysStat == 0; fprintf('FRF for z-axis forces done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Extract frequency vector
if isempty(Freq)
    Freq = ReadFreq(Inp_Harm3);
    ObjFun.FRF3 = zeros(length(Freq),ini.MaxIter);
    fprintf('ODB file processed.\n');
end

% Solution
fprintf('Reading results file...  ')
tic; HarmDisp3 = ReadFILResults(Inp_Harm3,ini.recordKeysEig);
fprintf('Done\n'); toc

% Extract FRF on the Aniresonance nodes only
FRF3 = abs( HarmDisp3(ini.AntiResBCnodes,:,:) );
Uabs3 = reshape( FRF3(:,3-ini.DOFkey3,:),length(ini.AntiResBCnodes),[] );
Uavg3 = mean(Uabs3,1)';
ObjFun.FRF3(1:length(Uavg3),it) = Uavg3;
fprintf('FRF for z-axis force processed.\n\n')
end

%% Modal FEA-1 for antiresonance eigenfrequency f_Ax
% Run ABAQUS simulation for antiresonance eigenfrequencies
if ini.w5 ~= 0
Inp_AntiRes1 = ['AntiEigFX_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_AntiRes1 ' cpus=20 interactive ask_delete=OFF']);
if SysStat == 0; fprintf('AntiresonanceEigenFreqX simulation done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
tic; 
[AntiEigV1,~,~,AntiEigF1] = ReadFILResults(Inp_AntiRes1,ini.recordKeysEig);
AntiEigF1 = AntiEigF1( 1:size(AntiEigV1,3) ); % Cut EigF if longer than EigV
AntiEigV1(:,:,AntiEigF1<1) = []; % Delete rigid body modes (1Hz threshold)
AntiEigF1(AntiEigF1<1) = [];     % Delete rigid body modes (1Hz threshold)
fprintf('Results file loaded.\n'); toc

% Remove localized modes
idx_Localized = LocalizedModes(ini,FEM,rho_h(:,it),AntiEigV1);
AntiEigF1(idx_Localized) = [];
AntiEigV1(:,:,idx_Localized) = [];
fprintf('Localized modes removed: %d.\n', length(idx_Localized))

else
    AntiEigF1 = [];
    AntiEigV1 = [];
end

%% Antiresonance peak identification for f_Ax
if ini.w5 ~= 0
    [Peak,idP,Width,Prom] = findpeaks(1./Uavg1);
    Prox = abs(ini.Target1 - Freq(idP));
    % Normalize data
    Peak = Peak/max(Peak); 
    Width = Width/max(Width);
    Prom = Prom/max(Prom);
    Prox = Prox/max(Prox);
    % Evaluate peak metrics
    factor = Peak +Prom +Width -Prox; [~,idx] = max(factor);
    % Extract harmonic displacement field
    RefAntiFreq = Freq(idP(idx)); ObjFun.RefAntiFreq1(it) = RefAntiFreq;
    AntiHarmDisp1 = abs( HarmDisp1(:,:,idP(idx)) ); % Displacement field
    fprintf('Antiresonance peak #1 identified.\n')
end

%% Modal Assurance Criteria (MAC) for antiresonance f_Ax
if ini.w5 ~= 0
    U_harm = AntiHarmDisp1(:);
    MAC = zeros(length(AntiEigF1),1);
    for idx = 1:size(AntiEigV1,3)
        AEigV = abs( AntiEigV1(:,:,idx) );
        MAC(idx) = ( (U_harm'*AEigV(:))^2 )/...
                   ( (U_harm'*U_harm)*(AEigV(:)'*AEigV(:)) );
    end
    % MAC = MAC./( abs(AntiEigF - RefAntiFreq)/1e3 );
    [~,idMAC] = max(MAC);
    AntiEigF1 = AntiEigF1(idMAC); ObjFun.AntiEigF1(it) = AntiEigF1;
    AntiEigV1 = AntiEigV1(:,:,idMAC); % ObjFun.AntiEigV(:,it) = AntiEigV(:);
    fprintf('MAC for AntiRes1 completed.\n')
    fprintf('Antiresonance Eigefrequencies #1 processed.\n\n')
end

%% Modal FEA-2 for antiresonance eigenfrequency f_Ay
% Run ABAQUS simulation for antiresonance eigenfrequencies
if ini.w6 ~= 0
Inp_AntiRes2 = ['AntiEigFY_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_AntiRes2 ' cpus=20 interactive ask_delete=OFF']);
if SysStat == 0; fprintf('AntiresonanceEigenFreq2 simulation done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
tic; 
[AntiEigV2,~,~,AntiEigF2] = ReadFILResults(Inp_AntiRes2,ini.recordKeysEig);
AntiEigF2 = AntiEigF2( 1:size(AntiEigV2,3) ); % Cut EigF if longer than EigV
AntiEigV2(:,:,AntiEigF2<1) = []; % Delete rigid body modes (1Hz threshold)
AntiEigF2(AntiEigF2<1) = [];     % Delete rigid body modes (1Hz threshold)
fprintf('Results file loaded.\n'); toc

% Remove localized modes
idx_Localized = LocalizedModes(ini,FEM,rho_h(:,it),AntiEigV2);
AntiEigF2(idx_Localized) = [];
AntiEigV2(:,:,idx_Localized) = [];
fprintf('Localized modes removed: %d.\n', length(idx_Localized))

else
    AntiEigF2 = [];
    AntiEigV2 = [];
end

%% Antiresonance peak identification for f_Ay
if ini.w6 ~= 0
    [Peak,idP,Width,Prom] = findpeaks(1./Uavg2);
    Prox = abs(ini.Target2 - Freq(idP));
    % Normalize data
    Peak = Peak/max(Peak); 
    Width = Width/max(Width);
    Prom = Prom/max(Prom);
    Prox = Prox/max(Prox);
    % Evaluate peak metrics
    factor = Peak +Prom +Width -Prox; [~,idx] = max(factor);
    % Extract harmonic displacement field
    RefAntiFreq = Freq(idP(idx)); ObjFun.RefAntiFreq2(it) = RefAntiFreq;
    AntiHarmDisp2 = abs( HarmDisp2(:,:,idP(idx)) ); % Displacement field
    fprintf('Antiresonance peak #2 identified.\n')
end

%% Modal Assurance Criteria (MAC) for antiresonance f_Ay
if ini.w6 ~= 0
    U_harm = AntiHarmDisp2(:);
    MAC = zeros(length(AntiEigF2),1);
    for idx = 1:size(AntiEigV2,3)
        AEigV = abs( AntiEigV2(:,:,idx) );
        MAC(idx) = ( (U_harm'*AEigV(:))^2 )/...
                   ( (U_harm'*U_harm)*(AEigV(:)'*AEigV(:)) );
    end
    % MAC = MAC./( abs(AntiEigF - RefAntiFreq)/1e3 );
    [~,idMAC] = max(MAC);
    AntiEigF2 = AntiEigF2(idMAC); ObjFun.AntiEigF2(it) = AntiEigF2;
    AntiEigV2 = AntiEigV2(:,:,idMAC); % ObjFun.AntiEigV(:,it) = AntiEigV(:);
    fprintf('MAC for AntiRes2 completed.\n')
    fprintf('Antiresonance Eigefrequencies #2 processed.\n\n')
end

%% Modal FEA-3 for antiresonance eigenfrequency f_Az
% Run ABAQUS simulation for antiresonance eigenfrequencies
if ini.w7 ~= 0
Inp_AntiRes3 = ['AntiEigFZ_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_AntiRes3 ' cpus=20 interactive ask_delete=OFF']);
if SysStat == 0; fprintf('AntiresonanceEigenFreq3 simulation done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
tic; 
[AntiEigV3,~,~,AntiEigF3] = ReadFILResults(Inp_AntiRes3,ini.recordKeysEig);
AntiEigF3 = AntiEigF3( 1:size(AntiEigV3,3) ); % Cut EigF if longer than EigV
AntiEigV3(:,:,AntiEigF3<1) = []; % Delete rigid body modes (1Hz threshold)
AntiEigF3(AntiEigF3<1) = [];     % Delete rigid body modes (1Hz threshold)
fprintf('Results file loaded.\n'); toc

% Remove localized modes
idx_Localized = LocalizedModes(ini,FEM,rho_h(:,it),AntiEigV3);
AntiEigF3(idx_Localized) = [];
AntiEigV3(:,:,idx_Localized) = [];
fprintf('Localized modes removed: %d.\n', length(idx_Localized))

else
    AntiEigF3 = [];
    AntiEigV3 = [];
end

%% Antiresonance peak identification for f_Az
if ini.w7 ~= 0
    [Peak,idP,Width,Prom] = findpeaks(1./Uavg3);
    Prox = abs(ini.Target3 - Freq(idP));
    % Normalize data
    Peak = Peak/max(Peak); 
    Width = Width/max(Width);
    Prom = Prom/max(Prom);
    Prox = Prox/max(Prox);
    % Evaluate peak metrics
    factor = Peak +Prom +Width -Prox; [~,idx] = max(factor);
    % Extract harmonic displacement field
    RefAntiFreq = Freq(idP(idx)); ObjFun.RefAntiFreq3(it) = RefAntiFreq;
    AntiHarmDisp3 = abs( HarmDisp3(:,:,idP(idx)) ); % Displacement field
    fprintf('Antiresonance peak #3 identified.\n')
end

%% Modal Assurance Criteria (MAC) for antiresonance f_Az
if ini.w7 ~= 0
    U_harm = AntiHarmDisp3(:);
    MAC = zeros(length(AntiEigF3),1);
    for idx = 1:size(AntiEigV3,3)
        AEigV = abs( AntiEigV3(:,:,idx) );
        MAC(idx) = ( (U_harm'*AEigV(:))^2 )/...
                   ( (U_harm'*U_harm)*(AEigV(:)'*AEigV(:)) );
    end
    % MAC = MAC./( abs(AntiEigF - RefAntiFreq)/1e3 );
    [~,idMAC] = max(MAC);
    AntiEigF3 = AntiEigF3(idMAC); ObjFun.AntiEigF3(it) = AntiEigF3;
    AntiEigV3 = AntiEigV3(:,:,idMAC); % ObjFun.AntiEigV(:,it) = AntiEigV(:);
    fprintf('MAC for AntiRes3 completed.\n')
    fprintf('Antiresonance Eigefrequencies #3 processed.\n\n')
end

%% Modal FEA for resonance eigenfrequencies
% Run ABAQUS simulation for resonance eigenfrequencies
if ini.w8 ~= 0 || ini.w9 ~= 0
Inp_Modal = ['ResEigF_' ini.ID]; tic
SysStat = ...
system(['abaqus job=' Inp_Modal ' cpus=20 interactive ask_delete=OFF']);
if SysStat == 0; fprintf('\nEigenfrequencies simulation done.\n')
else; error('Abaqus Job finished with errors'); end; toc

% Read Results file (.fil)
tic; [ResEigV,~,~,ResEigF] = ReadFILResults(Inp_Modal,ini.recordKeysEig);
if any(ResEigF)
    ResEigF = ResEigF( 1:size(ResEigV,3) ); % Cut EigF if longer than EigV
    ResEigV(:,:,ResEigF<1) = []; % Delete rigid body modes (1Hz threshold)
    ResEigF(ResEigF<1) = [];     % Delete rigid body modes (1Hz threshold)
    fprintf('Results file loaded.\n'); toc
    % Remove localized modes
    idx_Localized = LocalizedModes(ini,FEM,rho_h(:,it),ResEigV);
    ResEigF(idx_Localized) = [];
    ResEigV(:,:,idx_Localized) = [];
    fprintf('Localized modes removed: %d.\n', length(idx_Localized))
    % ObjFun.ResEigF(1:length(ResEigF),it) = ResEigF; % Store data
    fprintf('Resonance Eigemodes processed.\n\n')
else
    fprintf('No resonance frequencies found.\n')
end

end

%% Resonance peak identification
if ini.w8 ~= 0
    [~,idP,~,~] = findpeaks(Uavg1);
    if any(idP);HarmField1=abs(HarmDisp1(:,:,idP)); else; HarmField1=[];end
    fprintf('Resonance peaks #1 identified.\n')
end

if ini.w9 ~= 0
    [~,idP,~,~] = findpeaks(Uavg3);
    if any(idP);HarmField3=abs(HarmDisp3(:,:,idP)); else; HarmField3=[];end
    fprintf('Resonance peaks #3 identified.\n')
end

%% Modal Assurance Criteria (MAC) for resonances R1
if ini.w8 ~= 0
if ~isempty(HarmField1)
    EigF1 = zeros(size(HarmField1,3),1);
    EigV1 = zeros(size(ResEigV,1),size(ResEigV,2),size(HarmField1,3));
    for Freq_idx = 1:size(HarmField1,3)
        U_harm = HarmField1(:,:,Freq_idx);
        MAC = zeros(length(ResEigF),1);
        for idx = 1:size(ResEigV,3)
            RefEigV = abs( ResEigV(:,:,idx) );
            MAC(idx) = ( (U_harm(:)'*RefEigV(:))^2 )/...
                       ( (U_harm(:)'*U_harm(:))*(RefEigV(:)'*RefEigV(:)) );
        end
        [~,idMAC] = max(MAC);
        EigF1(Freq_idx) = ResEigF(idMAC);
        ObjFun.ResEigF1(it,Freq_idx) = EigF1(Freq_idx);
        EigV1(:,:,Freq_idx) = ResEigV(:,:,idMAC);
    end
else
    EigF1 = [];
    EigV1 = [];
end

% % Bypass MAC for resonances, so ALL EigF are selected
% EigF1 = ResEigF;
% EigV1 = ResEigV;
% ObjFun.ResEigF1(it,1:length(ResEigF)) = EigF1; % save frequencies

else
    EigF1 = [];
    EigV1 = [];
end

%% Modal Assurance Criteria (MAC) for resonances R3
if ini.w9 ~= 0
if ~isempty(HarmField3)
    EigF3 = zeros(size(HarmField3,3),1);
    EigV3 = zeros(size(ResEigV,1),size(ResEigV,2),size(HarmField3,3));
    for Freq_idx = 1:size(HarmField3,3)
        U_harm = HarmField3(:,:,Freq_idx);
        MAC = zeros(length(ResEigF),1);
        for idx = 1:size(ResEigV,3)
            RefEigV = abs( ResEigV(:,:,idx) );
            MAC(idx) = ( (U_harm(:)'*RefEigV(:))^2 )/...
                       ( (U_harm(:)'*U_harm(:))*(RefEigV(:)'*RefEigV(:)) );
        end
        [~,idMAC] = max(MAC);
        EigF3(Freq_idx) = ResEigF(idMAC);
        ObjFun.ResEigF3(it,Freq_idx) = EigF3(Freq_idx);
        EigV3(:,:,Freq_idx) = ResEigV(:,:,idMAC);
    end
else
    EigF3 = [];
    EigV3 = [];
end

% % Bypass MAC for resonances, so ALL EigF are selected
% EigF3 = ResEigF;
% EigV3 = ResEigV;
% ObjFun.ResEigF3(it,1:length(ResEigF)) = EigF3; % save frequencies

else
    EigF3 = [];
    EigV3 = [];
end