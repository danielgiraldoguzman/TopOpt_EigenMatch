%% Harmonic FEA (Frequency Resposnse Function)
% This FEA run is to check the harmonic response, looking for minima
load('TopOpt21.mat') % load data from optimization
load('FEMdata1372.mat')
it = 100; % choose it to run

% Create FEA matrices
fprintf('Writing updated FE Matrices...  ')
FEM = WriteUpdatedMTX(ini,rho_h(:,it),FEM);
fprintf('Done\n\n')

%% Freq response in x-direction
if ini.w2 ~= 0
    Inp_Harm = ['FreqResX_' ini.ID]; tic
    SysStat = ...
    system(['abaqus job=' Inp_Harm ini.cpus ' interactive ask_delete=OFF']);
    if SysStat == 0; fprintf('FRF done.\n')
    else; error('Abaqus Job finished with errors'); end; toc
    % read results
    fprintf('Reading results file...  ')
    tic; HarmDisp = ReadFILResults(Inp_Harm,ini.recordKeysEig);
    fprintf('Done\n'); toc
    % Extract FRF on the Aniresonance nodes only
    FRF_abs = abs( HarmDisp(ini.AntiResBCnodes,:,:) );
    Uabs = reshape(FRF_abs(:,3-ini.DOFkey1,:),length(ini.AntiResBCnodes),[]);
    FRF1 = mean(Uabs,1)';
    fprintf('FRF for z-axis force processed.\n\n')
else
    FRF1 = 0;
end

% Extract frequency vector
if isempty(Freq)
    Freq = ReadFreq(Inp_Harm);
    fprintf('ODB file processed.\n');
end

% save results
save('FEAonly','FRF1','Freq')

%% Freq response in z-direction
if ini.w4 ~= 0
    Inp_Harm = ['FreqResZ_' ini.ID]; tic
    SysStat = ...
    system(['abaqus job=' Inp_Harm ini.cpus ' interactive ask_delete=OFF']);
    if SysStat == 0; fprintf('FRF done.\n')
    else; error('Abaqus Job finished with errors'); end; toc
    % read results
    fprintf('Reading results file...  ')
    tic; HarmDisp = ReadFILResults(Inp_Harm,ini.recordKeysEig);
    fprintf('Done\n'); toc
    % Extract FRF on the Aniresonance nodes only
    FRF_abs = abs( HarmDisp(ini.AntiResBCnodes,:,:) );
    Uabs = reshape(FRF_abs(:,3-ini.DOFkey3,:),length(ini.AntiResBCnodes),[]);
    FRF3 = mean(Uabs,1)';
    fprintf('FRF for z-axis force processed.\n\n')
else
    FRF3 = 0;
end

% Extract frequency vector
if isempty(Freq)
    Freq = ReadFreq(Inp_Harm);
    fprintf('ODB file processed.\n');
end

% save results
save('FEAonly','FRF1','FRF3','Freq')