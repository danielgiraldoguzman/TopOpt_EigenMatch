% Import FEA data

% Import Global matrix data
StepMTX = '2'; % Step in which the matrices are requested to be exported
StifMTXname = [Inp_file '_STIF' StepMTX '.mtx'];
MassMTXname = [Inp_file '_MASS' StepMTX '.mtx'];
DampMTXname = [Inp_file '_DMPV' StepMTX '.mtx'];
LoadMTXname = [Inp_file '_LOAD' StepMTX '.mtx'];
Kmif = load(StifMTXname); % Matrix Input Format
Mmif = load(MassMTXname); % Matrix Input Format
Cmif = load(DampMTXname); % Matrix Input Format
% FEM.K = MIF2sparse(Kmif);
% FEM.M = MIF2sparse(Mmif);
% FEM.C = MIF2sparse(Cmif);

% Create Force vector
LoadMIF = ImportLoad(LoadMTXname);
LoadDOFs = 3*LoadMIF(:,1) + LoadMIF(:,2) - 3; % Calculate DOFs position
% Think about including complex forcing functions
FEM.Force = sparse(FEM.mesh.dof,1); % Allocate force vector
FEM.Force(LoadDOFs) = LoadMIF(:,3);

% *MATRIX EXPORT ELEMENT BY ELEMENT.
StepMTX = '3'; % Step in which the matrices are requested to be exported
StifMTXname = [Inp_file '_STIF' StepMTX '.mtx'];
MassMTXname = [Inp_file '_MASS' StepMTX '.mtx'];
DampMTXname = [Inp_file '_DMPV' StepMTX '.mtx'];
% LoadMTXname = [Inp_file '_LOAD' StepMTX '.mtx'];
k_sp = load(StifMTXname); % Matrix Sparse Format
m_sp = load(MassMTXname); % Matrix Sparse Format
c_sp = load(DampMTXname); % Matrix Sparse Format
% Element by element matrix assembly
DOF = size(FEM.mesh.edof,2); % 3 dof per node
k_all = zeros(DOF*DOF,FEM.mesh.nel);
m_all = zeros(DOF*DOF,FEM.mesh.nel);
c_all = zeros(DOF*DOF,FEM.mesh.nel);
parfor n = 1:FEM.mesh.nel
    idx_k = find(k_sp(:,1)==(n-1)); %#ok
    idx_m = find(m_sp(:,1)==(n-1)); %#ok
    idx_c = find(c_sp(:,1)==(n-1)); %#ok
    k_el = spconvert(k_sp(idx_k,2:4)); k_el = triu(k_el) + tril(k_el',-1);
    m_el = spconvert(m_sp(idx_m,2:4)); m_el = triu(m_el) + tril(m_el',-1);
    c_el = spconvert(c_sp(idx_c,2:4)); c_el = triu(c_el) + tril(c_el',-1);
    
    k_all(:,n) = k_el(:);
    m_all(:,n) = m_el(:);
    c_all(:,n) = c_el(:);
end
FEM.k = k_all;
FEM.m = m_all;
FEM.c = c_all;
fprintf('FEA Matrix data loaded.\n')

% Global Matrix Assembly
% n = FEM.mesh.dof;                 % Total DOF's number
% Ai = FEM.mesh.Ai;
% Aj = FEM.mesh.Aj;
% FEM.K = sparse(Ai,Aj,FEM.k,n,n); % Stiffness global matrix
% FEM.M = sparse(Ai,Aj,FEM.m,n,n); % Mass global matrix
% FEM.C = sparse(Ai,Aj,FEM.c,n,n); % Damping global matrix

% Import sparse matrix data
% StepMTX = '3'; % Step in which the matrices are requested to be exported
% StifMTXname = [Inp_file '_STIF' StepMTX '.mtx'];
% MassMTXname = [Inp_file '_MASS' StepMTX '.mtx'];
% Ksp = load(StifMTXname); % Matrix Sparse Format
% Msp = load(MassMTXname); % Matrix Sparse Format
% Ksp = spconvert(Ksp);
% Msp = spconvert(Msp);

% Save imported data
% save('FEMdata','FEM');