%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Developed by Daniel Giraldo Guzman %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Penn State University %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Contact: dzg5526@psu.edu %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% danielgiraldoguzman@gmail.com %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; format compact; format shortG; clc
fprintf('Topology Optimization of Resonant Structures\n\n')

%% Profile execution
% profile on
% Make sure to save Profiler results at the end of the code

%% Initialization
ini = InitialParameters(1);

switch ini.start
case 'Initial'
    % Generate the base model. Initial Abaqus FEA simulation
    Inp_file = ['Job-ini_' ini.ID]; % Input file name
    SysStat = ...
    system(['abaqus job=' Inp_file ' cpus=20 interactive ask_delete=OFF']);
    if SysStat == 0; fprintf('\nInitial FE simulation done.\n')
    else; error('Abaqus Job finished with errors'); end     
    
    % Process FEA Matrix data
    Mesh;    LoadFEAdata;    NegativeFEAmatrices
    clearvars -EXCEPT ini FEM % Clear data
    % load('FEMdata'); fprintf('Initial data retrieved \n')
    
    % Pre-allocation
    rho = zeros(FEM.mesh.nel,ini.MaxIter);   % Design variables
    % rho_d = zeros(FEM.mesh.nel,ini.MaxIter); % Density-filtered density
    rho_h = zeros(FEM.mesh.nel,ini.MaxIter); % Heaviside-filtered density
    Freq = [];
    
    ObjFun.ObjFunVal = zeros(ini.MaxIter,1); % Objective Function value
    ObjFun.Szz = zeros(ini.MaxIter,1);       % Stress values
    ObjFun.Ux = zeros(ini.MaxIter,1);        % Displacement values
    ObjFun.Uy = zeros(ini.MaxIter,1);        % Displacement values
    ObjFun.Uz = zeros(ini.MaxIter,1);        % Displacement values
    ObjFun.FRF1 = zeros(1,ini.MaxIter); % Frequency Response Function 1
    ObjFun.FRF2 = zeros(1,ini.MaxIter); % Frequency Response Function 2
    ObjFun.FRF3 = zeros(1,ini.MaxIter); % Frequency Response Function 3
    ObjFun.AntiEigF1 = zeros(ini.MaxIter,1); % Antiresonance Eigenfreq1
    ObjFun.AntiEigF2 = zeros(ini.MaxIter,1); % Antiresonance Eigenfreq2
    ObjFun.AntiEigF3 = zeros(ini.MaxIter,1); % Antiresonance Eigenfreq3
    ObjFun.RefAntiFreq1 = zeros(ini.MaxIter,1); % AntiRes1 from FRF1
    ObjFun.RefAntiFreq2 = zeros(ini.MaxIter,1); % AntiRes2 from FRF2
    ObjFun.RefAntiFreq3 = zeros(ini.MaxIter,1); % AntiRes2 from FRF3
    ObjFun.ResEigF1 = zeros(ini.MaxIter,1); % Resonance Eigenfrequencies1
    ObjFun.ResEigF3 = zeros(ini.MaxIter,1); % Resonance Eigenfrequencies3    
    if ini.VerifySensitivities == 1
        ObjFun.dFun = zeros(FEM.mesh.nel,ini.MaxIter);
    end
    
    % Initialization
    ini = InitialParameters(2,ini,FEM);
    rho(:,1) = ini.rho;
    ini.Iter = 1;
    ini.beta = ini.beta_ini;
    fprintf('Optimization to start at iteration: %d \n', ini.Iter)
    save('FEMdata','FEM','-v7.3')
    fprintf('Initialization completed.\n\n')
    
case 'Continue'
    load('DataToContinue')
    load('FEMdata')
    ini.Iter = find(ObjFun.ObjFunVal == 0, 1, 'first');
    fprintf('\nData to continue loaded \n')
    fprintf('Optimization to start at iteration: %d \n', ini.Iter)
    
otherwise
    error('Invalid starting point. It must be "Initial" or "Continue"')
end

%% OPTIMIZATION START
fprintf('Optimization Begins ------------------------------------------\n')
Vr = [ini.MaxVol*FEM.TotalVolume; -ini.MinVol*FEM.TotalVolume];
Ve = [FEM.mesh.Ve'; -FEM.mesh.Ve'];
display(ini)
% figure('position', [0, 0, 1600, 700])

for it = ini.Iter:ini.MaxIter
tStart = tic;
fprintf('\n\n\nIteration: %d \n',it)
fprintf('--------------------------------------------------------------\n')

%% Filter
if mod(it,ini.beta_it) == 1 % Increase beta every X iterations
    ini.beta = ini.beta + ini.beta_step; end
[~,rho_h(:,it),drho_d,drho_h] = ProjectionFilter(ini,FEM,rho(:,it));

%% Analysis module
FEanalysis

%% Objective function
% BC terms
if ini.w1 ~= 0; ObjFun.Szz(it) = norm(FEM.Szz(:),2); end
if ini.w2 ~= 0; ObjFun.Ux(it) = norm(FEM.Ux(:),2); end
if ini.w3 ~= 0; ObjFun.Uy(it) = norm(FEM.Uy(:),2); end
if ini.w4 ~= 0; ObjFun.Uz(it) = norm(FEM.Uz(:),2); end

% Antiresonance terms
if any(AntiEigF1); AntiTerm1 = ( (AntiEigF1 -ini.Target1)/ini.Target1 )^2;
else; AntiTerm1 = 0; end

if any(AntiEigF2); AntiTerm2 = ( (AntiEigF2 -ini.Target2)/ini.Target2 )^2;
else; AntiTerm2 = 0; end

if any(AntiEigF3); AntiTerm3 = ( (AntiEigF3 -ini.Target3)/ini.Target3 )^2;
else; AntiTerm3 = 0; end

% Resonance terms
if any(EigF1)
    if any(AntiEigF1); EigTerm1 = sqrt((AntiEigF1./(EigF1-AntiEigF1)).^2);
    else; EigTerm1 = sqrt((ini.Target1./(EigF1-ini.Target1)).^2); end
else
    EigTerm1 = 0;
end

if any(EigF3)
    if any(AntiEigF3); EigTerm3 = sqrt((AntiEigF3./(EigF3-AntiEigF3)).^2);
    else; EigTerm3 = sqrt((ini.Target3./(EigF3-ini.Target3)).^2); end
else
    EigTerm3 = 0;
end

% Objective function
ObjFun.ObjFunVal(it) = ini.w1* ObjFun.Szz(it) + ...
                       ini.w2* ObjFun.Ux(it) + ...
                       ini.w3* ObjFun.Uy(it) + ...
                       ini.w4* ObjFun.Uz(it) + ...
                       ini.w5* AntiTerm1 + ...
                       ini.w6* AntiTerm2 + ...
                       ini.w7* AntiTerm3 + ...
                       ini.w8* sum(EigTerm1) + ...
                       ini.w9* sum(EigTerm3);
fprintf('Objective function: %.5g \n',ObjFun.ObjFunVal(it))
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%% SHOW LAST 3 ITERATIONS %%%%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% Sensitivity analysis
tic; fprintf('\nSentivity analysis started.\n')
dObj = Sensitivity(ini,FEM,U_stX,U_stY,U_stZ,E_stZ,...
    AntiEigF1,AntiEigF2,AntiEigF3,EigF1,EigF3,...
    AntiEigV1,AntiEigV2,AntiEigV3,EigV1,EigV3,...
    rho_h(:,it),drho_d,drho_h);
fprintf('Sentivity analysis finished.\n'); toc

if ini.VerifySensitivities == 1; tic; ObjFun.dFun(:,it) = dObj; %#ok
    dFDiff(:,it) = FiniteDiff(ini,FEM,rho(:,it),it); %#ok
    save('SensVerify','dFun','dFDiff'); toc; end

% This part has been moved to the Sensitivity function
% % Check for void non-design elements.
% if any(ini.Sets.VoidNonDesign)
%     dObj(ini.Sets.VoidNonDesign) = 0;
%     fprintf('Void non-design elements sensitivity enforced to zero.\n')
% end
% % Check for solid non-design elements
% if any(ini.Sets.SolidNonDesign)
%     dObj(ini.Sets.SolidNonDesign) = 0;
%     fprintf('Solid non-design elements sensitivity enforced to zero.\n')
% end
% ObjFun.dFun(:,it) = dObj; % save the current sensitivity

%% Optimization solver
% [LB, UB] = MoveLimits(ini,FEM,rho(:,it)); % legacy function
LB = rho(:,it) - ini.step;    % Lower Move Limit
UB = rho(:,it) + ini.step;    % Upper Move Limit
[LB, UB] = MoveLimitsCheck(UB,LB,ini,FEM,rho(:,it)); % Check limits

try
    rho(:,it+1) = linprog(dObj,Ve,Vr,[],[],LB,UB,ini.opts);
catch
    try
    [LB, UB] = MoveLimitsCheck(UB,LB,ini,FEM,rho(:,it)); % Check limits
    rho(:,it+1) = linprog(dObj,Ve,Vr,[],[],LB,UB,ini.opts);
    catch
        % Relax move limits and run it again        
        LB = LB/1.5; LB(LB<1e-9) = 1e-9; UB = UB*1.5; UB(UB>1) = 1;
        [LB, UB] = MoveLimitsCheck(UB,LB,ini,FEM,rho(:,it)); % Check limits
        fprintf('Move Limits relaxed.\n')
        rho(:,it+1) = linprog(dObj,Ve,Vr,[],[],LB,UB,ini.opts);
    end
end

% Check for forced elements
if any(ini.Sets.VoidNonDesign); rho(ini.Sets.VoidNonDesign,it+1) = 1e-9;end
if any(ini.Sets.SolidNonDesign); rho(ini.Sets.SolidNonDesign,it+1) = 1; end
fprintf('SLP solution done.\n')

%% Save temporary data
if exist('Freq','var'); save('TopOpt','ini','ObjFun','rho','rho_h','Freq')
else; save('TopOpt','ini','ObjFun','rho','rho_h'); end
fprintf('End of Iteration: %d \n',it)
tEnd = toc(tStart); display(tEnd)

end % end of optimization loop

%% Save final results
if exist('Freq','var'); save('TopOpt','ini','ObjFun','rho','rho_h','Freq')
else; save('TopOpt','ini','ObjFun','rho','rho_h'); end
fprintf('\n\nEnd of Optimization \n\n')

%% Profiler save
% profsave