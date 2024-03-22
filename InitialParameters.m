function ini = InitialParameters(step,ini,FEM)

%% Run different step accordingly
switch step
case 1
ini.ID = '1372'; % ID is the number of elements and frequency
ini.cpus = ' cpus=40'; % Number of cores to be used in Abaqus
parpool('Processes',[1 100],'IdleTimeout', 120) % Start a parallel pool

%% General initial parameters
ini.start = 'Initial';  % 'Initial' or 'Continue'
ini.MaxIter = 200;      % Maximum number of iterations
ini.MatModel = 'SIMP';  % Material Model definition: "SIMP" or "RAMP"
ini.p = 1;              % Penalization factor
% ini.FreqBC1 = 30e3;     % Frequency of BC approach
% ini.FreqBC2 = 00e3;     % Frequency of BC approach
% ini.FreqBC3 = 40e3;     % Frequency of BC approach
ini.Target1 = 30e3;     % Target frequency #1
ini.Target2 = 30e3;     % Target frequency #2
ini.Target3 = 30e3;     % Target frequency #3
ini.MaxVol = 0.9;       % Max volume
ini.MinVol = 0.1;       % Min volume
ini.w1 = 1e-3;  % Stress-zz BC
ini.w2 = 1e9;   % Displacement BC x-disp
ini.w3 = 1e9;   % Displacement BC y-disp
ini.w4 = 0e9;   % Displacement BC z-disp
ini.w5 = 0e0;   % Antiresonance x-axis
ini.w6 = 0e0;   % Antiresonance y-axis
ini.w7 = 0e0;   % Antiresonance z-axis
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Pending: add y-axis ResFreq. make weights to start at w0
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ini.w8 = 0;     % Resonance gap #1 for AntiRes in x-axis
ini.w9 = 0;     % Resonance gap #3 for AntiRes in z-axis
ini.ThresholdUP = 0.7; % Upper threshold for localized modes (solid)
ini.ThreshoDOWN = 0.1; % Lower threshold for localized modes (soft)
ini.SymmBC = 'SymmY'; % 'SymmY', 'SymXY', 'none'
ini.AntiResBC = 'YesXYZ';
ini.DOFkey1 = 2; % 2--> x
ini.DOFkey2 = 1; % 1--> y
ini.DOFkey3 = 0; % 0--> z

ini.recordKeysBC = [101 11 21]; % 101->U, 11->S, 21->E,
ini.recordKeysEig = [101 1980]; % 101->U, 1980->EigF
ini.Sensitivity = 'Analytical'; % 'Analytical','Numerical','Verify'
ini.SolidNonDesign = 'OutElements';
ini.VoidNonDesign = 'FirstLayer';
ini.SolidElements = 'none';
ini.VoidElements = 'none';
ini.ElementType = 20;   % 8-node or 20-node elements
ini.VerifySensitivities = 0; % Verify sensitivies, 1 to activate

%% Filters
ini.Rdf = 0.004; % Filter radius for Heaviside and Sensitivity filter
ini.beta_step = 2; % Value to increase beta every # iterations
ini.beta_ini = 30 - ini.beta_step; % Starting beta value
ini.beta_it = 2; % # of iterations to increase Beta

%% Sequential linear programming (SLP) solver options
ini.step = 0.01; % move limits step size
ini.opts = optimset;
% ini.opts = optimset(ini.opts,'Algorithm','interior-point-legacy');
% ini.opts = optimset(ini.opts,'Algorithm','interior-point');
ini.opts = optimset(ini.opts,'Algorithm','dual-simplex');
ini.opts = optimset(ini.opts,'Display','iter');
% ini.opts = optimset(ini.opts,'LargeScale','on');
ini.opts = optimset(ini.opts,'AlwaysHonorConstraints','bounds');
ini.opts = optimset(ini.opts,'TolCon',1e-9);
% ini.opts = optimset(ini.opts,'TolFun',1e-6);

%% Material Properties
ini.MatProperties = 'ABS-like';
switch ini.MatProperties
    case 'ABS-like'
        Y = 3.9e9;  % Young's modulus
        v = 0.33;   % Poisson's ratio
        % Densitiy = 1214.650481 [kg/m3]
    case 'Aluminum'
        Y = 69e9;   % Young's modulus
        v = 0.33;   % Poisson's ratio
        % Densitiy = 2700 [kg/m3]
end

% Lame Constants
lambda = (Y*v)/((1 + v)*(1 - 2*v));
mu = Y/(2*(1+v));

% Isotropic Stiffness Matrix
cE = [lambda+2*mu lambda      lambda      0   0   0;
      lambda      lambda+2*mu lambda      0   0   0;
      lambda      lambda      lambda+2*mu 0   0   0;
      0           0           0           mu  0   0;
      0           0           0           0   mu  0;
      0           0           0           0   0   mu];

ini.cE = cE; % use for stress sensitivity
%
%
%
%
%
%
%
%
%
% End of Initial step 1

case 2
% Initial point
ini.rho = 0.5*ones(FEM.mesh.nel,1); % rho must be written column-wise

%% Output set to separate information (from inp file!!!)
ini.Sets.TotalNodes = 1:FEM.mesh.nnod;     % Set containing all nodes
ini.Sets.TotalElements = 1:FEM.mesh.nel;   % Set containing all elements

switch ini.ID
case '8' % Iter time: 1 sec (old)
    ini.Sets.OutNodes = [7:9,16:18,25:27,49,52,54:56,59,60,75,77:79,81];
    ini.Sets.OutElements = [3, 4, 7, 8];

case '64' % Iter time: 20 sec (old)
    ini.Sets.OutNodes = [47:49, 72:74, 97:99, 216, 221, 285, 287, 289,...
        290, 292, 350, 352, 354, 355, 357];
    ini.Sets.OutElements = [30, 31, 46, 47];

case '500' % Iter time: from ~2 hours, down to ~2 minutes
    ini.Sets.OutNodes = [ 325,  326,  391,  392,  457,  458, ...
                         1550, 1731, 1733, 1734, 1912, 1914, 1915];
    ini.Sets.OutElements = [246, 296];

case '1372'
    ini.Sets.OutNodes = [713, 714, 715, 833, 834, 835, 953, 954, 955, ...
        1073, 1074, 1075, 1193, 1194, 1195, 3681, 3685, 4018, 4020,...
        4021, 4022, 4024, 4355, 4357, 4358, 4359, 4361, 4692, 4694,...
        4695, 4696, 4698, 5029, 5031, 5032, 5033, 5035];
    ini.Sets.OutElements = [582, 583, 680, 681, 778, 779, 876, 877];

case '4000'
    ini.Sets.OutNodes = [1838, 1839, 1840, 1841, 2069, 2070, 2071, 2072,...
        2300, 2301, 2302, 2303, 2531, 2532, 2533, 2534, 2762, 2763, 2764,...
        2765, 2993, 2994, 2995, 2996, 3224, 3225, 3226, 3227, 9878, 9882,...
        9885, 10539, 10541, 10542, 10543, 10545, 10546, 10548, 11200,...
        11202, 11203, 11204, 11206, 11207, 11209, 11861, 11863, 11864,...
        11865, 11867, 11868, 11870, 12522, 12524, 12525, 12526, 12528,...
        12528, 12531, 13183, 13185, 13186, 13187, 13189, 13190, 13192,...
        13844, 13846, 13847, 13848, 13850, 13851, 13853];
    ini.Sets.OutElements = [ 1591, 1592, 1593, 1791, 1792, 1793, ...
                             1991, 1992, 1993, 2191, 2192, 2193, ...
                             2391, 2392, 2393, 2591, 2592, 2593];
otherwise
    error('ID case no available to define initial sets')
end
ini.rho(ini.Sets.OutElements) = 1; % Set out elements to solid

%% Non-design elements
BottomElements = find( FEM.mesh.Ce(:,3) <= 1.02*min(FEM.mesh.Ce(:,3)) )';
TopElements = find( FEM.mesh.Ce(:,3) >= 0.98*max(FEM.mesh.Ce(:,3)) );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Make the nondesign elements compatible with Set.OutElements
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
switch ini.SolidNonDesign
    case 'Bottom'
        ini.Sets.SolidNonDesign = BottomElements;
        ini.rho(ini.Sets.SolidNonDesign) = 1; % Set elements to solid 
    case 'OutElements'
        ini.Sets.SolidNonDesign = ini.Sets.OutElements;
        ini.rho(ini.Sets.SolidNonDesign) = 1; % Set elements to solid
    case 'none'
        ini.Sets.SolidNonDesign = [];
    otherwise
        ini.Sets.SolidNonDesign = [];
        warning('Incorrect definition of solid non-design elements')
end

switch ini.VoidNonDesign
    case 'FirstLayer'
        ini.Sets.VoidNonDesign = BottomElements;
        ini.Sets.VoidNonDesign = ...
            setdiff(ini.Sets.VoidNonDesign,ini.Sets.SolidNonDesign);
        ini.rho(ini.Sets.VoidNonDesign) = 1e-6; % Set elements to void   
    case 'none'
        ini.Sets.VoidNonDesign = [];
    otherwise
        ini.Sets.VoidNonDesign = [];
        warning('Incorrect definition of void non-design elements')
end

%% Starting point
switch ini.VoidElements
    case'First1mm'
        idx_1mmElements = FEM.mesh.Ce(:,3) <= 1e-3;
        ini.rho(idx_1mmElements) = 1e-6;
        fprintf('First 1mm layers set to void\n');
end

switch ini.SolidElements
    case 'Bottom'
        ini.rho(BottomElements) = 1;
    case 'Top'
        ini.rho(TopElements) = 1;
    case 'Bottom_Top'
        ini.rho(BottomElements) = 1;
        ini.rho(TopElements) = 1;
    case {'CenterStem','5mmCenterStem','11mmCenterStem'}
        % coordinates for center elements at the bottom
        xcoord = FEM.mesh.Ce([ini.CenterElements],1);
        ycoord = FEM.mesh.Ce([ini.CenterElements],2);
        % find elements with the same coordinates
        ElementsX = FEM.mesh.Ce(:,1) == xcoord';
        ElementsY = FEM.mesh.Ce(:,2) == ycoord';
        StemElements = and(ElementsX,ElementsY);
        Stem = sum(StemElements,2);
        switch ini.SolidElements
            case 'CenterStem'
                ini.rho(Stem==1) = 1;
                fprintf('Center stem added to starting point\n');
            case '5mmCenterStem'
                ElementsZ = FEM.mesh.Ce(:,3) <= 5e-3; % [mm]
                HalfStem = and(Stem,ElementsZ);
                ini.rho(HalfStem==1) = 1;
                fprintf('Half center stem added to starting point\n');
            case '11mmCenterStem'
                ElementsZ = FEM.mesh.Ce(:,3) <= 11e-3; % [mm]
                HalfStem = and(Stem,ElementsZ);
                ini.rho(HalfStem==1) = 1;
                fprintf('Half center stem added to starting point\n');
        end
end

if ini.VerifySensitivities == 1; ini.MaxVol = 1;
   ini.rho = 0.95*ini.MaxVol*ones(FEM.mesh.nel,1); end

%% Volume control
VolumeFraction = ( (FEM.mesh.Ve')*ini.rho )/sum(FEM.mesh.Ve);
VolCount = 0;
while (VolumeFraction>ini.MaxVol || VolumeFraction<ini.MinVol) && VolCount<50
    eta = VolumePreserving(ini.rho,2,FEM.mesh.Ve,ini.MaxVol,ini.MinVol);
    ini.rho = HeavisideFilter(2,eta,ini.rho);
    % Check for forced elements
    if any(ini.Sets.VoidNonDesign); ini.rho(ini.Sets.VoidNonDesign) = 1e-6; end
    if any(ini.Sets.SolidNonDesign); ini.rho(ini.Sets.SolidNonDesign) = 1; end
    VolumeFraction = ( (FEM.mesh.Ve')*ini.rho )/sum(FEM.mesh.Ve);
    VolCount = VolCount + 1;
end
VolumeFraction = ( (FEM.mesh.Ve')*ini.rho )/sum(FEM.mesh.Ve);
fprintf('Initial Volume Fraction: %.4g \n',VolumeFraction)
ini.rho(ini.rho<1e-6) = 1e-6; ini.rho(ini.rho>1) = 1;

%% Boundary conditions
% make sure the BCnodes file is filled up with zeros

ini.Sets.BCdof = []; % for BC matchign aproach
% If antiresonance eigenfrequencies are to be computed:
switch ini.AntiResBC
    case 'YesXYZ'
        ini.AntiResBCnodes = load(['AntiResBCnodes' ini.ID '.txt'])';
        ini.AntiResBCnodes = ini.AntiResBCnodes(:);
        ini.AntiResBCnodes(ini.AntiResBCnodes==0) = [];
        ini.Sets.AntiRes1BCdof = ini.AntiResBCnodes*3 - ini.DOFkey1;
        ini.Sets.AntiRes2BCdof = ini.AntiResBCnodes*3 - ini.DOFkey2;
        ini.Sets.AntiRes3BCdof = ini.AntiResBCnodes*3 - ini.DOFkey3;
    otherwise
        error('Incorrect definition of Boundary conditions')
end

% Symmetry conditions
switch ini.SymmBC
    case 'SymXY'
    % fill up with zeros manully
    SymmXBCnodes = load(['SymmXBCnodes' ini.ID '.txt'])';
    SymmXBCnodes = SymmXBCnodes(:); SymmXBCnodes(SymmXBCnodes==0) = [];
    SymmYBCnodes = load(['SymmYBCnodes' ini.ID '.txt'])';
    SymmYBCnodes = SymmYBCnodes(:); SymmYBCnodes(SymmYBCnodes==0) = [];
    ini.Sets.SymmBCdof = [SymmXBCnodes*3 - 2; SymmYBCnodes*3 - 1];

    case 'SymmY'
    % fill up with zeros manully
    SymmYBCnodes = load(['SymmYBCnodes' ini.ID '.txt'])';
    SymmYBCnodes = SymmYBCnodes(:); SymmYBCnodes(SymmYBCnodes==0) = [];
    ini.Sets.SymmBCdof = SymmYBCnodes*3 - 1;

    case 'none'
    ini.Sets.SymmBCdof = [];
end

% End of Initial step 2
end

% End of InitialParameters function
end