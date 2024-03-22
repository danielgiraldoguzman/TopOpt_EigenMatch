% Generate mesh (assuming a structured mesh with hexahedral elements)

%% Mesh generation
eType = ini.ElementType; % Pre-allocation
switch eType
    case 8
        enodes = load(['enodes' ini.ID '.txt']);
        enodes(:,1) = []; % Remove element numbers
    case 20
        enodes = enodes20(['enodes' ini.ID '.txt']);
        enodes = reshape(enodes',size(enodes,2)*2,size(enodes,1)/2)';
        enodes(:,~any(enodes,1)) = []; % Remove blank spaces
        enodes(:,1) = []; % Remove element numbers
    otherwise
        error('Wrong ini.ElementType definition')
end

FEM.mesh.nel = size(enodes,1);
DOF = size(enodes,2)*3;
edof = zeros(DOF,FEM.mesh.nel);
Ai = zeros(DOF*DOF,FEM.mesh.nel);
Aj = zeros(DOF*DOF,FEM.mesh.nel);
parfor n = 1:FEM.mesh.nel
    edof(:,n) = medof(enodes(n,:),eType); % DOFs arrange
    u = repmat(edof(:,n),1,DOF); % DOFs row index to matrix assembly
    v = u';                      % DOFs columns index to matrix assembly
    Ai(:,n) = u(:);              % Row index for sparse matrix creation
    Aj(:,n) = v(:);              % Columns index for matrix creation
end

FEM.mesh.conect = enodes;         % Store connectivity mtrix
FEM.mesh.edof = edof';            % Store element dof
FEM.mesh.nnod = max(max(enodes)); % Total number of nodes
FEM.mesh.dof = 3*FEM.mesh.nnod;   % Total Number of DOFs
FEM.mesh.Ai = Ai;
FEM.mesh.Aj = Aj;

FEM.mesh.ncoord = load(['ncoord' ini.ID '.txt']);% Import nodal coordinates
FEM.mesh.ncoord(:,1) = []; % Remove node numbers

%% Element volume and centroids
Ce = zeros(FEM.mesh.nel,3);
Ve = zeros(FEM.mesh.nel,1);
ncoord = FEM.mesh.ncoord;
conect = FEM.mesh.conect;
parfor i = 1:FEM.mesh.nel
    n = ncoord(conect(i,1:8),:); %#ok
    element_center = (max(n) - min(n))/2;   % Relative center
    Ce(i,:) = min(n) + element_center;      % Absolute center
    Ve(i,1) = prod(max(n) - min(n));        % Element Volume
end

% Distance between centroids of element i and element j
Rij = zeros(FEM.mesh.nel,FEM.mesh.nel);
nel = FEM.mesh.nel;
parfor i = 1:nel
    for j = 1:nel    
        if i ~= j     
            Rij(i,j) = sqrt((Ce(i,1)-Ce(j,1))^2 +...
                            (Ce(i,2)-Ce(j,2))^2 +...
                            (Ce(i,3)-Ce(j,3))^2); %#ok
        end
    end
end

% Store data
FEM.mesh.Ce = Ce;
FEM.mesh.Ve = Ve;
FEM.mesh.Rei = Rij;
FEM.TotalVolume  = sum(FEM.mesh.Ve);    % Total volume
fprintf('Mesh created.\n')

function edof = medof(enodes,dof)
    DOFs = [3*enodes(1:dof) - 2;
            3*enodes(1:dof) - 1;
            3*enodes(1:dof)];
    edof = DOFs(:);
end