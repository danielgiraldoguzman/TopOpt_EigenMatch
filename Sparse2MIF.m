function [MTX] = Sparse2MIF(mesh,K)
% Convert MATLAB Sparse Matrix into ABAQUS Matrix Input Format

% Extract indexing information
K = tril(K,0);     % Delete upper triangular data (due to symmetry)
[I,J,V] = find(K); % Extract indexes and data

% Create Matrix Input Format indexes from mesh
nodes = [1:mesh.nnod; 1:mesh.nnod; 1:mesh.nnod];
dofs = repmat([1;2;3],mesh.nnod,1);
idxMIF = [nodes(:) dofs];

% Convert idex from Sparse indexing to Matrix Input Format
I1 = idxMIF(I,:);
I2 = idxMIF(J,:);

% Matrix Input Format indexing
MTX = [I1 I2 V]; % New Matrix, ready to be input in ABAQUS

end