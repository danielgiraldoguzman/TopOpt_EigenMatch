function [FEM,MatFactorK,MatFactorM,MatFactorC] = WriteUpdatedMTX(ini,ro,FEM)

% Interpolate material properties
[MatFactorK,MatFactorM,MatFactorC] = MaterialModel(ini.MatModel,ro,ini.p);

% Apply material interpolation model to every element matrix
k = MatFactorK'.*FEM.k;
m = MatFactorM'.*FEM.m;
c = MatFactorC'.*FEM.c;

% Global Matrix Assembly
FEM.K = sparse(FEM.mesh.Ai,FEM.mesh.Aj,k,FEM.mesh.dof,FEM.mesh.dof);
FEM.M = sparse(FEM.mesh.Ai,FEM.mesh.Aj,m,FEM.mesh.dof,FEM.mesh.dof);
FEM.C = sparse(FEM.mesh.Ai,FEM.mesh.Aj,c,FEM.mesh.dof,FEM.mesh.dof);

% Convert Sparse matrices to ABAQUS Matrix Input Format
STIF = Sparse2MIF(FEM.mesh,FEM.K);
MASS = Sparse2MIF(FEM.mesh,FEM.M);
DAMP = Sparse2MIF(FEM.mesh,FEM.C);

% Write FEA matrices to be used in ABAQUS
ID_Stiff = fopen('STIF.mtx','w');
fprintf(ID_Stiff,'%d, %d, %d, %d, %.15E \n', STIF');
fclose(ID_Stiff);

ID_Mass = fopen('MASS.mtx','w');
fprintf(ID_Mass,'%d, %d, %d, %d, %.15E \n', MASS');
fclose(ID_Mass);

ID_Damp = fopen('DAMP.mtx','w');
fprintf(ID_Damp,'%d, %d, %d, %d, %.15E \n', DAMP');
fclose(ID_Damp);

%% Generate interpolation models plots
% ro = 0:0.01:1;
% [roK,roM] = MaterialModel('RAMP',ro,3,1);
% figure; plot(ro,roM,ro,roK); grid; legend('Mass','Stiffness')

end

