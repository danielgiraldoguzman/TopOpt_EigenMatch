% Write negative FEA matrices
% This is to import matrices to Abaqus, such that:
% K = K_0 + K_ext - K_0 and M = M_0 + M_ext - M_0
% K_0 is the matrix generated automatically in ABAQUS, but its influence
% must be deleted, then, a negative -K_0 is introduced to cancel it out.

% Create negative matrices
STIF_NEG = [Kmif(:,1:4) -1*Kmif(:,5)];
MASS_NEG = [Mmif(:,1:4) -1*Mmif(:,5)];
DAMP_NEG = [Cmif(:,1:4) -1*Cmif(:,5)];

% Write FEA matrices to be used in ABAQUS
ID_Stiff = fopen('STIF_NEG.mtx','w');
fprintf(ID_Stiff,'%d, %d, %d, %d, %.15E \n', STIF_NEG');
fclose(ID_Stiff);

ID_Mass = fopen('MASS_NEG.mtx','w');
fprintf(ID_Mass,'%d, %d, %d, %d, %.15E \n', MASS_NEG');
fclose(ID_Mass);

ID_Damp = fopen('DAMP_NEG.mtx','w');
fprintf(ID_Damp,'%d, %d, %d, %d, %.15E \n', DAMP_NEG');
fclose(ID_Damp);

fprintf('Negative Matrices created.\n\n')