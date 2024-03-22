function [rho_d,rho_h,drho_d,drho_h] = ProjectionFilter(ini,FEM,rho)

fprintf('Filtering pseudo-densities...  \n')
%% Density Filter by Bruns and Tortorelli (2001) and Bourdin (2001)
[rho_d,drho_d] = DensityFilter(FEM.mesh,ini.Rdf,rho);
% Check for forced elements
if any(ini.Sets.VoidNonDesign); rho_d(ini.Sets.VoidNonDesign) = 1e-6; end
if any(ini.Sets.SolidNonDesign); rho_d(ini.Sets.SolidNonDesign) = 1; end
Vol_d = (FEM.mesh.Ve'*rho_d)/FEM.TotalVolume*100;
fprintf('Volume after Density:   %.4g%% \n',Vol_d)

% fprintf('Densitiy filter turned off \n')
% rho_d = rho;
% drho_d = speye(FEM.mesh.nel);

%% Projection Filter based on Heaviside functions
fprintf('Beta (Heaviside): %d \n',ini.beta)
eta = VolumePreserving(rho_d,ini.beta,FEM.mesh.Ve);
fprintf('eta parameter:    %.3g \n',eta)
[rho_h,drho_h] = HeavisideFilter(ini.beta,eta,rho_d);
% Check for forced elements
if any(ini.Sets.VoidNonDesign); rho_h(ini.Sets.VoidNonDesign) = 1e-6; end
if any(ini.Sets.SolidNonDesign); rho_h(ini.Sets.SolidNonDesign) = 1; end
Vol_h = (FEM.mesh.Ve'*rho_h)/FEM.TotalVolume*100;
fprintf('Volume after Heaviside: %.4g%% \n',Vol_h)

% fprintf('Heaviside filter turned off \n')
% rho_h = rho_d;
% drho_h = ones(FEM.mesh.nel,1);

fprintf('Filtering completed\n')
% End of function
end