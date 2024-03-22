function [rho_d,drho_d] = DensityFilter(mesh,Rf,rho)
% Density Filter by Bruns and Tortorelli (2001) and Bourdin (2001). 
%
% INPUT:
% mesh      % Structure with mesh information
% Rf        % Filter radius
% rho       % Pseudo-density of element i
%
% OUTPUT:
% rho_d     % Density-Filtered pseudo-density
% drho_d    % Derivative of density-filtered pseudo-density

rho_d = zeros(mesh.nel,1);  % Pre-allocation
drho_d = zeros(mesh.nel);   % Pre-allocation

parfor e = 1:mesh.nel
    % Pseudo-density filtered
    Wei = Rf - mesh.Rei(e,:);
    Wei(Wei<0) = 0; % Only elements within the neighborhood
%     Wei(1,e) = 0;   % Filtred density won't depend on its own density
    Num = sum(Wei.*mesh.Ve'.*rho');
    Den = sum(Wei.*mesh.Ve');
    rho_d(e) = Num/Den;
    
    % Derivative of element e w.r.t. elemenent k within the neighborhood
    WekVek = Wei.*mesh.Ve';
    drho_d(e,:) = WekVek/Den; % Derivate of element e w.r.t. k
end

% End of function
end