function xPhys = TOPslicerData(Mapping,rho)
% Export data for TOPslicer tool
% This function will only work for strcutured meshes with cubic elements

% Input data
if size(rho,2) > 1
    rho = rho(:,end); % If multiple data entries, select the last entry
end

Cx = unique(Mapping(:,1));
Cy = unique(Mapping(:,2));
Cz = unique(Mapping(:,3));
xPhys = zeros(length(Cx),length(Cy),length(Cz));

% Generate idexes according to each element center location
for loop = 1:length(Cx)
    Mapping(Mapping(:,1) == Cx(loop),1) = loop;
end
for loop = 1:length(Cy)
    Mapping(Mapping(:,2) == Cy(loop),2) = loop;
end
for loop = 1:length(Cz)
    Mapping(Mapping(:,3) == Cz(loop),3) = loop;
end

% Write data in row-major ordering
for loop = 1:length(rho)
    xPhys(Mapping(loop,1),Mapping(loop,2),Mapping(loop,3)) = rho(loop);
end