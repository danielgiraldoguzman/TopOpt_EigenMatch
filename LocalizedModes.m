function [idx_Localized,Soft2solid] = LocalizedModes(ini,FEM,ro,EigV)
% Localized modes identification module

% Apply penalization scheme to properly evaluate elements density
ro = MaterialModel(ini.MatModel,ro,ini.p);

% find solid and soft elements, and its nodes
idx_SolidElements = find(ro > ini.ThresholdUP);  % solid elements
idx_SoftElements = find(ro < ini.ThreshoDOWN);   % soft elements
SolidNodes = FEM.mesh.conect(idx_SolidElements,:); % solid nodes
SolidNodes = unique(SolidNodes); % organized nodal information
SoftNodes = FEM.mesh.conect(idx_SoftElements,:);   % soft nodes
SoftNodes = unique(SoftNodes);   % organized nodal information

% analize solid elements against soft elements
Soft2solid = zeros(size(EigV,3),2);
for Mode=1:size(EigV,3)
    SolidDisp = abs( EigV(SolidNodes,:,Mode) );
    SolidDisp = SolidDisp(:);
    SolidDisp(SolidDisp == 0) = []; % Remove nodes with zero-like B
    SoftDisp = abs( EigV(SoftNodes,:,Mode) );
    SoftDisp = SoftDisp(:);
    SoftDisp(SoftDisp == 0) = [];   % Remove nodes with zero-like BC
    
    if any(SoftDisp) && any(SolidDisp)
        Soft2solid(Mode,1) = prctile(SoftDisp,90)/prctile(SolidDisp,90);
        Soft2solid(Mode,2) = max(SoftDisp,[],'all')/max(SolidDisp,[],'all');
        Soft2solid(Mode,3) = mean(SoftDisp)/mean(SolidDisp);
    end
end
Soft2solid(:,end+1) = sum(Soft2solid,2);

% Soft2Solid ratios really really large, are localized modes
RatioCriteria = 1e3; % how much soft is allowwed to move compared to solid
idx_Localized = find( Soft2solid(:,end) > RatioCriteria );

% End of function
end