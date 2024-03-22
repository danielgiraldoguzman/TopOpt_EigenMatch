function PlotTopology(mesh,ro,fc,ec)

%% PLOT ELEMENTS
n = mesh.conect(:,1:8);
% Conterclockwise selection w.r.t the normal vector of each face
Faces = [n(:,1) n(:,2) n(:,3) n(:,4);   % Face 1 - Bottom
         n(:,5) n(:,6) n(:,7) n(:,8);   % Face 2 - Top
         n(:,1) n(:,2) n(:,6) n(:,5);   % Face 3 - Right Front
         n(:,3) n(:,4) n(:,8) n(:,7);   % Face 4 - Left Back
         n(:,2) n(:,3) n(:,7) n(:,6);   % Face 5 - Right Back
         n(:,1) n(:,4) n(:,8) n(:,5)];  % Face 6 - Left Front

nfaces = size(Faces,1);
alphavertex = zeros(nfaces,1);
for i=1:mesh.nel
    alphavertex(i:mesh.nel:nfaces)=ro(i);
end

if nargin < 4
    fc = 'k';   % Face Color
    ec = 'w';   % Edge Color
end

patch(  'Vertices',mesh.ncoord,...
        'Faces',Faces,...
        'FaceVertexAlphaData',alphavertex,...
        'AlphaDataMapping','none',...
        'FaceColor',fc,'FaceAlpha','flat',...
        'EdgeColor',ec,'EdgeAlpha','flat')

view(3);
axis equal; axis tight;
