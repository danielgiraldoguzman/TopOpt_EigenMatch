function []=FACEtoSTLbin(filename,NODE,FACE3)

%% === INPUT DEFAULTS =====================================================
if isempty(FACE3)
    error('Facet data cannot be empty.')
elseif ( min(min(FACE3))<1 || max(max(FACE3))>size(NODE,1) )
    error('Facet data out-of-range.')
elseif isempty(filename)
    error('Output filename cannot be empty.')
end

%% === INPUT PARSE ========================================================
[filepath,filename,fileext] = fileparts(filename);
if isempty(fileext)
    fileext = '.stl';
elseif ~strcmpi(fileext,'.stl')
    fprintf('Warning - Specified extension is not *.stl\n')
end
if size(NODE,2)~=3
    error('NODE must be a [Nn x 3] matrix.')
elseif any(~isfinite(NODE))
    error('NODE values cannot be NaN or Inf.')
end
if isempty(FACE3)
    error('Facet information cannot be empty')
elseif size(FACE3,2)~=3
    error('FACE3 must be a [Nf x 3] matrix.')
end

%% === BEGIN OUTPUT =======================================================
fid = fopen(strcat(filepath,filename,fileext),'wb+');
if isequal(fid,-1), error('Unable to write to %s.', filename), end

Nt = size(FACE3,1);
time = clock;

title = sprintf('Exported from TOPslicer ~ %s %02d:%02d',date,time(4:5));
fprintf(fid,'%-80s',title);

fwrite(fid,Nt,'uint32'); % Number of trias

h = waitbar(0,'Please wait...','Name','Exporting to binary STL');
hw=findobj(h,'Type','Patch');
set(hw,'EdgeColor',[1 0.84 0.05],'FaceColor',[1 0.84 0.05]) % Just for fun
barlvl = 0.1; barnum = floor(barlvl*Nt);
for i=1:size(FACE3,1)
    V1 = diff(NODE(FACE3(i,[1 2]),:),1);
    V2 = diff(NODE(FACE3(i,[1 3]),:),1);
    N = cross(V1,V2);
    N = N / norm(N);
    WriteTriaBIN(fid,N,NODE(FACE3(i,:),:));
    if (i>barnum)
        waitbar(barlvl,h);
        while (i>barnum)
            barlvl = barlvl + 0.1;
            barnum = floor(barlvl*Nt);
        end
    end
end
close(h); % Close the waitbar

% Close the file stream
fclose(fid);

function []=WriteTriaBIN(fid,N,NODE)
% Write tria to fid with normal D and node coordinates
% NODE = [x1 y1 z1; x2 y2 z2; x3 y3 z3];
fwrite(fid,N,'float32');
fwrite(fid,NODE','float32');
fwrite(fid,0,'int16'); % RGB color
return