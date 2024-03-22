function varargout = TOPslicer(varargin)
%TOPSLICER M-file for TOPslicer.fig
%      TOPSLICER, by itself, creates a new TOPSLICER or raises the existing
%      singleton*.
%
%      H = TOPSLICER returns the handle to a new TOPSLICER or the handle to
%      the existing singleton*.
%
%      TOPSLICER('Property','Value',...) creates a new TOPSLICER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to TOPslicer_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TOPSLICER('CALLBACK') and TOPSLICER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TOPSLICER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TOPslicer

% Last Modified by GUIDE v2.5 02-May-2015 19:27:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TOPslicer_OpeningFcn, ...
                   'gui_OutputFcn',  @TOPslicer_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TOPslicer is made visible.
function TOPslicer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for TOPslicer
handles.output = hObject;

% Defaults
handles.pathname = '';
handles.filename = '';
handles.exportprefix = 'TOPslicer_out';
handles.dataformat = 1;
handles.cutoff = 0.5;
handles.scale = 1;
handles.transposeslice = 0;
handles.mirror = [0 0 0];
handles.isosurf = [];
handles.slicesurf = [];
handles.slicetop = [];
handles.sliceborder = [];
handles.slicecontour = [];

% Set GUI objects to their defaults
set(handles.editfilepath,'String',strcat(handles.pathname,handles.filename));
set(handles.editcutoff,'String',num2str(handles.cutoff,'%.3f'));
set(handles.editscale,'String',num2str(handles.scale,'%.3f'));
set(handles.editfileprefix,'String',handles.exportprefix);

% Parse input arguments if passed
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'input'
          set(handles.editfilepath,'String',varargin{index+1});
         case 'outputprefix'
          set(handles.editfileprefix,'String',varargin{index+1});
        end
    end
end

axes(handles.axes1)
hold on, axis equal, view(30,20), box on, axis([0 1 0 1 0 1]), rotate3d(gca,'on')
handles.slicesurf = patch('Faces',1:4,'Vertices',[0 0 1; 1 0 1; 1 1 1; 0 1 1],...
                          'FaceColor',[0.2 0.4 1],'EdgeColor','k','FaceAlpha',0.6,'Visible','off');
camlight

axes(handles.axes2)
hold on, axis equal, box on, axis([-1 1 -1 1 -1 1])
colormap(flipud(gray)), caxis([0 1])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TOPslicer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TOPslicer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function editfilepath_Callback(hObject, eventdata, handles)


function editfilepath_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbuttonbrowse_Callback(hObject, eventdata, handles)
[handles.filename,handles.pathname] = uigetfile({'*.mat'},'Load data...');
if (handles.filename~=0)
    set(handles.editfilepath,'String',strcat(handles.pathname,handles.filename));
    guidata(hObject, handles);
end


function editcutoff_Callback(hObject, eventdata, handles)
aux=sscanf(get(hObject,'String'),'%f',1);
if (isempty(aux) || aux>1 || aux<0), aux = handles.cutoff; end
set(hObject,'String',num2str(aux,'%.3f'));
handles.cutoff = aux;
guidata(hObject, handles);


function editcutoff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenuslicer_Callback(hObject, eventdata, handles)
handles.sliceplane = get(hObject,'Value');
handles = PlotTOP3slice(handles);
guidata(hObject, handles);


function popupmenuslicer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editslice_Callback(hObject, eventdata, handles)
aux=sscanf(get(hObject,'String'),'%d',1);
if isempty(aux), aux = handles.point(handles.sliceplane-1); end
if handles.sliceplane>1
    switch handles.sliceplane
        case 2
            aux = max([1 min([aux handles.size(2)]) ]);
        case 3
            aux = max([1 min([aux handles.size(1)]) ]);
        case 4
            aux = max([1 min([aux handles.size(3)]) ]);
    end
    handles.point(handles.sliceplane-1) = aux;
    set(hObject,'String',num2str(aux,'%.0f'))
    handles = PlotTOP3slice(handles);
else
    set(hObject,'String','')
end
guidata(hObject, handles);


function editslice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbuttonload_Callback(hObject, eventdata, handles)
filename = get(handles.editfilepath,'String');
aux = dir(filename);
if (length(aux)==0)
    errordlg('File not found','Load file error...')
else
    handles.data = load(filename);
    handles.vars = fieldnames(handles.data);
    set(handles.popupmenuvars,'Value',1);
    set(handles.popupmenuvars,'String',handles.vars);
end
guidata(hObject, handles);


function pushbuttonSM_Callback(hObject, eventdata, handles)
if handles.sliceplane>1
    aux = handles.point(handles.sliceplane-1) - 1;
    switch handles.sliceplane
        case 2
            aux = max([1 min([aux handles.size(2)]) ]);
        case 3
            aux = max([1 min([aux handles.size(1)]) ]);
        case 4
            aux = max([1 min([aux handles.size(3)]) ]);
    end
    handles.point(handles.sliceplane-1) = aux;
    set(handles.editslice,'String',num2str(aux,'%.0f'))
    handles = PlotTOP3slice(handles);
end
guidata(hObject, handles);


function pushbuttonSP_Callback(hObject, eventdata, handles)
if handles.sliceplane>1
    aux = handles.point(handles.sliceplane-1) + 1;
    switch handles.sliceplane
        case 2
            aux = max([1 min([aux handles.size(2)]) ]);
        case 3
            aux = max([1 min([aux handles.size(1)]) ]);
        case 4
            aux = max([1 min([aux handles.size(3)]) ]);
    end
    handles.point(handles.sliceplane-1) = aux;
    set(handles.editslice,'String',num2str(aux,'%.0f'))
    handles = PlotTOP3slice(handles);
end
guidata(hObject, handles);


function pushbuttonSMM_Callback(hObject, eventdata, handles)
if handles.sliceplane>1
    aux = handles.point(handles.sliceplane-1) - 10;
    switch handles.sliceplane
        case 2
            aux = max([1 min([aux handles.size(2)]) ]);
        case 3
            aux = max([1 min([aux handles.size(1)]) ]);
        case 4
            aux = max([1 min([aux handles.size(3)]) ]);
    end
    handles.point(handles.sliceplane-1) = aux;
    set(handles.editslice,'String',num2str(aux,'%.0f'))
    handles = PlotTOP3slice(handles);
end
guidata(hObject, handles);


function pushbuttonSPP_Callback(hObject, eventdata, handles)
if handles.sliceplane>1
    aux = handles.point(handles.sliceplane-1) + 10;
    switch handles.sliceplane
        case 2
            aux = max([1 min([aux handles.size(2)]) ]);
        case 3
            aux = max([1 min([aux handles.size(1)]) ]);
        case 4
            aux = max([1 min([aux handles.size(3)]) ]);
    end
    handles.point(handles.sliceplane-1) = aux;
    set(handles.editslice,'String',num2str(aux,'%.0f'))
    handles = PlotTOP3slice(handles);
end
guidata(hObject, handles);


function popupmenuvars_Callback(hObject, eventdata, handles)


function popupmenuvars_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenuDataFormat_Callback(hObject, eventdata, handles)
handles.dataformat = get(hObject,'Value');
guidata(hObject, handles);


function popupmenuDataFormat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbuttonplot_Callback(hObject, eventdata, handles)
fieldnumber = get(handles.popupmenuvars,'Value');
handles.x = getfield(handles.data,handles.vars{fieldnumber});
if length(size(handles.x))<3
    aux = sprintf('Variable must be a 3D density array\nwith values [0 - 1]');
    errordlg(aux,'Data error...')
else
    % MATLAB's default format is column-major, i.e. [Y X Z]
    if handles.dataformat==2 % row-major storage convention [X Y Z]
        handles.x = permute(handles.x,[2 1 3]); % MATLAB expects column-major
    elseif handles.dataformat==3 % TOP3D's storage convention [-Z X -Y]
        handles.x = flip(flip(permute(handles.x,[3 2 1]),3),1);
    end
    handles.x = ApplyMarginSymmetry(handles.x,handles.mirror);
    handles.size = size(handles.x) - 2;
    
    handles = PlotTOP3(handles);
    
    set(handles.popupmenuslicer,'Value',1)
    set(handles.textslice,'String','')
    set(handles.editslice,'String','')
    set(handles.checkboxtransposeslice,'Value',0)
    handles.transposeslice = 0;
    handles.sliceplane = 1;
    handles.point=[1 1 1];
    handles = PlotTOP3slice(handles);
end
guidata(hObject, handles);


function checkboxtransposeslice_Callback(hObject, eventdata, handles)
handles.transposeslice = get(hObject,'Value');
if (~isempty(handles.isosurf) && handles.sliceplane>1)
    handles = PlotTOP3slice(handles);
end
guidata(hObject, handles);


function checkboxMXl_Callback(hObject, eventdata, handles)
aux = get(hObject,'Value');
if (aux)
    set(handles.checkboxMXu,'Value',0);
    handles.mirror(1) = -1;
else
    handles.mirror(1) = 0;
end
guidata(hObject, handles);


function checkboxMXu_Callback(hObject, eventdata, handles)
aux = get(hObject,'Value');
if (aux)
    set(handles.checkboxMXl,'Value',0);
    handles.mirror(1) = 1;
else
    handles.mirror(1) = 0;
end
guidata(hObject, handles);


function checkboxMYl_Callback(hObject, eventdata, handles)
aux = get(hObject,'Value');
if (aux)
    set(handles.checkboxMYu,'Value',0);
    handles.mirror(2) = -1;
else
    handles.mirror(2) = 0;
end
guidata(hObject, handles);


function checkboxMYu_Callback(hObject, eventdata, handles)
aux = get(hObject,'Value');
if (aux)
    set(handles.checkboxMYl,'Value',0);
    handles.mirror(2) = 1;
else
    handles.mirror(2) = 0;
end
guidata(hObject, handles);


function checkboxMZl_Callback(hObject, eventdata, handles)
aux = get(hObject,'Value');
if (aux)
    set(handles.checkboxMZu,'Value',0);
    handles.mirror(3) = -1;
else
    handles.mirror(3) = 0;
end
guidata(hObject, handles);


function checkboxMZu_Callback(hObject, eventdata, handles)
aux = get(hObject,'Value');
if (aux)
    set(handles.checkboxMZl,'Value',0);
    handles.mirror(3) = 1;
else
    handles.mirror(3) = 0;
end
guidata(hObject, handles);


function editscale_Callback(hObject, eventdata, handles)
aux=sscanf(get(hObject,'String'),'%f',1);
if (isempty(aux) || aux<=0), aux = handles.scale; end
set(hObject,'String',num2str(aux,'%.4g'));
handles.scale = aux;
guidata(hObject, handles);


function editscale_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editdimX_Callback(hObject, eventdata, handles)


function editdimX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editdimY_Callback(hObject, eventdata, handles)


function editdimY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editdimZ_Callback(hObject, eventdata, handles)


function editdimZ_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editfileprefix_Callback(hObject, eventdata, handles)
fileprefix = get(hObject,'String');
if isempty(fileprefix)
    set(hObject,'String',handles.exportprefix);
else
    handles.exportprefix = fileprefix;
    guidata(hObject, handles);
end
    

function editfileprefix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbuttonexportX3D_Callback(hObject, eventdata, handles)
if isfield(handles,'x')
    h = waitbar(0,'Please wait...','Name','Exporting to X3D & X3DOM');
    hw=findobj(h,'Type','Patch');
    set(hw,'EdgeColor',[0 0 1],'FaceColor',[0 0 1]) % Just for fun
    color = [1 0 0]; % Default color
    [nely,nelx,nelz]=size(handles.x);
    Px = (0:nelx-1)-0.5;
    Py = (0:nely-1)-0.5;
    Pz = (0:nelz-1)-0.5;
    [FACE,NODE] = isosurface(Px,Py,Pz,handles.x,handles.cutoff);
    FACE = fliplr(FACE);
    waitbar(1/3,h);
    WriteX3D(NODE,FACE,strcat(handles.exportprefix,'.x3d'),color,handles.scale)
    waitbar(2/3,h);
    WriteX3DOM(NODE,FACE,strcat(handles.exportprefix,'.html'),color,handles.scale)
    waitbar(1,h); close(h);
else
    errordlg('No data is loaded','Export X3D/X3DOM error...')
end


function pushbuttonexportSTL_Callback(hObject, eventdata, handles)
if isfield(handles,'x')
    [nely,nelx,nelz] = size(handles.x);
    Px = handles.scale * ((0:nelx-1)-0.5);
    Py = handles.scale * ((0:nely-1)-0.5);
    Pz = handles.scale * ((0:nelz-1)-0.5);
    [FACE,NODE] = isosurface(Px,Py,Pz,handles.x,handles.cutoff);
    metadata = [nelx-2 nely-2 nelz-2 handles.cutoff];
    symm = [handles.mirror==-1; handles.mirror==1];
    FACEtoSTLbin(strcat(handles.exportprefix,'.stl'),NODE,fliplr(FACE));
else
    errordlg('No data is loaded','Export STL error...')
end

%% === ADDITIONAL FUNCTIONS ===============================================

function [handles]=PlotTOP3(handles)
% Display 3D topology in iso-view
[nely,nelx,nelz] = size(handles.x);
Px = handles.scale * ((0:nelx-1)-0.5);
Py = handles.scale * ((0:nely-1)-0.5);
Pz = handles.scale * ((0:nelz-1)-0.5);
axes(handles.axes1)
if ~isempty(handles.isosurf)
    delete(handles.isosurf)
    handles.isosurf = [];
end
axis(handles.scale * [-1 nelx+1 -1 nely+1 -1 nelz+1])
fv = isosurface(Px,Py,Pz,handles.x,handles.cutoff); % facet-vertex struct
handles.isosurf = patch(fv,'FaceColor','r','EdgeColor','none',...
                        'FaceLighting','gouraud','AmbientStrength',0.5);
drawnow
dims = max(fv.vertices) - min(fv.vertices);
if isempty(dims)
    set(handles.editdimX,'String','--');
    set(handles.editdimY,'String','--');
    set(handles.editdimZ,'String','--');
else
    set(handles.editdimX,'String',num2str(dims(1),'%.4g'));
    set(handles.editdimY,'String',num2str(dims(2),'%.4g'));
    set(handles.editdimZ,'String',num2str(dims(3),'%.4g'));
end


function [handles]=PlotTOP3slice(handles)
% Plot slice plane in 3D plot
if (isempty(handles.isosurf) || handles.sliceplane==1)
    set(handles.slicesurf,'Visible','off')
    axes(handles.axes2), cla
    set(handles.textslice,'String','');
    set(handles.editslice,'String','');
else
    set(handles.slicesurf,'Visible','on')
    if (handles.sliceplane==2)
        set(handles.editslice,'String',num2str(handles.point(1),'%.0f'));
        textslice = sprintf('[001 - %03.0f]',handles.size(2));
        aux = squeeze(handles.x(:,handles.point(1)+1,:))'; % Default is ^T
        if handles.transposeslice, aux=aux'; end
        x = handles.point(1)-0.5; y = handles.size(1); z = handles.size(3);
        set(handles.slicesurf,'Vertices',[x 0 0; x y 0; x y z; x 0 z]*handles.scale);
    elseif (handles.sliceplane==3)
        set(handles.editslice,'String',num2str(handles.point(2),'%.0f'));
        textslice = sprintf('[001 - %03.0f]',handles.size(1));
        aux = squeeze(handles.x(handles.point(2)+1,:,:))'; % Default is ^T
        if handles.transposeslice, aux=aux'; end
        x = handles.size(2); y = handles.point(2)-0.5; z = handles.size(3);
        set(handles.slicesurf,'Vertices',[0 y 0; x y 0; x y z; 0 y z]*handles.scale);
    elseif (handles.sliceplane==4)
        set(handles.editslice,'String',num2str(handles.point(3),'%.0f'));
        textslice = sprintf('[001 - %03.0f]',handles.size(3));
        aux = squeeze(handles.x(:,:,handles.point(3)+1));
        if handles.transposeslice, aux=aux'; end
        x = handles.size(2); y = handles.size(1); z = handles.point(3)-0.5;
        set(handles.slicesurf,'Vertices',[0 0 z; x 0 z; x y z; 0 y z]*handles.scale);
    end
    set(handles.textslice,'String',textslice);
    handles = PlotSlice(handles,aux);
end


function [handles]=PlotSlice(handles,aux)
% Plots the slice and contour
Nx = size(aux,2); Ny = size(aux,1);
axes(handles.axes2), cla
handles.slicetop = imagesc(-0.5,-0.5,aux); % Shift image to account for margin
handles.sliceborder = plot3((Nx-2)*[0 1 1 0 0],(Ny-2)*[0 0 1 1 0],...
                            0.25*[1 1 1 1 1],'Color',[0.0 0.6 0.90],'LineWidth',1.5);
if (handles.cutoff>=min(min(aux)) && handles.cutoff < max(max(aux)))
    [~,handles.slicecontour] = contour(-0.5:Nx-1.5,-0.5:Ny-1.5,aux,...
                                       handles.cutoff*[1 1],'r','LineWidth',1.5);
end
axis([0 Nx-2 0 Ny-2]+0.5*[-1 1 -1 1])


function [rho]=ApplyMarginSymmetry(rho,mirror)
% Applies symmetry (if specified) and adds a zero-density margin
[nely,nelx,nelz]=size(rho);
if mirror(1)<0 % Mirror on Xmin
    aux=rho;
    rho=zeros(nely,2*nelx,nelz);
    rho(:,1:nelx,:)=flipdim(aux,2);
    rho(:,nelx+1:end,:)=aux;
    nelx=2*nelx;
elseif mirror(1)>0 % Mirror on Xmax
    aux=rho;
    rho=zeros(nely,2*nelx,nelz);
    rho(:,1:nelx,:)=aux;
    rho(:,nelx+1:end,:)=flipdim(aux,2);
    nelx=2*nelx;
end
if mirror(2)<0 % Mirror on Ymin
    aux=rho;
    rho=zeros(2*nely,nelx,nelz);
    rho(1:nely,:,:)=flipdim(aux,1);
    rho(nely+1:end,:,:)=aux;
    nely=2*nely;
elseif mirror(2)>0 % Mirror on Ymax
    aux=rho;
    rho=zeros(2*nely,nelx,nelz);
    rho(1:nely,:,:)=aux;
    rho(nely+1:end,:,:)=flipdim(aux,1);
    nely=2*nely;
end
if mirror(3)<0 % Mirror on Zmin
    aux=rho;
    rho=zeros(nely,nelx,2*nelz);
    rho(:,:,1:nelz)=flipdim(aux,3);
    rho(:,:,nelz+1:end)=aux;
    nelz=2*nelz;
elseif mirror(3)>0 % Mirror on Zmax
    aux=rho;
    rho=zeros(nely,nelx,2*nelz);
    rho(:,:,1:nelz)=aux;
    rho(:,:,nelz+1:end)=flipdim(aux,3);
    nelz=2*nelz;
end
% Add margin layer
aux=rho;
rho=zeros(nely+2,nelx+2,nelz+2);
rho(2:end-1,2:end-1,2:end-1)=aux;