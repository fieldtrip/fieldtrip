function varargout = metchgui(varargin)
%    alldata = metchgui(node,elem,points,pface) or metchgui(volume,points,pface)
%
%    A GUI to register a point cloud to a mesh or volumetric image
%
%    author: Qianqian Fang <q.fang at neu.edu>
%    date: 12/16/2008
%
%   parameters:
%        node: node coordinate of the surface mesh (nn x 3)
%        elem: element list of the surface mesh (3 columns for
%              triangular mesh, 4 columns for cubic surface mesh)
%        points: the coordinates (3 columns for x/y/z) of the
%              point cloud which you want to register
%        pface:trianglular surface defined on the point cloud.
%              pface is optional; if presents, metch will display
%              a surface object instead of a point cloud.
%
%   the input can also be two parameters in form of metchgui(volume,points),
%    where volume is a 3D image (array).
%
%   outputs:
%        alldata: a structrure containing all processing outputs
%        the fields include:
%         .node: the input node
%         .elem: the input surface mesh elements
%         .volume: if the input volumetric image
%         .A0: the affine rotation for selected point pairs (after Initialize)
%         .b0: the affine translation for selected point pairs (after Initialize)
%         .A: the affine rotation for the point cloud (after Optimize)
%         .b: the affine translation for the point cloud (after Optimize)
%         .points: the input point cloud
%         .pointsinit: the point cloud after initialization
%         .pointsopt: the point cloud after optimization
%         .pointsproj: the point cloud after projecting to the surface
%         .initplot: the handle to the point cloud plot after init
%         .optplot: the handle to the point cloud plot after optimization
%         .projplot: the handle to the point cloud plot after projection
%
%   If user supplys an output variable, the GUI will not return until the
%   user hits the "close" button or close the window; if user does not
%   supply any output, the call will return immediately; any data user
%   intends to save, he has to click on "Save Session" button and provides
%   a mat-file file name. A single structure named "metchsession" will be
%   stored in this file.
%
%   example: (meshasphere/meshunitsphere are defined in iso2mesh http://iso2mesh.sf.net)
%
%       [noderef,faceref,elemref]=meshunitsphere(0.08,10);
%       [no,fc]=removeisolatednode(noderef(:,1:3),faceref(:,1:3));
%       [node,face,elem]=meshasphere([10 20 15],3,0.5,10);
%       [no2,fc2]=removeisolatednode(node(:,1:3),face(:,1:3));
%       alldata = metchgui(no,fc,no2);
%       % or alldata = metchgui(no,fc,no2,fc2);
%
%   Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
%   this function is part of "metch" toobox, see COPYING for license

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @metchgui_OpeningFcn, ...
                   'gui_OutputFcn',  @metchgui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}, 'hasoutput');
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before metchgui is made visible.
function metchgui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.hasoutput = 0;
if (isempty(varargin))
    fprintf(1, 'Metch GUI must be called with parameters:\nFormat: alldata = metchgui(node,elem,points,pface);\n');
    close(handles.MetchGUI);
    return
end
if (ischar(varargin{end}) && strcmp(varargin{end}, 'hasoutput'))
    handles.hasoutput = 1;
    varargin(end) = [];
end
handles.output = hObject;

set(handles.btAddMeshPt, 'userdata', [handles.axMesh, handles.btAddMeshPt, handles.btAddCloudPt]);
set(handles.btAddCloudPt, 'userdata', [handles.axMesh, handles.btAddCloudPt, handles.btAddMeshPt]);

% if uses supplied 2 input variables, assume a volume image and a point cloud
if (isnumeric(varargin{1}) && length(size(varargin{1})) == 3)
    vol = varargin{1};
    pt = varargin{2};

    dat.volume = vol;
    dat.points = pt;

    slice = round(size(vol, 3) / 2);
    hs = imagesc(vol(:, :, slice), 'parent', handles.axMesh);
    set(handles.slPos, 'max', size(vol, 3), 'min', 1, 'value', slice);

    if (length(varargin) >= 3)
        pface = varargin{3};
        dat.pface = pface;
        trisurf(pface(:, 1:3), pt(:, 1), pt(:, 2), pt(:, 3), 'parent', handles.axPoints);
    else
        % plot3(pt(:,1),pt(:,2),pt(:,3),'.','parent',handles.axPoints);
        ptcolor = (pt - repmat(min(pt), size(pt, 1), 1)) ./ repmat(max(pt) - min(pt), size(pt, 1), 1);
        drawnow;
        scatter3(pt(:, 1), pt(:, 2), pt(:, 3), 3, ptcolor, 'filled');
    end
    set(handles.MetchGUI, 'userdata', dat);

    axis(handles.axMesh, 'equal');
    axis(handles.axPoints, 'equal');
    axis(handles.axMesh, 'off');
    grid(handles.axPoints, 'on');
    % axis(handles.axPoints,'off');

    set(handles.axMesh, 'tag', 'axMesh');
    set(handles.axPoints, 'tag', 'axPoints');
    set(handles.slPos, 'visible', 'on');
    set(handles.lbZPos, 'visible', 'on');
    rotate3d(handles.axPoints, 'on');
    rotate3d(gcf, 'on');
end

% if uses supplied 3 input variables, assume a surface mesh and a point cloud/surface
if (length(varargin) >= 3 && length(size(varargin{1})) == 2)
    node = varargin{1};
    elem = varargin{2};
    pt = varargin{3};

    dat.node = node;
    dat.elem = elem;
    dat.points = pt;
    if (length(varargin) >= 4)
        pface = varargin{4};
        dat.pface = pface;
    end
    set(handles.MetchGUI, 'userdata', dat);

    drawinit(handles, dat);
end
guidata(hObject, handles);

if (handles.hasoutput)
    uiwait(handles.MetchGUI);
end

% ---------------------------------------------------------------------------
function drawinit(handles, dat)

node = dat.node;
elem = dat.elem;
pt = dat.points;

hs = trisurf(elem, node(:, 1), node(:, 2), node(:, 3), 'parent', handles.axMesh);
% set(hs,'linestyle','none');
% set(hs,'facecolor','b','facealpha',0.8);

% plot3(pt(:,1),pt(:,2),pt(:,3),'.','parent',handles.axPoints);
if (~isfield(dat, 'pface'))
    ptcolor = (pt - repmat(min(pt), size(pt, 1), 1)) ./ repmat(max(pt) - min(pt), size(pt, 1), 1);
    drawnow;
    scatter3(pt(:, 1), pt(:, 2), pt(:, 3), 3, ptcolor, 'filled');
else
    trisurf(dat.pface, pt(:, 1), pt(:, 2), pt(:, 3), 'parent', handles.axPoints);
end
% hold(handles.axPoints,'on');
% plot3(pt(5:7,1),pt(5:7,2),pt(5:7,3),'ro','parent',handles.axPoints);

axis(handles.axMesh, 'equal');
axis(handles.axPoints, 'equal');
axis(handles.axMesh, 'off');
% axis(handles.axPoints,'off');

set(handles.axMesh, 'tag', 'axMesh');
set(handles.axPoints, 'tag', 'axPoints');

rotate3d(handles.axPoints, 'on');
rotate3d(handles.axMesh, 'on');
rotate3d(gcf, 'on');

% ---------------------------------------------------------------------------
function varargout = metchgui_OutputFcn(hObject, eventdata, handles)
if (length(handles) && isfield(handles, 'output') && isfield(handles, 'hasoutput') && handles.hasoutput)
    handles.output = get(handles.MetchGUI, 'userdata');
    varargout{1} = handles.output;
    close(handles.MetchGUI);
end

% ---------------------------------------------------------------------------
function lbMesh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

% ---------------------------------------------------------------------------
function lbPoints_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

% ---------------------------------------------------------------------------
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

% ---------------------------------------------------------------------------
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

% ---------------------------------------------------------------------------
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

% ---------------------------------------------------------------------------
function isSelect_Callback(hObject, eventdata, handles)
if (get(hObject, 'Value'))
    datacursormode(gcf, 'on');
    set(datacursormode(gcf), 'UpdateFcn', @myupdatefcn);
else
    datacursormode(gcf, 'off');
    rotate3d(gcf, 'on');
end

% ---------------------------------------------------------------------------
% the respond function when there is a data-tip to popup
% ---------------------------------------------------------------------------
function txt = myupdatefcn(empt, event_obj)
pos = get(event_obj, 'Position');
idx =  get(event_obj, 'DataIndex');

if (length(pos) == 3)
    txt = {['x: ', num2str(pos(1))], ...
           ['y: ', num2str(pos(2))], ['z: ', num2str(pos(3))], ['index:', num2str(idx)]};
elseif (length(pos) == 2)
    txt = {['x: ', num2str(pos(1))], ['y: ', num2str(pos(2))]};
end
targetup = get(get(event_obj, 'Target'), 'parent');
set(targetup, 'userdata', struct('pos', pos, 'idx', idx));
if (targetup == findobj('tag', 'axMesh'))
    set(findobj('tag', 'btAddMeshPt'), 'enable', 'on');
    set(findobj('tag', 'btAddCloudPt'), 'enable', 'off');
elseif (targetup == findobj('tag', 'axPoints'))
    set(findobj('tag', 'btAddMeshPt'), 'enable', 'off');
    set(findobj('tag', 'btAddCloudPt'), 'enable', 'on');
end

% ---------------------------------------------------------------------------
function bInit_Callback(hObject, eventdata, handles)
mapto = get(handles.lbMesh, 'userdata');
mapfrom = get(handles.lbPoints, 'userdata');

maptoidx = get(handles.txMapTo, 'userdata');
mapfromidx = get(handles.txMapFrom, 'userdata');

if (length(mapto) < 4 | length(mapfrom) < 4)
    msgbox('You have to select >3 points from the point cloud plot and corresponding points from the mesh', 'Error', 'error');
    return
end
[A, b] = affinemap(mapfrom, mapto);
dat = get(handles.MetchGUI, 'userdata');
dat.A0 = A;
dat.b0 = b;

newpt = (A * dat.points' + repmat(b(:), 1, size(dat.points, 1)))';
dat.pointsinit = newpt;

hold(handles.axMesh, 'on');
if (isfield(dat, 'initplot'))
    delete dat.initplot;
    dat.initplot = 0;
end
dat.initplot = plot3(newpt(:, 1), newpt(:, 2), newpt(:, 3), 'r.', 'parent', handles.axMesh);

dat.fromidx = mapfromidx;
dat.toidx = maptoidx;
set(handles.MetchGUI, 'userdata', dat);

% ---------------------------------------------------------------------------
function axPoints_ButtonDownFcn(hObject, eventdata, handles)
if (get(handles.btSelect, 'value') == 1)
    pp = getCursorInfo(datacursormode(gcf));
    str = [get(handles.lbPoints, 'string'); mat2str(pp.Position)];
    set(handles.lbPoints, 'string', str);
end

% ---------------------------------------------------------------------------
function addselectedpt(pos, idx, lb)
if (isempty(get(lb, 'value')))
    set(lb, 'value', 1);
end
listpt = get(lb, 'string');
listpt{end + 1} = [num2str(idx) ':' mat2str(pos)];
set(lb, 'string', listpt);
set(lb, 'userdata', [get(lb, 'userdata'); pos]);

% ---------------------------------------------------------------------------
function btAddMeshPt_Callback(hObject, eventdata, handles)
dat = get(handles.axMesh, 'userdata');

dat0 = get(handles.MetchGUI, 'userdata');
if (isfield(dat0, 'volume'))
    dat.pos(:, 3) = get(handles.slPos, 'value');
end

if (isfield(dat, 'pos'))
    addselectedpt(dat.pos, dat.idx, handles.lbMesh);
    set(handles.txMapTo, 'userdata', [get(handles.txMapTo, 'userdata'); dat.idx]);
    dat0.mapto = get(handles.lbMesh, 'userdata');
    dat0.maptoidx = get(handles.txMapTo, 'userdata');
    set(handles.MetchGUI, 'userdata', dat0);
else
    msgbox('No point was selected. Please click on "Select" and select a point on the mesh or point cloud', 'Error', 'error');
    return
end

% ---------------------------------------------------------------------------
function btAddCloudPt_Callback(hObject, eventdata, handles)
dat = get(handles.axPoints, 'userdata');
dat0 = get(handles.MetchGUI, 'userdata');
if (isfield(dat, 'pos'))
    addselectedpt(dat.pos, dat.idx, handles.lbPoints);
    set(handles.txMapFrom, 'userdata', [get(handles.txMapFrom, 'userdata'); dat.idx]);
    dat0.mapfrom = get(handles.lbPoints, 'userdata');
    dat0.mapfromidx = get(handles.txMapFrom, 'userdata');
    set(handles.MetchGUI, 'userdata', dat0);
else
    msgbox('No point was selected. Please click on "Select" and select a point on the mesh or point cloud', 'Error', 'error');
    return
end

% ---------------------------------------------------------------------------
function btOptimize_Callback(hObject, eventdata, handles)
dat = get(handles.MetchGUI, 'userdata');
if (isfield(dat, 'A0') && isfield(dat, 'b0') & isfield(dat, 'node') & isfield(dat, 'elem') & ...
    isfield(dat, 'pointsinit') & isfield(dat, 'toidx') & isfield(dat, 'fromidx'))
    pmask = -1 * ones(size(dat.pointsinit, 1), 1);
    pmask(dat.fromidx) = dat.toidx;
    if (isfield(dat, 'A') && isfield(dat, 'b'))
        [Anew, bnew, posnew] = regpt2surf(dat.node, dat.elem, dat.points, pmask, dat.A, dat.b, ones(12, 1), 20);
    else
        [Anew, bnew, posnew] = regpt2surf(dat.node, dat.elem, dat.points, pmask, dat.A0, dat.b0, ones(12, 1), 20);
    end
    dat.A = Anew;
    dat.b = bnew;
    dat.pointsopt = posnew;
    if (isfield(dat, 'optplot') && dat.optplot)
        delete dat.optplot;
        dat.optplot = 0;
    end
    dat.optplot = plot3(posnew(:, 1), posnew(:, 2), posnew(:, 3), 'g+', 'parent', handles.axMesh);
else
    msgbox('You have to select 4 points and click "Initialize" button first', 'Error', 'error');
    return
end
set(handles.MetchGUI, 'userdata', dat);

% ---------------------------------------------------------------------------
function slPos_Callback(hObject, eventdata, handles)
dat = get(handles.MetchGUI, 'userdata');
if (isfield(dat, 'volume'))
    hold(handles.axMesh, 'off');
    imagesc(dat.volume(:, :, round(get(hObject, 'value'))), 'parent', handles.axMesh);
    set(handles.axMesh, 'tag', 'axMesh');
    set(handles.lbZPos, 'string', num2str(round(get(hObject, 'value'))));
end
set(handles.MetchGUI, 'userdata', dat);

% ---------------------------------------------------------------------------
function slPos_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', [.9 .9 .9]);
end

% --------------------------------------------------------------------------
function btProj_Callback(hObject, eventdata, handles)
dat = get(handles.MetchGUI, 'userdata');
if (isfield(dat, 'pointsopt'))
    if (isfield(dat, 'projplot') && dat.projplot)
        delete dat.projplot;
        dat.projplot = 0;
    end
    nv = nodesurfnorm(dat.node, dat.elem);
    [d2surf, cn] = dist2surf(dat.node, nv, dat.pointsopt);
    [dat.pointsproj dat.elemid dat.weight] = proj2mesh(dat.node, dat.elem, dat.pointsopt, nv, cn);
    hold(handles.axMesh, 'on');
    dat.projplot = plot3(dat.pointsproj(:, 1), dat.pointsproj(:, 2), dat.pointsproj(:, 3), 'c.', 'parent', handles.axMesh);
else
    msgbox('You have to first select 4 points, then click "Initialize" and "Optimize" button', 'Error', 'error');
    return
end
set(handles.MetchGUI, 'userdata', dat);

% --------------------------------------------------------------------------
function btSaveRes_Callback(hObject, eventdata, handles)
[filename, pathname] = uiputfile('*.mat', 'Save Metch Workspace as');
mapto = get(handles.lbMesh, 'userdata');
mapfrom = get(handles.lbPoints, 'userdata');

maptoidx = get(handles.txMapTo, 'userdata');
mapfromidx = get(handles.txMapFrom, 'userdata');

metchsession = get(handles.MetchGUI, 'userdata');
if (~isempty(mapto))
    metchsession.mapto = mapto;
end
if (~isempty(mapfrom))
    metchsession.mapfrom = mapfrom;
end
if (~isempty(maptoidx))
    metchsession.maptoidx = maptoidx;
end
if (~isempty(mapfromidx))
    metchsession.mapfromidx = mapfromidx;
end
fname = [pathname filename];
save(fname, 'metchsession');

% --------------------------------------------------------------------------
function btLoadSession_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.mat', 'Load Metch Workspace from');
fname = [pathname filename];
load(fname);
handle.output = metchsession;
cla(handles.axMesh);
cla(handles.axPoints);
drawinit(handles, metchsession.node, metchsession.elem, metchsession.points);

% --------------------------------------------------------------------------
function btPlotResults_Callback(hObject, eventdata, handles)

function btClose_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume;

function btHelp_Callback(hObject, eventdata, handles)
helpmsg = {
           'Metch GUI: A mesh/volume registration toolbox'
           ''
           'Author: Qianqian Fang <q.fang at neu.edu>'
           '        Martinos Center for Biomedical Imaging'
           '        Charlestown, MA 02129, USA'
           ''
           '== Description of the workflow =='
           ''
           ' 1. when the GUI pops up, it will display the mesh and the points,'
           '    you can rotate both plots so that you can identify the matching '
           '    features'
           ' 2. switch on "Select" mode, then, click on a land-mark point on the point'
           '    plot, when a data-tip shows up, click "Add Selected" button'
           ' 3. click on the corresponding position on the mesh, and click'
           '    "Add Selected"      '
           ' 4. repeat the above for at least 4 point pairs (you can select more);'
           '    if you want to change views, switch off "Select" box and rotate;'
           '    after rotation, switch on "Select" box again'
           ' 5. click "Initialize": this will create the initial mapping using the'
           '    selected point pairs'
           ' 6. click "Optimize": this will fit the surface with the whole point cloud'
           ' 7. click "Proj2Mesh": this will project the fitted point clouds onto the'
           '    mesh'
           ' 8. you can quit the GUI by hit "Close", your results will be saved to reg'
           ' 9. close the window '};

helpdlg(helpmsg);
