function varargout = img2mesh(varargin)
%
%  Format:
%      newworkspace = img2mesh or imgmesh(workspace)
%
%  A GUI for Iso2Mesh for streamlined mesh data processing
%
%  Author: Qianqian Fang <q.fang at neu.edu>
%
%  Input:
%        workspace (optional): a struct containing the below fields
%           .graph: a digraph object containing the i2m workspace data
%  Output:
%        newworkspace (optional): the updated workspace, with the same
%        subfields as the input.
%
%   If a user supplys an output variable, the GUI will not return until
%   the user closes the window; if a user does not provide any output,
%   the call will return immediately.
%
%   Please find more information at http://iso2mesh.sf.net/
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @i2m_OpeningFcn, ...
                   'gui_OutputFcn',  @i2m_OutputFcn, ...
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

% --- Executes just before i2m is made visible.
function i2m_OpeningFcn(hObject, eventdata, handles, varargin)

cm = uicontextmenu;
uimenu(cm, 'Label', 'Plot', 'CallBack', {@processdata, handles});
uimenu(cm, 'Label', 'Rename', 'CallBack', {@processdata, handles});
uimenu(cm, 'Label', 'Delete', 'CallBack', {@processdata, handles});

mimeshing = uimenu(cm, 'Label', 'Meshing');
miv2s = uimenu(mimeshing, 'Label', 'Volume to surface', 'CallBack', {@processdata, handles});
uimenu(mimeshing, 'Label', 'Volume to mesh', 'CallBack', {@processdata, handles});
uimenu(mimeshing, 'Label', 'Surface to mesh', 'CallBack', {@processdata, handles});
uimenu(mimeshing, 'Label', 'Surface to volume', 'CallBack', {@processdata, handles});
uimenu(mimeshing, 'Label', 'Close and fill volume', 'CallBack', {@processdata, handles});
uimenu(mimeshing, 'Label', 'Extract surface', 'CallBack', {@processdata, handles});
uimenu(mimeshing, 'Label', 'Tessellate surface', 'CallBack', {@processdata, handles});

mirepair = uimenu(cm, 'Label', 'Surface repair');
uimenu(mirepair, 'Label', 'Clean surface', 'CallBack', {@processdata, handles});
uimenu(mirepair, 'Label', 'Repair surface', 'CallBack', {@processdata, handles});
uimenu(mirepair, 'Label', 'Smooth surface', 'CallBack', {@processdata, handles});
uimenu(mirepair, 'Label', 'Simplify surface', 'CallBack', {@processdata, handles});
uimenu(mirepair, 'Label', 'Remesh surface', 'CallBack', {@processdata, handles});
uimenu(mirepair, 'Label', 'Simplify surface', 'CallBack', {@processdata, handles});
uimenu(mirepair, 'Label', 'Reorient mesh elements', 'CallBack', {@processdata, handles});

mibool = uimenu(cm, 'Label', 'Surface boolean');
uimenu(mibool, 'Label', 'Or', 'CallBack', {@processdata, handles});
uimenu(mibool, 'Label', 'And', 'CallBack', {@processdata, handles});
uimenu(mibool, 'Label', 'All', 'CallBack', {@processdata, handles});
uimenu(mibool, 'Label', 'Diff', 'CallBack', {@processdata, handles});
uimenu(mibool, 'Label', 'First', 'CallBack', {@processdata, handles});
uimenu(mibool, 'Label', 'Second', 'CallBack', {@processdata, handles});

mireport = uimenu(cm, 'Label', 'Report');
uimenu(mireport, 'Label', 'Containing data', 'CallBack', {@processdata, handles});
uimenu(mireport, 'Label', 'Mesh quality histogram', 'CallBack', {@processdata, handles});
uimenu(mireport, 'Label', 'Element volume histogram', 'CallBack', {@processdata, handles});
uimenu(mireport, 'Label', 'Total volume', 'CallBack', {@processdata, handles});

mirefresh = uimenu(cm, 'Label', 'Refresh', 'CallBack', {@miRefresh_Callback, handles});
uimenu(cm, 'Label', 'Save as', 'CallBack', {@processdata, handles});

miv2s.Separator = 'on';
mimeshing.Separator = 'on';
mirefresh.Separator = 'on';

root = get(handles.fgI2M, 'userdata');

if (isempty(root))
    root = struct('graph', digraph, 'menu', cm);
end
set(handles.fgI2M, 'userdata', root);
set(handles.axFlow, 'position', [0 0 1 1]);
set(handles.axPreview, 'position', [0 0 1 1]);

set(handles.fgI2M, 'UIContextMenu', handles.meCreate);

axis(handles.axFlow, 'off');
axis(handles.axPreview, 'off');

% Choose default command line output for i2m
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes i2m wait for user response (see UIRESUME)
% uiwait(handles.fgI2M);

function processdata(source, callbackdata, handles)

try
    obj = get(handles.fgI2M, 'currentobject');
    if (strcmp(class(obj), 'matlab.graphics.chart.primitive.GraphPlot') == 0)
        if (~isempty(obj.UserData) && strcmp(class(obj.UserData), 'matlab.graphics.chart.primitive.GraphPlot'))
            obj = obj.UserData;
        else
            return
        end
    end
    root = get(handles.fgI2M, 'userdata');
    pos = get(handles.axFlow, 'currentpoint');
    [nodedata, nodetype, nodeid] = getnodeat(root, obj, pos);

    newtype = dummytype;
    prefix = 'x';
    switch source.Label
        case 'Volume to surface'
            if (isstruct(nodetype) && isfield(nodetype, 'hasvol') && nodetype.hasvol)
                [newdata, newtype] = v2sgui(nodedata);
                prefix = 'Vol2Surf';
            else
                error('no volume data found');
            end
        case 'Volume to mesh'
            if (nodetype.hasvol)
                [newdata, newtype] = v2mgui(nodedata);
                prefix = 'Vol2Mesh';
            else
                error('no volume data found');
            end
        case 'Surface to mesh'
            if (nodetype.hasnode && nodetype.hasface)
                [newdata, newtype] = s2mgui(nodedata);
                prefix = 'Surf2Mesh';
            else
                error('no surface data found');
            end
        case 'Surface to volume'
            if (nodetype.hasnode && nodetype.hasface)
                ndiv = inputdlg('Division number along the shortest dimension:', ...
                                'surf2vol - rasterizing a surface mesh', 1, {'50'});
                if (isempty(ndiv))
                    return
                end
                newdata.vol = s2v(nodedata.node, nodedata.face, str2num(ndiv{1}));
                newtype.hasvol = 1;
                prefix = 'Surf2Vol';
            else
                error('no surface data found');
            end
        case 'Close and fill volume'
            if (nodetype.hasvol)
                rad = inputdlg('maximum gap length in voxel (scalar)', 'Close and fill a volume', 1, {'1'});
                newdata.vol = fillholes3d(nodedata.vol, str2num(rad{1}));
                newtype = nodetype;
                prefix = 'FillVol';
            else
                error('no volume data found');
            end
        case 'Extract surface'
            if (isstruct(nodetype) && isfield(nodetype, 'hasvol') && nodetype.hasvol)
                [newdata.node, newdata.face] = binsurface(nodedata.vol);
                prefix = 'BinSurf';
            elseif (nodetype.hasnode && nodetype.haselem)
                newdata.node = nodedata.node;
                newdata.face = volface(nodedata.elem(:, 1:min(4, size(nodedata.elem, 2))));
                [newdata.node, newdata.face] = removeisolatednode(newdata.node, newdata.face);
                prefix = 'VolSurf';
            elseif (nodetype.hasnode && nodetype.hasface)
                [newdata.node, newdata.face] = deal(nodedata.node, nodedata.face);
                prefix = 'CopySurf';
            end
            newtype.hasnode = 1;
            newtype.hasface = 1;
        case 'Clean surface'
            if (nodetype.hasnode && nodetype.hasface)
                [newdata.node, newdata.face] = meshcheckrepair(nodedata.node, nodedata.face, 'deep');
                newtype = nodetype;
                prefix = 'CleanSurf';
            else
                error('no surface data found');
            end
        case 'Repair surface'
            if (nodetype.hasnode && nodetype.hasface)
                [newdata.node, newdata.face] = meshcheckrepair(nodedata.node, nodedata.face, 'meshfix');
                newtype = nodetype;
                prefix = 'RepairSurf';
            else
                error('no surface data found');
            end
        case 'Smooth surface'
            if (nodetype.hasnode && nodetype.hasface)
                res = inputdlg({'Method (laplacian,laplacianhc,lowpass):', 'Iteration (integer):', 'Alpha (scalar):'}, ...
                               'sms - smoothing a surface mesh', [1, 1, 1], {'lowpass', '20', '0.5'});
                if (isempty(res))
                    return
                end
                newdata = nodedata;
                newtype = nodetype;
                newdata.node = sms(nodedata.node, nodedata.face, str2num(res{2}), str2num(res{3}), res{1});
                prefix = 'SmoothSurf';
            else
                error('no surface data found');
            end
        case 'Simplify surface'
            if (nodetype.hasnode && nodetype.hasface)
                res = inputdlg('Percentage of edges to keep (0-1):', ...
                               'Simplify mesh', [1], {'1'});
                if (isempty(res))
                    return
                end
                [newdata.node, newdata.face] = meshresample(nodedata.node, nodedata.face, str2num(res{1}));
                newtype = nodetype;
                prefix = 'SimplifySurf';
            else
                error('no surface data found');
            end
        case 'Remesh surface'
            if (nodetype.hasnode && nodetype.hasface)
                res = inputdlg({'Rasterization voxel size (scalar):', 'Max gap to close (scalar):', 'Max surface element radis (scalar):'}, ...
                               'Remesh surface', [1, 1, 1], {'1', '20', '3'});
                if (isempty(res))
                    return
                end
                opt.gridsize = str2num(res{1});
                opt.closesize = str2num(res{2});
                opt.elemsize = str2num(res{3});
                [newdata.node, newdata.face] = remeshsurf(nodedata.node, nodedata.face, opt);
                newtype.hasnode = 1;
                newtype.hasface = 1;
                prefix = 'RemeshSurf';
            else
                error('no surface data found');
            end
        case 'Tessellate surface'
            if (nodetype.hasnode && nodetype.hasface)
                [newdata.node, newdata.elem] = fillsurf(nodedata.node, nodedata.face);
                newdata.face = volface(newdata.elem(:, 1:min(4, size(newdata.elem, 2))));
                newtype = nodetype;
                newtype.haselem = 1;
                prefix = 'TessMesh';
            else
                error('no surface data found');
            end
        case 'Reorient mesh elements'
            if (nodetype.hasnode && nodetype.haselem)
                [newdata.node, newdata.elem] = meshreorient(nodedata.node, nodedata.elem);
                newtype = nodetype;
                prefix = 'ReorientMesh';
            elseif (nodetype.hasnode && nodetype.hasface)
                [newdata.node, newdata.face] = surfreorient(nodedata.node, nodedata.face);
                newtype = nodetype;
                prefix = 'ReorientSurf';
            else
                error('no surface or tetrahedral mesh data found');
            end
        case 'Mesh quality histogram'
            if (nodetype.hasnode && nodetype.haselem)
                quality = meshquality(nodedata.node, nodedata.elem);
            elseif (nodetype.hasnode && nodetype.hasface)
                quality = meshquality(nodedata.node, nodedata.face);
            end
            if (exist('quality', 'var'))
                figure;
                hist(quality, 50);
            end
            return
        case 'Element volume histogram'
            if (nodetype.hasnode && nodetype.haselem)
                evol = elemvolume(nodedata.node, nodedata.elem);
            elseif (nodetype.hasnode && nodetype.hasface)
                evol = elemvolume(nodedata.node, nodedata.face);
            end
            if (exist('evol', 'var'))
                figure;
                hist(evol, 50);
            end
            return
        case 'Total volume'
            if (nodetype.hasnode && nodetype.haselem)
                evol = elemvolume(nodedata.node, nodedata.elem);
            elseif (nodetype.hasnode && nodetype.hasface)
                [no, el] = fillsurf(nodedata.node, nodedata.face);
                evol = elemvolume(no, el);
            end
            if (exist('evol', 'var'))
                msgbox(sprintf('Total volume is %f cubic voxel', sum(evol)), 'Total volume');
            end
            return
        case 'Containing data'
            msg = '';
            if (nodetype.hasnode)
                msg = [msg sprintf('\nContaining %d nodes (%d columns)', size(nodedata.node, 1), size(nodedata.node, 2))];
            end
            if (nodetype.hasface)
                msg = [msg sprintf('\nContaining %d triangles (%d columns)', size(nodedata.face, 1), size(nodedata.face, 2))];
            end
            if (nodetype.haselem)
                msg = [msg sprintf('\nContaining %d tetrehedra (%d columns)', size(nodedata.elem, 1), size(nodedata.elem, 2))];
            end
            if (nodetype.hasvol)
                msg = [msg sprintf('\nContaining [%d x %d x %d ] volume', size(nodedata.vol, 1), size(nodedata.vol, 2), size(nodedata.vol, 3))];
            end
            msgbox(msg, 'Mesh data report');
            return
        case 'Plot'
            if (isstruct(nodetype) && isfield(nodetype, 'hasnode') && nodetype.hasnode)
                if (isfield(nodetype, 'haselem') && nodetype.haselem)
                    figure('keypressfcn', @plotfigevent, 'userdata', struct('node', nodedata.node, 'face', [], 'elem', nodedata.elem));
                    plotmesh(nodedata.node, [], nodedata.elem);
                else
                    figure('keypressfcn', @plotfigevent, 'userdata', struct('node', nodedata.node, 'elem', [], 'face', nodedata.face));
                    plotmesh(nodedata.node, nodedata.face);
                end
            else
                figure('keypressfcn', @plotfigevent, 'userdata', struct('vol', nodedata.vol));
                hs = slice(double(nodedata.vol), [], [ceil(size(nodedata.vol, 2) * 0.5)], ceil(size(nodedata.vol, 3) * 0.5));
                set(hs, 'linestyle', 'none');
            end
        case {'Or', 'And', 'Diff', 'All', 'First', 'Second'}
            if (isstruct(nodetype) && isfield(nodetype, 'hasnode'))
                if ((nodetype.hasface || nodetype.haselem) && nodetype.hasnode)
                    pt = ginput(1);
                    [nodedata2, nodetype2, nodeid2] = getnodeat(root, obj, pt);
                    if (~nodetype2.hasnode || ~nodetype2.hasface)
                        if (nodetype2.hasnode && nodetype2.haselem)
                            nodedata2.face = volface(nodedata2.elem);
                        else
                            error('Second operand does not contain a surface');
                        end
                    end
                    op = source.Label;
                    if (strcmp(op, 'Intersect'))
                        op = 'inter';
                    end
                    if (~nodetype.hasface)
                        nodedata.face = volface(nodedata.elem);
                    end
                    [newdata.node, newdata.face] = surfboolean(nodedata.node, nodedata.face, lower(op), nodedata2.node, nodedata2.face(:, [1 3 2]));
                    newtype.hasnode = 1;
                    newtype.hasface = 1;
                    prefix = source.Label;
                else
                    warndlg('Selected node does not contain a surface mesh');
                end
            end
        case 'Delete'
            button = questdlg('Are you sure to delete the selected node?', 'Confirm', 'No');
            if strcmpi(button, 'No') ||  strcmpi(button, 'Cancel')
                return
            end
            root.graph = rmnode(root.graph, root.graph.Nodes.Name{nodeid});
            updategraph(root, handles);
        case 'Rename'
            newname = inputdlg('Define a new name:', ...
                               'Rename', 1, {root.graph.Nodes.Name{nodeid}});
            if (isempty(newname))
                return
            end
            if (isempty(newname{1}) || ~isempty(cell2mat(regexp(root.graph.Nodes.Name, ['^' newname{1} '$']))))
                error('empty or duplicated node name');
            end
            root.graph.Nodes.Name{nodeid} = newname{1};
            updategraph(root, handles);
        case 'Save as'
            if (~nodetype.haselem && ~nodetype.hasnode)
                error('selected data does not have a mesh');
                return
            end
            filter = {'*.jmesh'; '*.*'};
            [file, path] = uiputfile(filter, 'Export mesh');
            if ~isequal(file, 0) && ~isequal(path, 0)
                if (nodetype.haselem)
                    savejmesh(nodedata.node, nodedata.face, nodedata.elem, fullfile(path, file));
                else
                    savejmesh(nodedata.node, nodedata.face, fullfile(path, file));
                end
            end
    end

    if (exist('newdata', 'var') && exist('newtype', 'var'))
        cla(handles.axPreview);
        cla(handles.axFlow);
        newdata.preview = getpreview(newdata, newtype, [400, 400]);
        [newkey, root.graph] = addnodewithdata(handles, newdata, newtype, prefix);
        root.graph = addedge(root.graph, {root.graph.Nodes.Name{nodeid}}, {newkey});
        if (strcmp(source.Parent.Type, 'uimenu') && strcmp(source.Parent.Label, 'Surface boolean'))
            root.graph = addedge(root.graph, {root.graph.Nodes.Name{nodeid2}}, {newkey});
        end
        updategraph(root, handles);
    end

catch ME
    msg = sprintf('Error: \n%s\n', ME.message);
    for e = 1:length(ME.stack)
        msg = sprintf('%s\nFile: %s\nFunction: %s\nLine: %d\n\n', msg, ME.stack(e).file, ME.stack(e).name, ME.stack(e).line);
    end
    uiwait(warndlg(msg, 'I2M ERROR'));
    updategraph(root, handles);
end

function plotfigevent(hobject, event)
data = get(hobject, 'userdata');
plotpos = [];
if (isfield(data, 'plotpos'))
    plotpos = data.plotpos;
end
if (isempty(plotpos))
    switch event.Key
        case {'rightarrow', 'uparrow'}
            plotpos = 0;
        case {'leftarrow', 'downarrow'}
            plotpos = 9;
        otherwise
    end
end
if (isfield(data, 'node'))
    pmax = max(data.node);
    pmin = min(data.node);
    cla;
    switch event.Key
        case 'rightarrow'
            plotpos = min(plotpos + 1, 9);
            plotmesh(data.node, data.face, data.elem, sprintf('x>%f', (pmax(1) - pmin(1)) * plotpos * 0.1 + pmin(1)));
        case 'leftarrow'
            plotpos = max(plotpos - 1, 0);
            plotmesh(data.node, data.face, data.elem, sprintf('x>%f', (pmax(1) - pmin(1)) * plotpos * 0.1 + pmin(1)));
        case 'uparrow'
            plotpos = min(plotpos + 1, 9);
            plotmesh(data.node, data.face, data.elem, sprintf('y>%f', (pmax(2) - pmin(2)) * plotpos * 0.1 + pmin(2)));
        case 'downarrow'
            plotpos = max(plotpos - 1, 0);
            plotmesh(data.node, data.face, data.elem, sprintf('y>%f', (pmax(2) - pmin(2)) * plotpos * 0.1 + pmin(2)));
        otherwise
    end
end
data.plotpos = plotpos;
set(hobject, 'userdata', data);

% ----------------------------------------------------------------
function [nodedata, nodetype, nodeid] = getnodeat(root, obj, pos)
nodedist = [obj.XData(:) - pos(1, 1) obj.YData(:) - pos(1, 2)];
nodedist = sum(nodedist .* nodedist, 2);
[mindist, nodeid] = min(nodedist);
nodedata = root.graph.Nodes.Data{nodeid};
nodetype = root.graph.Nodes.Type{nodeid};

% ----------------------------------------------------------------
function mytype = dummytype
mytype.hasnode = 0;
mytype.hasface = 0;
mytype.haselem = 0;
mytype.hasvol = 0;

% ----------------------------------------------------------------
function [newdata, newtype] = v2sgui(data)
prompt = {'Threshold (scalar or array):', ...
          'Surface element radius bound (scalar):', ...
          'Surface element distance bound (scalar)', 'Method: (cgalsurf,cgalmesh,simplify)'};
title = 'vol2surf - extracting surface mesh from volume';
dims = [1 1 1 1];
definput = {'0.5', '5', '1', 'cgalsurf'};
res = inputdlg(prompt, title, dims, definput);
newdata = [];
newtype = dummytype;
if (isempty(res))
    return
end

opt = struct('radbound', str2num(res{2}), 'distbound', str2num(res{3}));
[newdata.node, newdata.face] = v2s(data.vol, eval(res{1}), opt, res{4});
newtype.hasnode = 1;
newtype.hasface = 1;

% ----------------------------------------------------------------
function [newdata, newtype] = v2mgui(data)
prompt = {'Threshold (scalar or []):', ...
          'Surface element radius bound (scalar):', ...
          'Surface element distance bound (scalar)', ...
          'Max element volume (scalar):', ...
          'Method (cgalsurf,cgalmesh,simplify):'};
title = 'vol2mesh - extracting tet mesh from volume';
dims = [1 1 1 1 1];
definput = {'[]', '5', '1', '30', 'cgalmesh'};
res = inputdlg(prompt, title, dims, definput);
newdata = [];
newtype = dummytype;
if (isempty(res))
    return
end

opt = struct('radbound', str2num(res{2}), 'distbound', str2num(res{3}));
[newdata.node, newdata.elem, newdata.face] = v2m(data.vol, eval(res{1}), ...
                                                 opt, res{4});
newtype.hasnode = 1;
newtype.hasface = 1;
newtype.haselem = 1;

% ----------------------------------------------------------------
function img = getpreview(nodedata, nodetype, imsize)
hfpreview = figure('visible', 'off');
ax = axes('parent', hfpreview, 'Units', 'pixels', 'position', [1, 1, imsize(1), imsize(2)]);
if (isfield(nodetype, 'haselem') && nodetype.haselem)
    plotmesh(nodedata.node, [], nodedata.elem, 'linestyle', ':', 'edgealpha', 0.3, 'parent', ax);
elseif (isfield(nodetype, 'hasface') && nodetype.hasface)
    plotmesh(nodedata.node, nodedata.face, 'linestyle', '-', 'parent', ax);
elseif (isfield(nodetype, 'hasvol') && nodetype.hasvol)
    hs = slice(double(nodedata.vol), [], [ceil(size(nodedata.vol, 2) * 0.5)], ceil(size(nodedata.vol, 3) * 0.5), 'parent', ax);
    set(hs, 'linestyle', 'none');
elseif (isfield(nodetype, 'hasnode') && nodetype.hasnode)
    plotmesh(nodedata.node, '.', 'parent', ax);
end
set(ax, 'color', 'none');
axis(ax, 'equal');
axis(ax, 'off');
img = getframe(gca);
img = flipud(img.cdata);
delete(ax);
close(hfpreview);

% ----------------------------------------------------------------
function [newdata, newtype] = s2mgui(data)
prompt = {'Simplification ratio (%edges to keep, 0-1):', ...
          'Max element volume (scalar):', ...
          'Method (tetgen,tetgen1.5,cgalpoly):', ...
          'Region seeds (N x 3 array):', ...
          'Hole seeds (N x 3 array):'};
title = 'surf2mesh - creating tet mesh from surfaces';
dims = [1 1 1 1 1];
definput = {'1', '30', 'tetgen', '[]', '[]'};
res = inputdlg(prompt, title, dims, definput);
newdata = [];
newtype = dummytype;
if (isempty(res))
    return
end

[newdata.node, newdata.elem, newdata.face] = ...
   s2m(data.node, data.face, str2num(res{1}), ...
       str2num(res{2}), res{3}, eval(res{4}), eval(res{5}));
newtype.hasnode = 1;
newtype.hasface = 1;
newtype.haselem = 1;

% --- Outputs from this function are returned to the command line.
function varargout = i2m_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.output = get(handles.fgI2M, 'userdata');
varargout{1} = handles.output;

% --------------------------------------------------------------------
function miWeb_Callback(hObject, eventdata, handles)
web('http://iso2mesh.sourceforge.net');

% --------------------------------------------------------------------
function miDoc_Callback(hObject, eventdata, handles)
web('http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc');

% --------------------------------------------------------------------
function miAbout_Callback(hObject, eventdata, handles)
helpmsg = {
           '\bf\fontsize{12}I2M: An Integrated GUI for Iso2Mesh Meshing Toolbox\rm\fontsize{10}'
           ''
           'Copyright (c) 2018 Qianqian Fang <q.fang at neu.edu>'
           ''
           'Computational Optics&Translational Imaging Lab (http://fanglab.org)'
           'Department of Bioengineering'
           'Northeastern University'
           '360 Huntington Ave, Boston, MA 02115, USA'
           ''
           'URL:    http://iso2mesh.sourceforge.net'
           ''};

opt.Interpreter = 'tex';
opt.WindowStyle = 'modal';

msgbox(helpmsg, 'About', 'help', opt);

% --------------------------------------------------------------------
function miSphere_Callback(hObject, eventdata, handles)
prompt = {'Center:', 'Radius (scalar):', ...
          'Surface element radius bound (scalar)', ...
          'Max element volume (scalar, 0 - only create surface):'};
title = 'Create Mesh';
dims = [1 1 1 1];
definput = {'[0 0 0]', '50', '6', '20'};
res = inputdlg(prompt, title, dims, definput);
if (isempty(res))
    return
end
newtype = dummytype;
opt = str2num(res{3});
if (str2num(res{4}) == 0)
    [newdata.node, newdata.face] = meshasphere(eval(res{1}), ...
                                               str2num(res{2}), opt);
else
    [newdata.node, newdata.face, newdata.elem] = meshasphere(eval(res{1}), ...
                                                             str2num(res{2}), opt, str2num(res{4}));
    newtype.haselem = 1;
end
newtype.hasnode = 1;
newtype.hasface = 1;

if (exist('newdata', 'var') && exist('newtype', 'var'))
    newkey = addnodewithdata(handles, newdata, newtype, 'Sphere');
end

% --------------------------------------------------------------------
function miBox_Callback(hObject, eventdata, handles)
prompt = {'Diagonal end point 1 (1x3 vector):', 'Diagonal end point 2 (1x3 vector):', ...
          'Surface element radius bound (scalar)', ...
          'Max element volume (scalar, 0 - only create surface):'};
title = 'Create Mesh';
dims = [1 1 1 1];
definput = {'[0 0 0]', '[100 60 30]', '6', '30'};
res = inputdlg(prompt, title, dims, definput);
if (isempty(res))
    return
end
newtype = dummytype;
opt = str2num(res{3});
if (str2num(res{4}) == 0)
    [newdata.node, newdata.face] = meshabox(eval(res{1}), eval(res{2}), opt);
else
    [newdata.node, newdata.face, newdata.elem] = meshabox(eval(res{1}), ...
                                                          eval(res{2}), opt, str2num(res{4}));
    newtype.haselem = 1;
end
newtype.hasnode = 1;
newtype.hasface = 1;

if (exist('newdata', 'var') && exist('newtype', 'var'))
    newkey = addnodewithdata(handles, newdata, newtype, 'Box');
end

% --------------------------------------------------------------------
function miCylinder_Callback(hObject, eventdata, handles)
prompt = {'Axis end-point 1', 'Axis end-point 2', 'Radius (scalar):', ...
          'Surface element radius bound (scalar)', ...
          'Max element volume (scalar, 0 - only create surface):', ...
          'Circle division:'};
title = 'Create Mesh';
dims = [1 1 1 1 1 1];
definput = {'[0 0 0]', '[0 0 50]', '10', '3', '20', '20'};
res = inputdlg(prompt, title, dims, definput);
if (isempty(res))
    return
end
newtype = dummytype;
opt = str2num(res{4});
maxvol = str2num(res{5});

if (maxvol == 0)
    [newdata.node, newdata.face] = meshacylinder(eval(res{1}), eval(res{2}), ...
                                                 str2num(res{3}), opt);
else
    [newdata.node, newdata.face, newdata.elem] = meshacylinder(eval(res{1}), eval(res{2}), ...
                                                               str2num(res{3}), opt, maxvol);
    newtype.haselem = 1;
end
newtype.hasnode = 1;
newtype.hasface = 1;

if (exist('newdata', 'var') && exist('newtype', 'var'))
    newkey = addnodewithdata(handles, newdata, newtype, 'Cyl');
end

% --------------------------------------------------------------------
function miLoadVol_Callback(hObject, eventdata, handles)

nodedata = struct;
nodetype = dummytype;
filters = {'*.nii;*.hdr;*.img;*.tif;*.tiff;*.inr;*.bin;*.ubj', '3D volume file (*.nii;*.hdr;*.img;*.tif;*.tiff;*.inr;*.bin;*.ubj)'; ...
           '*.nii', 'Nifti file (*.nii)'; ...
           '*.hdr;*.img', 'Analyze 7.5 file (*.hdr;*.img)'; ...
           '*.tif;*.tiff', 'Multipage TIFF file (*.tif)'; ...
           '*.inr', 'INR image (*.inr)'; ...
           '*.bin', 'Binary file (*.bin)'; ...
           '*.ubj', 'Universal JSON (*.ubj)'; ...
           '*.*', 'All (*.*)'};
[file, path, idx] = uigetfile(filters);
if isequal(file, 0)
    return
else
    if (regexp(file, '\.[Nn][Ii][Ii]$'))
        im = readnifti(fullfile(path, file));
        nodedata.vol = im.img;
        nodetype.hasvol = 1;
    elseif (regexp(file, '(\.[Hh][Dd][Rr]$|\.[Ii][Mm][Gg]$)'))
        im = readnifti(fullfile(path, file));
        nodedata.vol = im.img;
        nodetype.hasvol = 1;
    elseif (regexp(file, '\.[Tt][Ii][Ff][Ff]*$'))
        nodedata.vol = readmptiff(fullfile(path, file));
        nodetype.hasvol = 1;
    elseif (regexp(file, '\.[Ii][Nn][Rr]$'))
        nodedata.vol = readinr(fullfile(path, file));
        nodetype.hasvol = 1;
    elseif (regexp(file, '\.[Bb][Ii][Nn]$'))
        prompt = {'Dimension (1x3 vector):', ...
                  'Datatype (short,float,double,integer,...):'};
        title = 'Load generic binary file';
        dims = [1 1];
        definput = {'[]', 'short'};
        [res, isok] = inputdlg(prompt, title, dims, definput);
        if (isok == 0)
            return
        end
        nodedata.vol = loadmc2(fullfile(path, file), eval(res{1}), res{2});
        nodetype.hasvol = 1;
    elseif (regexp(file, '\.[Uu][Bb][Jj]$'))
        nodedata = loadbj(fullfile(path, file));
        if (isstruct(nodedata) && isfield(nodedata, 'vol'))
            nodetype.hasvol = 1;
        end
    end
end

if (exist('nodedata', 'var'))
    nodetype = getnodetype(nodedata);
    if (nodetype.hasvol)
        addnodewithdata(handles, nodedata, nodetype, 'Vol');
    end
else
    warndlg('no valid mesh data found', 'Warning');
end

% --------------------------------------------------------------------
function miLoadMesh_Callback(hObject, eventdata, handles)

nodedata = struct;
nodetype = dummytype;
filters = {'*.jmesh;*.off;*.medit;*.smf;*.json', '3D Mesh files (*.jmesh;*.off;*.medit;*.smf;*.json)'; ...
           '*.jmesh', 'JSON mesh (*.jmesh)'; ...
           '*.off', 'OFF file (*.off)'; ...
           '*.medit', 'Medit file (*.medit)'; ...
           '*.ele', 'Tetgen element mesh file (*.ele)'; ...
           '*.json', 'JSON file (*.json)'; '*.*', 'All (*.*)'};
[file, path, idx] = uigetfile(filters);
if isequal(file, 0)
    return
else
    if (regexp(file, '\.[Oo][Ff][Ff]$'))
        [nodedata.node, nodedata.face] = readoff(fullfile(path, file));
    elseif (regexp(file, '\.[Mm][Ee][Dd][Ii][Tt]$'))
        [nodedata.node, nodedata.elem] = readmedit(fullfile(path, file));
    elseif (regexp(file, '\.[Ee][Ll][Ee]$'))
        [pathstr, name, ext] = fileparts(fullfile(path, file));
        [nodedata.node, nodedata.elem] = readtetgen(fullfile(pathstr, name));
    elseif (regexp(file, '\.[Jj][Mm][Ee][Ss][Hh]$'))
        nodedata = importjmesh(fullfile(path, file));
    elseif (regexp(file, '\.[Jj][Ss][Oo][Nn]$'))
        nodedata = loadjson(fullfile(path, file));
    end
end
if (exist('nodedata', 'var'))
    adddatatograph(handles, nodedata);
else
    warndlg('no valid mesh data found', 'Warning');
end

% --------------------------------------------------------------------
function miLoadSurf_Callback(hObject, eventdata, handles)
nodedata = struct;
nodetype = dummytype;
filters = {'*.jmesh;*.off;*.asc;*.smf;*.smf;*.json', '3D Mesh files (*.jmesh;*.off;*.asc;*.smf;*.smf;*.json)'; ...
           '*.jmesh', 'JSON mesh (*.jmesh)'; ...
           '*.off', 'OFF file (*.off)'; ...
           '*.asc', 'ASC file (*.asc)'; ...
           '*.gts', 'GNU Trangulated Surface file (*.gts)'; ...
           '*.smf', 'Simple Model Format (*.smf)'; ...
           '*.json', 'JSON file (*.json)'; '*.*', 'All (*.*)'};
[file, path, idx] = uigetfile(filters);
if isequal(file, 0)
    return
else
    if (regexp(file, '\.[Oo][Ff][Ff]$'))
        [nodedata.node, nodedata.face] = readoff(fullfile(path, file));
    elseif (regexp(file, '\.[Aa][Ss][Cc]$'))
        [nodedata.node, nodedata.face] = readasc(fullfile(path, file));
    elseif (regexp(file, '\.[Gg][Tt][Ss]$'))
        [nodedata.node, nodedata.face] = readgts(fullfile(path, file));
    elseif (regexp(file, '\.[Ss][Mm][Ff]$'))
        [nodedata.node, nodedata.face] = readsmf(fullfile(path, file));
    elseif (regexp(file, '\.[Jj][Mm][Ee][Ss][Hh]$'))
        nodedata = importjmesh(fullfile(path, file));
    elseif (regexp(file, '\.[Jj][Ss][Oo][Nn]$'))
        nodedata = loadjson(fullfile(path, file));
    end
end
if (exist('nodedata', 'var'))
    adddatatograph(handles, nodedata);
else
    warndlg('no valid mesh data found', 'Warning');
end

% --------------------------------------------------------------------
function nodedata = importjmesh(filename)
data = loadjson(filename);
nodedata = struct;
if (isfield(data, 'MeshNode'))
    nodedata.node = data.MeshNode;
end
if (isfield(data, 'MeshElem'))
    nodedata.elem = data.MeshElem;
end
if (isfield(data, 'MeshSurf'))
    nodedata.face = data.MeshSurf;
end
if (isfield(data, 'MeshNodeVal'))
    nodedata.node(:, end + 1:end + size(data.MeshNodeVal, 2)) = data.MeshNodeVal;
end
if (isfield(data, 'MeshTetraVal'))
    nodedata.elem(:, end + 1:end + size(data.MeshTetraVal, 2)) = data.MeshTetraVal;
end

% --------------------------------------------------------------------
function adddatatograph(handles, nodedata)
nodetype = getnodetype(nodedata);
if (nodetype.haselem)
    addnodewithdata(handles, nodedata, nodetype, 'Tet');
elseif (nodetype.hasface)
    addnodewithdata(handles, nodedata, nodetype, 'Surf');
elseif (nodetype.hasnode)
    addnodewithdata(handles, nodedata, nodetype, 'Point');
elseif (nodetype.hasvol)
    addnodewithdata(handles, nodedata, nodetype, 'Vol');
end

% --- Executes during object creation, after setting all properties.
function axFlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axFlow

% --- Executes during object creation, after setting all properties.
function fgI2M_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fgI2M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --------------------------------------------------------------------
function nodetype = getnodetype(nodedata)
nodetype = dummytype;
if (~isstruct(nodedata))
    return
end
names = fieldnames(nodedata);
for i = 1:length(names)
    switch names{i}
        case 'node'
            nodetype.hasnode = 1;
        case 'face'
            nodetype.hasface = 1;
        case 'elem'
            nodetype.haselem = 1;
        case 'vol'
            nodetype.hasvol = 1;
    end
    dat = nodedata.(names{i});
    if ((isnumeric(dat) || islogical(dat)) && ndims(dat) == 3)
        nodetype.hasvol = 1;
    end
end
% --------------------------------------------------------------------

function [key, newgraph] = addnodewithdata(handles, nodedata, nodetype, name)

root = get(handles.fgI2M, 'userdata');

if (isempty(root))
    root = struct('graph', digraph, 'menu', uicontextmenu);
end
if (nargin < 4)
    name = 'x';
end

id = 1;
if (~isempty(root.graph.Nodes))
    while (find(strcmp(root.graph.Nodes.Name, sprintf('%s%d', name, id))))
        id = id + 1;
    end
end
key = sprintf('%s%d', name, id);

cla(handles.axPreview);
cla(handles.axFlow);
nodedata.preview = getpreview(nodedata, nodetype, [400, 400]);

nodeprop = table({key}, {nodedata}, {nodetype}, 'VariableNames', {'Name', 'Data', 'Type'});
root.graph = addnode(root.graph, nodeprop);
if (nargout > 1)
    newgraph = root.graph;
end
updategraph(root, handles);

% --------------------------------------------------------------------
function hobj = updatepreview(root, obj, handles)
cla(handles.axPreview);
view(handles.axPreview, 2);
nx = obj.XData(:);
ny = obj.YData(:);
nn = length(nx);
hold(handles.axPreview, 'on');
dx = get(handles.axFlow, 'xlim');
dy = get(handles.axFlow, 'ylim');
wd = min([diff(dx) diff(dy)]) / 15;
hobj = zeros(1, nn);
set(handles.axPreview, 'xlim', dx);
set(handles.axPreview, 'ylim', dy);
dim = get(handles.axPreview, 'dataaspectratio');
[wfig, hfig] = getwindowsize(handles.fgI2M);
wfig = wfig * dim(2);
hfig = hfig * dim(1);

for i = 1:nn
    if (isfield(root.graph.Nodes.Data{i}, 'preview'))
        hobj(i) = imagesc(nx(i) + [-2 * wd 0], ny(i) + [-wd wd] * (wfig / hfig), imresize(root.graph.Nodes.Data{i}.preview, 0.25), ...
                          'parent', handles.axPreview);
    end
end
set(hobj, 'userdata', obj);
set(hobj, 'UIContextMenu', root.menu);
set(handles.axPreview, 'xlim', dx);
set(handles.axPreview, 'ylim', dy);
axis(handles.axPreview, 'off');
hold(handles.axPreview, 'off');
% --------------------------------------------------------------------

function [width, height] = getwindowsize(fig)
oldunits = get(fig, 'Units');
set(fig, 'Units', 'pixels');
figpos = get(fig, 'Position');
set(fig, 'Units', oldunits);
width = figpos(3);
height = figpos(4);

function [hg, hobj] = updategraph(root, handles)
set(handles.fgI2M, 'userdata', root);
hg = plot(root.graph, 'parent', handles.axFlow, 'ArrowSize', 15);
hobj = updatepreview(root, hg, handles);
% set(hg,'Selected','on');
axis(handles.axFlow, 'off');
set(handles.axFlow, 'XColor', 'none');
set(hg, 'UIContextMenu', root.menu);

% --------------------------------------------------------------------
function miEllipsoid_Callback(hObject, eventdata, handles)
prompt = {'Center (1x3 vector):', ...
          'Radii (a scalar, or 1x3 or 1x5 vector):', ...
          'Max element volume (scalar, 0 - only create surface):', ...
          'Circle division:'};
title = 'Create an Ellipsoid Mesh';
dims = [1 1 1 1];
definput = {'[0 0 0]', '[50 30 20]', '3', '30'};
res = inputdlg(prompt, title, dims, definput);
if (isempty(res))
    return
end
newtype = dummytype;
opt = str2num(res{3});
maxvol = str2num(res{4});

if (maxvol == 0)
    [newdata.node, newdata.face] = meshanellip(eval(res{1}), eval(res{2}), opt);
else
    [newdata.node, newdata.face, newdata.elem] = meshanellip(eval(res{1}), ...
                                                             eval(res{2}), opt, maxvol);
    newtype.haselem = 1;
end
newtype.hasnode = 1;
newtype.hasface = 1;

if (exist('newdata', 'var') && exist('newtype', 'var'))
    newkey = addnodewithdata(handles, newdata, newtype, 'Cyl');
end

% --------------------------------------------------------------------
function miLattice_Callback(hObject, eventdata, handles)
prompt = {'X-lattice range (a vector):', ...
          'Y-lattice range (a vector):', ...
          'Z-lattice range (a vector):', ...
          'Max element volume (scalar, 0 - only create surface):'};
title = 'Create Lattice Grid Mesh';
dims = [1 1 1 1];
definput = {'[1 100]', '[1 50]', '[1 10 30]', '30'};
res = inputdlg(prompt, title, dims, definput);
if (isempty(res))
    return
end
newtype = dummytype;
maxvol = str2num(res{4});
if (maxvol == 0)
    [newdata.node, newdata.face] = latticegrid(eval(res{1}), eval(res{2}), eval(res{3}));
else
    [no, fc, c0] = latticegrid(eval(res{1}), eval(res{2}), eval(res{3}));
    [newdata.node, newdata.elem, newdata.face] = surf2mesh(no, fc, [], [], 1, maxvol, c0);
    newtype.haselem = 1;
end
newtype.hasnode = 1;
newtype.hasface = 1;

if (exist('newdata', 'var') && exist('newtype', 'var'))
    newkey = addnodewithdata(handles, newdata, newtype, 'Lattice');
end

% --------------------------------------------------------------------
function miMeshgrid5_Callback(hObject, eventdata, handles)
prompt = {'X-lattice range (a vector):', ...
          'Y-lattice range (a vector):', ...
          'Z-lattice range (a vector):'};
title = 'Create Meshgrid (5 tet/cell) Mesh';
dims = [1 1 1];
definput = {'1:10', '1:8', '1:5'};
res = inputdlg(prompt, title, dims, definput);
if (isempty(res))
    return
end
newtype = dummytype;

[newdata.node, newdata.elem] = meshgrid5(eval(res{1}), eval(res{2}), eval(res{3}));
newdata.face = volface(newdata.elem);
newtype.hasnode = 1;
newtype.hasface = 1;
newtype.haselem = 1;

if (exist('newdata', 'var') && exist('newtype', 'var'))
    newkey = addnodewithdata(handles, newdata, newtype, 'Meshgrid5_');
end

% --------------------------------------------------------------------
function miMeshgrid6_Callback(hObject, eventdata, handles)
prompt = {'X-lattice range (a vector):', ...
          'Y-lattice range (a vector):', ...
          'Z-lattice range (a vector):'};
title = 'Create Meshgrid (6 tet/cell) Mesh';
dims = [1 1 1];
definput = {'1:10', '1:8', '1:5'};
res = inputdlg(prompt, title, dims, definput);
if (isempty(res))
    return
end
newtype = dummytype;

[newdata.node, newdata.elem] = meshgrid5(eval(res{1}), eval(res{2}), eval(res{3}));
newdata.face = volface(newdata.elem);
newtype.hasnode = 1;
newtype.hasface = 1;
newtype.haselem = 1;

if (exist('newdata', 'var') && exist('newtype', 'var'))
    newkey = addnodewithdata(handles, newdata, newtype, 'Meshgrid6_');
end

% --------------------------------------------------------------------
function miOpen_Callback(hObject, eventdata, handles)
filter = {'*.mat'; '*.*'};
root = get(handles.fgI2M, 'userdata');

[file, path] = uigetfile(filter, 'Load workspace');
if isequal(file, 0)
    return
else
    data = load(fullfile(path, file));
end
if (isfield(data, 'i2mworkspace'))
    root.graph = data.i2mworkspace;
    updategraph(root, handles);
else
    warndlg('no saved workspace found', 'Warning');
end

% --------------------------------------------------------------------
function miSaveAll_Callback(hObject, eventdata, handles)
root = get(handles.fgI2M, 'userdata');
i2mworkspace = root.graph;
filter = {'*.mat'; '*.*'};
[file, path] = uiputfile(filter, 'Save workspace');
if ~isequal(file, 0) && ~isequal(path, 0)
    save(fullfile(path, file), 'i2mworkspace');
end

% --------------------------------------------------------------------
function miLoadVar_Callback(hObject, eventdata, handles)

nodedata = struct;
[file, path] = uigetfile({'*.mat', 'MATLAB data (*.mat)'});
if isequal(file, 0)
    return
else
    data = load(fullfile(path, file));
    vars = fieldnames(data);
    [idx, isok] = listdlg('ListString', vars, ...
                          'PromptString', 'Select a 3D array:');
    if (~isok)
        return
    end
    for i = 1:length(idx)
        dat = data.(vars{idx(i)});
        if (ndims(dat) == 3 && ~isfield(nodedata, 'vol'))
            nodedata.vol = dat;
        else
            nodedata.(vars{idx(i)}) = dat;
        end
    end
end
adddatatograph(handles, nodedata);

% --------------------------------------------------------------------
function miExit_Callback(hObject, eventdata, handles)
close(handles.fgI2M);

% --------------------------------------------------------------------
function miRefresh_Callback(hObject, eventdata, handles)
root = get(handles.fgI2M, 'userdata');
updategraph(root, handles);
