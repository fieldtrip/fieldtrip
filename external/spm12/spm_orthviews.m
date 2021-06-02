function varargout = spm_orthviews(action,varargin)
% Display orthogonal views of a set of images
% FORMAT H = spm_orthviews('Image',filename[,area[,F]])
% filename - name of image to display
% area     - position of image {relative}
%            [left, bottom, width, height]
% F        - figure handle
% H        - handle for orthogonal sections
%
% FORMAT spm_orthviews('Reposition',centre)
% centre   - X, Y & Z coordinates of centre voxel {mm}
%
% FORMAT spm_orthviews('Space'[,handle[,M,dim]])
% handle   - the view to define the space by, optionally with extra
%            transformation matrix and dimensions (e.g. one of the blobs
%            of a view)
% with no arguments - puts things into mm space
%
% FORMAT H = spm_orthviews('Caption', handle, string, [Property, Value])
% handle   - the view to which a caption should be added
% string   - the caption text to add
% Property,Value - optional, e.g. 'FontWeight', 'Bold'
% H        - the handle to the object whose String property has the caption
%
% FORMAT spm_orthviews('BB',bb)
% bb       - bounding box
%            [loX loY loZ
%             hiX hiY hiZ]
%
% FORMAT spm_orthviews('MaxBB')
% Set the bounding box big enough to display the whole of all images.
%
% FORMAT spm_orthviews('Resolution'[,res])
% res      - resolution (mm)
% Set the sampling resolution for all images. The effective resolution will
% be the minimum of res and the voxel sizes of all images. If no resolution
% is specified, the minimum of 1mm and the voxel sizes of the images is
% used.
%
% FORMAT spm_orthviews('Zoom'[,fov[,res]])
% fov      - half width of field of view (mm)
% res      - resolution (mm)
% Set the displayed part and sampling resolution for all images. The image
% display will be centered at the current crosshair position. The image
% region [xhairs-fov xhairs+fov] will be shown.
% If no argument is given or fov == Inf, the image display will be reset to
% "Full Volume". If fov == 0, the image will be zoomed to the bounding box
% from spm_get_bbox for the non-zero voxels of the image. If fov is NaN,
% then a threshold can be entered, and spm_get_bbox will be used to derive
% the bounding box of the voxels above this threshold.
% Optionally, the display resolution can be set as well.
%
% FORMAT spm_orthviews('Redraw')
% Redraw the images.
%
% FORMAT spm_orthviews('Reload_mats')
% Reload the voxel-world mapping matrices from the headers stored on disk,
% e.g. following reorientation of some images.
%
% FORMAT spm_orthviews('Delete', handle)
% handle   - image number to delete
%
% FORMAT spm_orthviews('Reset')
% Clear the orthogonal views
%
% FORMAT spm_orthviews('Pos')
% Return the co-ordinate of the crosshairs in millimetres in the standard
% space.
%
% FORMAT spm_orthviews('Pos', i)
% Return the voxel co-ordinate of the crosshairs in the image in the ith
% orthogonal section.
%
% FORMAT spm_orthviews('Xhairs','off') OR spm_orthviews('Xhairs')
% Disable the cross-hairs on the display.
%
% FORMAT spm_orthviews('Xhairs','on')
% Enable the cross-hairs.
%
% FORMAT spm_orthviews('Interp',hld)
% Set the hold value to hld (see spm_slice_vol).
%
% FORMAT spm_orthviews('AddBlobs',handle,XYZ,Z,mat,name)
% Add blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations
% Z        - blob voxel intensities
% mat      - matrix from voxels to millimetres of blob
% name     - a name for this blob
% This method only adds one set of blobs, and displays them using a split
% colour table.
%
% FORMAT spm_orthviews('SetBlobsMax', vn, bn, mx)
% Set maximum value for blobs overlay number bn of view number vn to mx.
%
% FORMAT spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour,name)
% Add blobs from a pointlist to the image specified by the handle(s)
% handle   - image number to add blobs to
% XYZ      - blob voxel locations
% Z        - blob voxel intensities
% mat      - matrix from voxels to millimeters of blob.
% colour   - the 3 vector containing the colour that the blobs should be
% name     - a name for this blob
% Several sets of blobs can be added in this way, and it uses full colour.
% Although it may not be particularly attractive on the screen, the colour
% blobs print well.
%
% FORMAT spm_orthviews('AddColourBar',handle,blobno)
% Add colourbar for a specified blob set.
% handle    - image number
% blobno    - blob number
%
% FORMAT spm_orthviews('RemoveBlobs',handle)
% Remove all blobs from the image specified by the handle(s).
%
% FORMAT spm_orthviews('Addtruecolourimage',handle,filename,colourmap,prop,mx,mn)
% Add blobs from an image in true colour.
% handle    - image number to add blobs to [Default: 1]
% filename  - image containing blob data [Default: GUI input]
% colourmap - colormap to display blobs in [Default: GUI input]
% prop      - intensity proportion of activation cf grayscale [default: 0.4]
% mx        - maximum intensity to scale to [maximum value in activation image]
% mn        - minimum intensity to scale to [minimum value in activation image]
%
% FORMAT spm_orthviews('Register',hReg)
% hReg      - Handle of HandleGraphics object to build registry in
% See spm_XYZreg for more information.
%
% FORMAT spm_orthviews('AddContext',handle)
% FORMAT spm_orthviews('RemoveContext',handle)
% handle    - image number to add/remove context menu to
%
% FORMAT spm_orthviews('ZoomMenu',zoom,res)
% FORMAT [zoom, res] = spm_orthviews('ZoomMenu')
% zoom      - A list of predefined zoom values
% res       - A list of predefined resolutions
% This list is used by spm_image and spm_orthviews('addcontext',...) to
% create the 'Zoom' menu. The values can be retrieved by calling
% spm_orthviews('ZoomMenu') with 2 output arguments. Values of 0, NaN and
% Inf are treated specially, see the help for spm_orthviews('Zoom' ...).
%__________________________________________________________________________
%
% PLUGINS
% The display capabilities of spm_orthviews can be extended with plugins.
% These are located in the spm_orthviews subdirectory of the SPM
% distribution.
% The functionality of plugins can be accessed via calls to
% spm_orthviews('plugin_name', plugin_arguments). For detailed descriptions
% of each plugin see help spm_orthviews/spm_ov_'plugin_name'.
%__________________________________________________________________________
% Copyright (C) 1996-2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner et al
% $Id: spm_orthviews.m 7376 2018-07-20 10:30:59Z guillaume $


% The basic fields of st are:
%         n        - the number of images currently being displayed.
%         vols     - a cell array containing the data on each of the
%                    displayed images.
%         Space    - a mapping between the displayed images and the
%                    mm space of each image.
%         bb       - the bounding box of the displayed images.
%         centre   - the current centre of the orthogonal views.
%         callback - a callback to be evaluated on a button-click.
%         xhairs   - crosshairs off/on.
%         hld      - the interpolation method.
%         fig      - the figure that everything is displayed in.
%         mode     - the position/orientation of the sagittal view.
%                    (currently always 1)
%
%         st.registry.hReg \_ See spm_XYZreg for documentation
%         st.registry.hMe  /
%
% For each of the displayed images, there is a non-empty entry in the
% vols cell array.  Handles returned by "spm_orthviews('Image',.....)"
% indicate the position in the cell array of the newly created ortho-view.
% Operations on each ortho-view require the handle to be passed.
%
% When a new image is displayed, the cell entry contains the information
% returned by spm_vol (type help spm_vol for more info).  In addition,
% there are a few other fields, some of which are documented here:
%
%         premul  - a matrix to premultiply the .mat field by.  Useful
%                   for re-orienting images.
%         window  - either 'auto' or an intensity range to display the
%                   image with.
%         mapping - Mapping of image intensities to grey values. Currently
%                   one of 'linear', 'histeq', loghisteq',
%                   'quadhisteq'. Default is 'linear'.
%                   Histogram equalisation depends on the image toolbox
%                   and is only available if there is a license available
%                   for it.
%         ax      - a cell array containing an element for the three
%                   views.  The fields of each element are handles for
%                   the axis, image and crosshairs.
%
%         blobs   - optional.  Is there for using to superimpose blobs.
%                   vol     - 3D array of image data
%                   mat     - a mapping from vox-to-mm (see spm_vol, or
%                             help on image formats).
%                   max     - maximum intensity for scaling to.  If it
%                             does not exist, then images are auto-scaled.
%
%                   There are two colouring modes: full colour, and split
%                   colour.  When using full colour, there should be a
%                   'colour' field for each cell element.  When using
%                   split colourscale, there is a handle for the colorbar
%                   axis.
%
%                   colour  - if it exists it contains the
%                             red,green,blue that the blobs should be
%                             displayed in.
%                   cbar    - handle for colorbar (for split colourscale).
%
% PLUGINS
% The plugin concept has been developed to extend the display capabilities
% of spm_orthviews without the need to rewrite parts of it. Interaction
% between spm_orthviews and plugins takes place
% a) at startup: The subfunction 'reset_st' looks for folders
%                'spm_orthviews' in spm('Dir') and each toolbox
%                folder. Files with a name spm_ov_PLUGINNAME.m in any of
%                these folders will be treated as plugins.
%                For each such file, PLUGINNAME will be added to the list
%                st.plugins{:}.
%                The subfunction 'add_context' calls each plugin with
%                feval(['spm_ov_', st.plugins{k}], ...
%                  'context_menu', i, parent_menu)
%                Each plugin may add its own submenu to the context
%                menu.
% b) at redraw:  After images and blobs of st.vols{i} are drawn, the
%                struct st.vols{i} is checked for field names that occur in
%                the plugin list st.plugins{:}. For each matching entry, the
%                corresponding plugin is called with the command 'redraw':
%                feval(['spm_ov_', st.plugins{k}], ...
%                  'redraw', i, TM0, TD, CM0, CD, SM0, SD);
%                The values of TM0, TD, CM0, CD, SM0, SD are defined in the
%                same way as in the redraw subfunction of spm_orthviews.
%                It is up to the plugin to do all necessary redraw
%                operations for its display contents. Each displayed item
%                must have set its property 'HitTest' to 'off' to let events
%                go through to the underlying axis, which is responsible for
%                callback handling. The order in which plugins are called is
%                undefined.


global st

persistent zoomlist reslist

if isempty(st), reset_st; end

if ~nargin, action = ''; end

switch lower(action)
    case 'image'
        if numel(varargin) > 2
            F = varargin{3};
        else
            F = st.fig;
        end
        H = specify_image(varargin{1},F);
        if ~isempty(H)
            if numel(varargin)>=2 && ~isempty(varargin{2})
                st.vols{H}.area = varargin{2};
            else
                st.vols{H}.area = [0 0 1 1];
            end
            if isempty(st.bb), st.bb = maxbb; end
            resolution;
            bbox(H);
        end
        varargout{1} = H;
        mmcentre     = mean(st.Space*[st.bb';1 1],2)';
        st.centre    = mmcentre(1:3);
        redraw_all;
        try, set(F, 'WindowScrollWheelFcn', @scroll_wheel); end


    case 'caption'
        if ~isnumeric(varargin{1})
            varargin{1} = cellstr(varargin{1});
            xlh = NaN(numel(varargin{1}),1);
            for i=1:numel(varargin{1})
                h = spm_orthviews('Caption',i,varargin{1}{i},varargin{3:end});
                if ~isempty(h), xlh(i) = h; end
            end
            varargout{1} = xlh;
            return;
        end
        
        vh = valid_handles(varargin{1});
        nh = numel(vh);
        
        xlh = nan(nh, 1);
        for i = 1:nh
            xlh(i) = get(st.vols{vh(i)}.ax{3}.ax, 'XLabel');
            if iscell(varargin{2})
                if i <= length(varargin{2})
                    set(xlh(i), 'String', varargin{2}{i});
                end
            else
                set(xlh(i), 'String', varargin{2});
            end
            for np = 4:2:nargin
                property = varargin{np-1};
                value = varargin{np};
                set(xlh(i), property, value);
            end
        end
        varargout{1} = xlh;
        
    case 'bb'
        if ~isempty(varargin) && all(size(varargin{1})==[2 3]), st.bb = varargin{1}; end
        bbox;
        redraw_all;
        
    case 'redraw'
        if isempty(varargin)
            redraw_all;
        else
            redraw(varargin{1});
        end
        callback;
        if isfield(st,'registry')
            spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
        end
        
    case 'reload_mats'
        if nargin > 1
            handles = valid_handles(varargin{1});
        else
            handles = valid_handles;
        end
        for i = handles
            fnm = spm_file(st.vols{i}.fname, 'number', st.vols{i}.n);
            st.vols{i}.mat = spm_get_space(fnm);
        end
        % redraw_all (done in spm_orthviews('reorient','context_quit'))
        
    case 'reposition'
        if isempty(varargin), tmp = findcent;
        else tmp = varargin{1}; end
        if numel(tmp) == 3
            h = valid_handles(st.snap);
            if ~isempty(h)
                tmp = st.vols{h(1)}.mat * ...
                    round(st.vols{h(1)}.mat\[tmp(:); 1]);
            end
            st.centre = tmp(1:3);
        end
        redraw_all;
        callback;
        if isfield(st,'registry')
            spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
        end
        
    case 'setcoords'
        st.centre = varargin{1};
        st.centre = st.centre(:);
        redraw_all;
        callback;
        
    case 'space'
        if isempty(varargin)
            st.Space = eye(4);
            st.bb = maxbb;
        else
            space(varargin{:});
        end
        resolution;
        bbox;
        redraw_all;
        
    case 'maxbb'
        st.bb = maxbb;
        bbox;
        redraw_all;
        
    case 'resolution'
        resolution(varargin{:});
        bbox;
        redraw_all;
        
    case 'window'
        if numel(varargin)<2
            win = 'auto';
        elseif numel(varargin{2})==2
            win = varargin{2};
        end
        for i=valid_handles(varargin{1})
            st.vols{i}.window = win;
        end
        redraw(varargin{1});
        
    case 'delete'
        my_delete(varargin{1});
        
    case 'move'
        move(varargin{1},varargin{2});
        % redraw_all;
        
    case 'reset'
        my_reset;
        
    case 'pos'
        if isempty(varargin)
            H = st.centre(:);
        else
            H = pos(varargin{1});
        end
        varargout{1} = H;
        
    case 'interp'
        st.hld = varargin{1};
        redraw_all;
        
    case 'xhairs'
        xhairs(varargin{:});
        
    case 'register'
        register(varargin{1});
        
    case 'addblobs'
        addblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'setblobsmax'
        st.vols{varargin{1}}.blobs{varargin{2}}.max = varargin{3};
        spm_orthviews('redraw')
        
    case 'addcolouredblobs'
        addcolouredblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'addimage'
        addimage(varargin{1}, varargin{2});
        % redraw(varargin{1});
        
    case 'addcolouredimage'
        addcolouredimage(varargin{1}, varargin{2},varargin{3});
        % redraw(varargin{1});
        
    case 'addtruecolourimage'
        if nargin < 2
            varargin(1) = {1};
        end
        if nargin < 3
            varargin(2) = {spm_select(1, 'image', 'Image with activation signal')};
        end
        if nargin < 4
            actc = [];
            while isempty(actc)
                actc = getcmap(spm_input('Colourmap for activation image', '+1','s'));
            end
            varargin(3) = {actc};
        end
        if nargin < 5
            varargin(4) = {0.4};
        end
        if nargin < 6
            actv = spm_vol(varargin{2});
            varargin(5) = {max([eps maxval(actv)])};
        end
        if nargin < 7
            varargin(6) = {min([0 minval(actv)])};
        end
        
        addtruecolourimage(varargin{1}, varargin{2},varargin{3}, varargin{4}, ...
            varargin{5}, varargin{6});
        % redraw(varargin{1});
        
    case 'addcolourbar'
        addcolourbar(varargin{1}, varargin{2});
        
    case {'removeblobs','rmblobs'}
        rmblobs(varargin{1});
        redraw(varargin{1});
        
    case 'addcontext'
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        addcontexts(handles);
        
    case {'removecontext','rmcontext'}
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        rmcontexts(handles);
        
    case 'context_menu'
        c_menu(varargin{:});
        
    case 'valid_handles'
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        varargout{1} = valid_handles(handles);

    case 'zoom'
        zoom_op(varargin{:});
        
    case 'zoommenu'
        if isempty(zoomlist)
            zoomlist = [NaN 0 5    10  20 40 80 Inf];
            reslist  = [1   1 .125 .25 .5 .5 1  1  ];
        end
        if nargin >= 3
            if all(cellfun(@isnumeric,varargin(1:2))) && ...
                    numel(varargin{1})==numel(varargin{2})
                zoomlist = varargin{1}(:);
                reslist  = varargin{2}(:);
            else
                warning('spm_orthviews:zoom',...
                        'Invalid zoom or resolution list.')
            end
        end
        if nargout > 0
            varargout{1} = zoomlist;
        end
        if nargout > 1
            varargout{2} = reslist;
        end
        
    otherwise
        addonaction = strcmpi(st.plugins,action);
        if any(addonaction)
            feval(['spm_ov_' st.plugins{addonaction}],varargin{:});
        end
end


%==========================================================================
% function H = specify_image(img, h)
%==========================================================================
function H = specify_image(img,h)
global st
H = [];
if isstruct(img)
    V = img(1);
else
    try
        V = spm_vol(img);
    catch
        fprintf('Can not use image "%s"\n', img);
        return;
    end
end
if numel(V)>1, V=V(1); end

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end

DeleteFcn = ['spm_orthviews(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3
    ax = axes('Visible','off', 'Parent',h, ...
        'YDir','normal', 'DeleteFcn',DeleteFcn, 'ButtonDownFcn',@repos_start);
    d  = image(0, 'Tag','Transverse', 'Parent',ax, 'DeleteFcn',DeleteFcn);
    set(ax, 'Ydir','normal', 'ButtonDownFcn',@repos_start);
    
    lx = line(0,0, 'Parent',ax, 'DeleteFcn',DeleteFcn, 'Color',[0 0 1]);
    ly = line(0,0, 'Parent',ax, 'DeleteFcn',DeleteFcn, 'Color',[0 0 1]);
    if ~st.xhairs
        set(lx, 'Visible','off');
        set(ly, 'Visible','off');
    end
    V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end
V.premul    = eye(4);
V.window    = 'auto';
V.mapping   = 'linear';
st.vols{ii} = V;

H = ii;


%==========================================================================
% function addblobs(handle, xyz, t, mat, name)
%==========================================================================
function addblobs(handle, xyz, t, mat, name)
global st
if nargin < 5
    name = '';
end
for i=valid_handles(handle)
    if ~isempty(xyz)
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        st.vols{i}.blobs=cell(1,1);
        mx = max([eps max(t)]);
        mn = min([0 min(t)]);
        st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx, 'min',mn,'name',name);
        addcolourbar(handle,1);
    end
end


%==========================================================================
% function addimage(handle, fname)
%==========================================================================
function addimage(handle, fname)
global st
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    st.vols{i}.blobs=cell(1,1);
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx,'min',mn);
    addcolourbar(handle,1);
end


%==========================================================================
% function addcolouredblobs(handle, xyz, t, mat, colour, name)
%==========================================================================
function addcolouredblobs(handle, xyz, t, mat, colour, name)
if nargin < 6
    name = '';
end
global st
for i=valid_handles(handle)
    if ~isempty(xyz)
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        if ~isfield(st.vols{i},'blobs')
            st.vols{i}.blobs=cell(1,1);
            bset = 1;
        else
            bset = numel(st.vols{i}.blobs)+1;
        end
        mx = max([eps maxval(vol)]);
        mn = min([0 minval(vol)]);
        st.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
            'max',mx, 'min',mn, 'colour',colour, 'name',name);
    end
end


%==========================================================================
% function addcolouredimage(handle, fname,colour)
%==========================================================================
function addcolouredimage(handle, fname,colour)
global st
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    if ~isfield(st.vols{i},'blobs')
        st.vols{i}.blobs=cell(1,1);
        bset = 1;
    else
        bset = numel(st.vols{i}.blobs)+1;
    end
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
        'max',mx, 'min',mn, 'colour',colour);
end


%==========================================================================
% function addtruecolourimage(handle,fname,colourmap,prop,mx,mn)
%==========================================================================
function addtruecolourimage(handle,fname,colourmap,prop,mx,mn)
% adds true colour image to current displayed image
global st
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    if ~isfield(st.vols{i},'blobs')
        st.vols{i}.blobs=cell(1,1);
        bset = 1;
    else
        bset = numel(st.vols{i}.blobs)+1;
    end
    c = struct('cmap', colourmap,'prop',prop);
    st.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
        'max',mx, 'min',mn, 'colour',c);
    addcolourbar(handle,bset);
end


%==========================================================================
% function addcolourbar(vh,bh)
%==========================================================================
function addcolourbar(vh,bh)
global st
if st.mode == 0,
    axpos = get(st.vols{vh}.ax{2}.ax,'Position');
else
    axpos = get(st.vols{vh}.ax{1}.ax,'Position');
end
H = ancestor(st.vols{vh}.ax{1}.ax,{'figure','uipanel'});
st.vols{vh}.blobs{bh}.cbar = axes('Parent',H,...
    'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
    'Box','on', 'YDir','normal', 'XTickLabel',[], 'XTick',[]);
if isfield(st.vols{vh}.blobs{bh},'name')
    ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
end


%==========================================================================
% function rmblobs(handle)
%==========================================================================
function rmblobs(handle)
global st
for i=valid_handles(handle)
    if isfield(st.vols{i},'blobs')
        for j=1:numel(st.vols{i}.blobs)
            if isfield(st.vols{i}.blobs{j},'cbar') && ishandle(st.vols{i}.blobs{j}.cbar),
                delete(st.vols{i}.blobs{j}.cbar);
            end
        end
        st.vols{i} = rmfield(st.vols{i},'blobs');
    end
end


%==========================================================================
% function register(hreg)
%==========================================================================
function register(hreg)
global st
%tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles;
if ~isempty(h)
    tmp = st.vols{h(1)}.ax{1}.ax;
    st.registry = struct('hReg',hreg,'hMe', tmp);
    spm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'spm_orthviews');
else
    warning('Nothing to register with');
end
st.centre = spm_XYZreg('GetCoords',st.registry.hReg);
st.centre = st.centre(:);


%==========================================================================
% function callback
%==========================================================================
function callback
global st
if ~iscell(st.callback), st.callback = { st.callback }; end
for i=1:numel(st.callback)
    if isa(st.callback{i},'function_handle')
        feval(st.callback{i});
    else
        eval(st.callback{i});
    end
end


%==========================================================================
% function xhairs(state)
%==========================================================================
function xhairs(state)
global st
if ~nargin, if st.xhairs, state = 'off'; else state = 'on'; end; end
st.xhairs = 0;
opt = 'on';
if ~strcmpi(state,'on')
    opt = 'off';
else
    st.xhairs = 1;
end
for i=valid_handles
    for j=1:3
        set(st.vols{i}.ax{j}.lx,'Visible',opt);
        set(st.vols{i}.ax{j}.ly,'Visible',opt);
    end
end


%==========================================================================
% function H = pos(handle)
%==========================================================================
function H = pos(handle)
global st
H = [];
for i=valid_handles(handle)
    is = inv(st.vols{i}.premul*st.vols{i}.mat);
    H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
end


%==========================================================================
% function my_reset
%==========================================================================
function my_reset
global st
if ~isempty(st) && isfield(st,'registry') && ishandle(st.registry.hMe)
    delete(st.registry.hMe); st = rmfield(st,'registry');
end
my_delete(1:max_img);
reset_st;


%==========================================================================
% function my_delete(handle)
%==========================================================================
function my_delete(handle)
global st
% remove blobs (and colourbars, if any)
rmblobs(handle);
% remove displayed axes
for i=valid_handles(handle)
    for j=1:3
        try
            set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
            delete(st.vols{i}.ax{j}.ax);
        end
    end
    st.vols{i} = [];
end


%==========================================================================
% function resolution(res)
%==========================================================================
function resolution(res)
global st
if ~nargin, res = 1; end % Default minimum resolution 1mm
for i=valid_handles
    % adapt resolution to smallest voxel size of displayed images
    res  = min([res,sqrt(sum((st.vols{i}.mat(1:3,1:3)).^2))]);
end
res      = res/mean(svd(st.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st.Space = st.Space*Mat;
st.bb    = st.bb/res;


%==========================================================================
% function move(handle,pos)
%==========================================================================
function move(handle,pos)
global st
for i=valid_handles(handle)
    st.vols{i}.area = pos;
end
bbox;
% redraw(valid_handles(handle));


%==========================================================================
% function bb = maxbb
%==========================================================================
function bb = maxbb
global st
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles
    premul = st.Space \ st.vols{i}.premul;
    bb = spm_get_bbox(st.vols{i}, 'fv', premul);
    mx = max([bb ; mx]);
    mn = min([bb ; mn]);
end
bb = [mn ; mx];


%==========================================================================
% function space(handle,M,dim)
%==========================================================================
function space(handle,M,dim)
global st
if isempty(st.vols{handle}), return; end
if nargin < 2
    M    = st.vols{handle}.mat;
    dim  = st.vols{handle}.dim(1:3);
end
Mat      = st.vols{handle}.premul(1:3,1:3)*M(1:3,1:3);
vox      = sqrt(sum(Mat.^2));
if det(Mat(1:3,1:3))<0, vox(1) = -vox(1); end
Mat      = diag([vox 1]);
Space    = (M)/Mat;
bb       = [1 1 1; dim];
bb       = [bb [1;1]];
bb       = bb*Mat';
bb       = bb(:,1:3);
bb       = sort(bb);
st.Space = Space;
st.bb    = bb;


%==========================================================================
% function zoom_op(fov,res)
%==========================================================================
function zoom_op(fov,res)
global st
if nargin < 1, fov = Inf; end
if nargin < 2, res = Inf; end

if isinf(fov)
    st.bb = maxbb;
elseif isnan(fov) || fov == 0
    current_handle = valid_handles;
    if numel(current_handle) > 1 % called from check reg context menu
        current_handle = get_current_handle;
    end
    if fov == 0
        % zoom to bounding box of current image ~= 0
        thr = 'nz';
    else
        % zoom to bounding box of current image > chosen threshold
        thr = spm_input('Threshold (Y > ...)', '+1', 'r', '0', 1);
    end
    premul = st.Space \ st.vols{current_handle}.premul;
    st.bb = spm_get_bbox(st.vols{current_handle}, thr, premul);
else
    vx    = sqrt(sum(st.Space(1:3,1:3).^2));
    vx    = vx.^(-1);
    pos   = spm_orthviews('pos');
    pos   = st.Space\[pos ; 1];
    pos   = pos(1:3)';
    st.bb = [pos-fov*vx; pos+fov*vx];
end
resolution(res);
bbox;
redraw_all;
if isfield(st.vols{1},'sdip')
    spm_eeg_inv_vbecd_disp('RedrawDip');
end


%==========================================================================
% function repos_start(varargin)
% function repos_move(varargin)
% function repos_end(varargin)
%==========================================================================
function repos_start(varargin)
% don't use right mouse button to start reposition
if ~strcmpi(get(gcbf,'SelectionType'),'alt')
    set(gcbf,'windowbuttonmotionfcn',@repos_move, 'windowbuttonupfcn',@repos_end);
    spm_orthviews('reposition');
end

function repos_move(varargin)
spm_orthviews('reposition');

function repos_end(varargin)
set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');

function scroll_wheel(varargin)
ScrollWheelData = varargin{2};
% detect if mouse is over an axis and change 
centre = spm_orthviews('Pos');
spm_orthviews('setcoords',centre+[0 0 ScrollWheelData.VerticalScrollCount]');
%spm_orthviews('reposition');


%==========================================================================
% function bbox
%==========================================================================
function bbox(handles)
global st
if ~nargin, handles = valid_handles; end

Dims = diff(st.bb)'+1;

TD = Dims([1 2])';
CD = Dims([1 3])';
if st.mode == 0, SD = Dims([3 2])'; else SD = Dims([2 3])'; end

for i=valid_handles(handles)
    
    H     = ancestor(st.vols{i}.ax{1}.ax,{'figure','uipanel'});
    un    = get(H,'Units');set(H,'Units','Pixels');
    sz    = get(H,'Position');set(H,'Units',un);
    sz    = sz(3:4);
    sz(2) = sz(2)-40;
    
    area   = st.vols{i}.area(:);
    area   = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
    if st.mode == 0
        sx = area(3)/(Dims(1)+Dims(3))/1.02;
    else
        sx = area(3)/(Dims(1)+Dims(2))/1.02;
    end
    sy     = area(4)/(Dims(2)+Dims(3))/1.02;
    s      = min([sx sy]);
    
    offy   = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
    sky    = s*(Dims(2)+Dims(3))*0.02;
    if st.mode == 0
        offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
        skx  = s*(Dims(1)+Dims(3))*0.02;
    else
        offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
        skx  = s*(Dims(1)+Dims(2))*0.02;
    end
    
    % Transverse
    set(st.vols{i}.ax{1}.ax,'Units','pixels', ...
        'Position',[offx offy s*Dims(1) s*Dims(2)],...
        'Units','normalized','Xlim',[0 TD(1)]+0.5,'Ylim',[0 TD(2)]+0.5,...
        'Visible','on','XTick',[],'YTick',[]);
    
    % Coronal
    set(st.vols{i}.ax{2}.ax,'Units','Pixels',...
        'Position',[offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)],...
        'Units','normalized','Xlim',[0 CD(1)]+0.5,'Ylim',[0 CD(2)]+0.5,...
        'Visible','on','XTick',[],'YTick',[]);
    
    % Sagittal
    if st.mode == 0
        set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    else
        set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    end
end


%==========================================================================
% function mx = maxval(vol)
%==========================================================================
function mx = maxval(vol)
if isstruct(vol)
    mx = -Inf;
    for i=1:vol.dim(3)
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imx = max(tmp(isfinite(tmp)));
        if ~isempty(imx), mx = max(mx,imx); end
    end
else
    mx = max(vol(isfinite(vol)));
end


%==========================================================================
% function mn = minval(vol)
%==========================================================================
function mn = minval(vol)
if isstruct(vol)
    mn = Inf;
    for i=1:vol.dim(3)
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imn = min(tmp(isfinite(tmp)));
        if ~isempty(imn), mn = min(mn,imn); end
    end
else
    mn = min(vol(isfinite(vol)));
end


%==========================================================================
% function redraw(arg1)
%==========================================================================
function redraw(arg1)
global st
bb   = st.bb;
Dims = round(diff(bb)'+1);
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);

for i = valid_handles(arg1)
    M = st.Space\st.vols{i}.premul*st.vols{i}.mat;
    TM0 = [ 1 0 0 -bb(1,1)+1
            0 1 0 -bb(1,2)+1
            0 0 1 -cent(3)
            0 0 0 1];
    TM = inv(TM0*M);
    TD = Dims([1 2]);
    
    CM0 = [ 1 0 0 -bb(1,1)+1
            0 0 1 -bb(1,3)+1
            0 1 0 -cent(2)
            0 0 0 1];
    CM = inv(CM0*M);
    CD = Dims([1 3]);
    
    if st.mode == 0
        SM0 = [ 0 0 1 -bb(1,3)+1
                0 1 0 -bb(1,2)+1
                1 0 0 -cent(1)
                0 0 0 1];
        SM = inv(SM0*M); 
        SD = Dims([3 2]);
    else
        SM0 = [ 0 -1 0 +bb(2,2)+1
                0  0 1 -bb(1,3)+1
                1  0 0 -cent(1)
                0  0 0 1];
        SM = inv(SM0*M);
        SD = Dims([2 3]);
    end
    
    try
        imgt = spm_slice_vol(st.vols{i},TM,TD,st.hld)';
        imgc = spm_slice_vol(st.vols{i},CM,CD,st.hld)';
        imgs = spm_slice_vol(st.vols{i},SM,SD,st.hld)';
    catch
        fprintf('Cannot access file "%s".\n', st.vols{i}.fname);
        fprintf('%s\n',getfield(lasterror,'message'));
        continue;
    end
    
    % get min/max threshold
    if strcmp(st.vols{i}.window,'auto')
        mn = -Inf;
        mx = Inf;
    else
        mn = min(st.vols{i}.window);
        mx = max(st.vols{i}.window);
        % threshold images
        imgt = max(imgt,mn); imgt = min(imgt,mx);
        imgc = max(imgc,mn); imgc = min(imgc,mx);
        imgs = max(imgs,mn); imgs = min(imgs,mx);
    end
    
    % compute intensity mapping, if histeq is available
    if license('test','image_toolbox') == 0
        st.vols{i}.mapping = 'linear';
    end
    switch st.vols{i}.mapping
        case 'linear'
        case 'histeq'
            % scale images to a range between 0 and 1
            imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
            imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
            imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
            img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
            imgt = reshape(img(1:numel(imgt1)),size(imgt1));
            imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
            imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
            mn = 0;
            mx = 1;
        case 'quadhisteq'
            % scale images to a range between 0 and 1
            imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
            imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
            imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
            img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
            imgt = reshape(img(1:numel(imgt1)),size(imgt1));
            imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
            imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
            mn = 0;
            mx = 1;
        case 'loghisteq'
            sw = warning('off','MATLAB:log:logOfZero');
            imgt = log(imgt-min(imgt(:)));
            imgc = log(imgc-min(imgc(:)));
            imgs = log(imgs-min(imgs(:)));
            warning(sw);
            imgt(~isfinite(imgt)) = 0;
            imgc(~isfinite(imgc)) = 0;
            imgs(~isfinite(imgs)) = 0;
            % scale log images to a range between 0 and 1
            imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
            imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
            imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
            img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
            imgt = reshape(img(1:numel(imgt1)),size(imgt1));
            imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
            imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
            mn = 0;
            mx = 1;
    end
    % recompute min/max for display
    if strcmp(st.vols{i}.window,'auto')
        mx = -Inf;
        mn = Inf;
    end
    tmp = imgt(isfinite(imgt));
    mx  = max([mx max(tmp)]);
    mn  = min([mn min(tmp)]);
    tmp = imgc(isfinite(imgc));
    mx  = max([mx max(tmp)]);
    mn  = min([mn min(tmp)]);
    tmp = imgs(isfinite(imgs));
    mx  = max([mx max(tmp)]);
    mn  = min([mn min(tmp)]);
    if mx==mn, mx = mn + eps; end

    if isfield(st.vols{i},'blobs')
        if ~isfield(st.vols{i}.blobs{1},'colour')
            % Add blobs for display using the split colourmap
            scal = 64/(mx-mn);
            dcoff = -mn*scal;
            imgt = imgt*scal+dcoff;
            imgc = imgc*scal+dcoff;
            imgs = imgs*scal+dcoff;

            if isfield(st.vols{i}.blobs{1},'max')
                mx = st.vols{i}.blobs{1}.max;
            else
                mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
                st.vols{i}.blobs{1}.max = mx;
            end
            if isfield(st.vols{i}.blobs{1},'min')
                mn = st.vols{i}.blobs{1}.min;
            else
                mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
                st.vols{i}.blobs{1}.min = mn;
            end

            vol  = st.vols{i}.blobs{1}.vol;
            M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
            tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
            tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
            tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';

            %tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
            %tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
            %tmps_z = find(tmps==0);tmps(tmps_z) = NaN;

            sc   = 64/(mx-mn);
            off  = 65.51-mn*sc;
            msk  = find(isfinite(tmpt)); imgt(msk) = off+tmpt(msk)*sc;
            msk  = find(isfinite(tmpc)); imgc(msk) = off+tmpc(msk)*sc;
            msk  = find(isfinite(tmps)); imgs(msk) = off+tmps(msk)*sc;

            cmap = get(st.fig,'Colormap');
            if size(cmap,1)~=128
                figure(st.fig)
                spm_colourmap('gray-hot')
            end
            redraw_colourbar(i,1,[mn mx],(1:64)'+64);
        elseif isstruct(st.vols{i}.blobs{1}.colour)
            % Add blobs for display using a defined colourmap

            % colourmaps
            gryc = (0:63)'*ones(1,3)/63;

            % scale grayscale image, not isfinite -> black
            gimgt = scaletocmap(imgt,mn,mx,gryc,65);
            gimgc = scaletocmap(imgc,mn,mx,gryc,65);
            gimgs = scaletocmap(imgs,mn,mx,gryc,65);
            gryc  = [gryc; 0 0 0];
            cactp = 0;

            for j=1:numel(st.vols{i}.blobs)
                % colourmaps
                actc = st.vols{i}.blobs{j}.colour.cmap;
                actp = st.vols{i}.blobs{j}.colour.prop;

                % get min/max for blob image
                if isfield(st.vols{i}.blobs{j},'max')
                    cmx = st.vols{i}.blobs{j}.max;
                else
                    cmx = max([eps maxval(st.vols{i}.blobs{j}.vol)]);
                end
                if isfield(st.vols{i}.blobs{j},'min')
                    cmn = st.vols{i}.blobs{j}.min;
                else
                    cmn = -cmx;
                end

                % get blob data
                vol  = st.vols{i}.blobs{j}.vol;
                M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
                tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';

                % actimg scaled round 0, black NaNs
                topc = size(actc,1)+1;
                tmpt = scaletocmap(tmpt,cmn,cmx,actc,topc);
                tmpc = scaletocmap(tmpc,cmn,cmx,actc,topc);
                tmps = scaletocmap(tmps,cmn,cmx,actc,topc);
                actc = [actc; 0 0 0];

                % combine gray and blob data to truecolour
                if isnan(actp)
                    if j==1, imgt = gryc(gimgt(:),:); end
                    imgt(tmpt~=size(actc,1),:) = actc(tmpt(tmpt~=size(actc,1)),:);
                    if j==1, imgc = gryc(gimgc(:),:); end
                    imgc(tmpc~=size(actc,1),:) = actc(tmpc(tmpc~=size(actc,1)),:);
                    if j==1, imgs = gryc(gimgs(:),:); end
                    imgs(tmps~=size(actc,1),:) = actc(tmps(tmps~=size(actc,1)),:);
                else
                    cactp = cactp + actp;
                    if j==1, imgt = actc(tmpt(:),:)*actp; else imgt = imgt + actc(tmpt(:),:)*actp; end
                    if j==numel(st.vols{i}.blobs), imgt = imgt + gryc(gimgt(:),:)*(1-cactp); end
                    if j==1, imgc = actc(tmpc(:),:)*actp; else imgc = imgc + actc(tmpc(:),:)*actp; end
                    if j==numel(st.vols{i}.blobs), imgc = imgc + gryc(gimgc(:),:)*(1-cactp); end
                    if j==1, imgs = actc(tmps(:),:)*actp; else imgs = imgs + actc(tmps(:),:)*actp; end
                    if j==numel(st.vols{i}.blobs), imgs = imgs + gryc(gimgs(:),:)*(1-cactp); end
                end
                if j==numel(st.vols{i}.blobs)
                    imgt = reshape(imgt,[size(gimgt) 3]);
                    imgc = reshape(imgc,[size(gimgc) 3]);
                    imgs = reshape(imgs,[size(gimgs) 3]);
                end

                 % colourbar
                csz   = size(st.vols{i}.blobs{j}.colour.cmap);
                cdata = reshape(st.vols{i}.blobs{j}.colour.cmap, [csz(1) 1 csz(2)]);
                redraw_colourbar(i,j,[cmn cmx],cdata);
            end

        else
            % Add full colour blobs - several sets at once
            scal  = 1/(mx-mn);
            dcoff = -mn*scal;

            wt = zeros(size(imgt));
            wc = zeros(size(imgc));
            ws = zeros(size(imgs));

            imgt  = repmat(imgt*scal+dcoff,[1,1,3]);
            imgc  = repmat(imgc*scal+dcoff,[1,1,3]);
            imgs  = repmat(imgs*scal+dcoff,[1,1,3]);

            cimgt = zeros(size(imgt));
            cimgc = zeros(size(imgc));
            cimgs = zeros(size(imgs));

            colour = zeros(numel(st.vols{i}.blobs),3);
            for j=1:numel(st.vols{i}.blobs) % get colours of all images first
                if isfield(st.vols{i}.blobs{j},'colour')
                    colour(j,:) = reshape(st.vols{i}.blobs{j}.colour, [1 3]);
                else
                    colour(j,:) = [1 0 0];
                end
            end
            %colour = colour/max(sum(colour));

            for j=1:numel(st.vols{i}.blobs)
                if isfield(st.vols{i}.blobs{j},'max')
                    mx = st.vols{i}.blobs{j}.max;
                else
                    mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
                    st.vols{i}.blobs{j}.max = mx;
                end
                if isfield(st.vols{i}.blobs{j},'min')
                    mn = st.vols{i}.blobs{j}.min;
                else
                    mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
                    st.vols{i}.blobs{j}.min = mn;
                end

                vol  = st.vols{i}.blobs{j}.vol;
                M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
                tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                % check min/max of sampled image
                % against mn/mx as given in st
                tmpt(tmpt(:)<mn) = mn;
                tmpc(tmpc(:)<mn) = mn;
                tmps(tmps(:)<mn) = mn;
                tmpt(tmpt(:)>mx) = mx;
                tmpc(tmpc(:)>mx) = mx;
                tmps(tmps(:)>mx) = mx;
                tmpt = (tmpt-mn)/(mx-mn);
                tmpc = (tmpc-mn)/(mx-mn);
                tmps = (tmps-mn)/(mx-mn);
                tmpt(~isfinite(tmpt)) = 0;
                tmpc(~isfinite(tmpc)) = 0;
                tmps(~isfinite(tmps)) = 0;

                cimgt = cimgt + cat(3,tmpt*colour(j,1),tmpt*colour(j,2),tmpt*colour(j,3));
                cimgc = cimgc + cat(3,tmpc*colour(j,1),tmpc*colour(j,2),tmpc*colour(j,3));
                cimgs = cimgs + cat(3,tmps*colour(j,1),tmps*colour(j,2),tmps*colour(j,3));

                wt = wt + tmpt;
                wc = wc + tmpc;
                ws = ws + tmps;
                cdata = permute(shiftdim((1/64:1/64:1)'* ...
                    colour(j,:),-1),[2 1 3]);
                redraw_colourbar(i,j,[mn mx],cdata);
            end

            imgt = repmat(1-wt,[1 1 3]).*imgt+cimgt;
            imgc = repmat(1-wc,[1 1 3]).*imgc+cimgc;
            imgs = repmat(1-ws,[1 1 3]).*imgs+cimgs;

            imgt(imgt<0)=0; imgt(imgt>1)=1;
            imgc(imgc<0)=0; imgc(imgc>1)=1;
            imgs(imgs<0)=0; imgs(imgs>1)=1;
        end
    else
        scal = 64/(mx-mn);
        dcoff = -mn*scal;
        imgt = imgt*scal+dcoff;
        imgc = imgc*scal+dcoff;
        imgs = imgs*scal+dcoff;
    end

    set(st.vols{i}.ax{1}.d,'HitTest','off', 'Cdata',imgt);
    set(st.vols{i}.ax{1}.lx,'HitTest','off',...
        'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
    set(st.vols{i}.ax{1}.ly,'HitTest','off',...
        'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));

    set(st.vols{i}.ax{2}.d,'HitTest','off', 'Cdata',imgc);
    set(st.vols{i}.ax{2}.lx,'HitTest','off',...
        'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
    set(st.vols{i}.ax{2}.ly,'HitTest','off',...
        'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));

    set(st.vols{i}.ax{3}.d,'HitTest','off','Cdata',imgs);
    if st.mode == 0
        set(st.vols{i}.ax{3}.lx,'HitTest','off',...
            'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
        set(st.vols{i}.ax{3}.ly,'HitTest','off',...
            'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
    else
        set(st.vols{i}.ax{3}.lx,'HitTest','off',...
            'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
        set(st.vols{i}.ax{3}.ly,'HitTest','off',...
            'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
    end

    if ~isempty(st.plugins) % process any addons
        for k = 1:numel(st.plugins)
            if isfield(st.vols{i},st.plugins{k})
                feval(['spm_ov_', st.plugins{k}], ...
                    'redraw', i, TM0, TD, CM0, CD, SM0, SD);
            end
        end
    end
end
drawnow;


%==========================================================================
% function redraw_all
%==========================================================================
function redraw_all
redraw(1:max_img);


%==========================================================================
% function redraw_colourbar(vh,bh,interval,cdata)
%==========================================================================
function redraw_colourbar(vh,bh,interval,cdata)
global st
if isfield(st.vols{vh}.blobs{bh},'cbar')
    if st.mode == 0
        axpos = get(st.vols{vh}.ax{2}.ax,'Position');
    else
        axpos = get(st.vols{vh}.ax{1}.ax,'Position');
    end
    % only scale cdata if we have out-of-range truecolour values
    if ndims(cdata)==3 && max(cdata(:))>1
        cdata=cdata./max(cdata(:));
    end
    image([0 1],interval,cdata,'Parent',st.vols{vh}.blobs{bh}.cbar);
    set(st.vols{vh}.blobs{bh}.cbar, ...
        'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1)...
        (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
        'YDir','normal','XTickLabel',[],'XTick',[]);
    if isfield(st.vols{vh}.blobs{bh},'name')
        ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
    end
end


%==========================================================================
% function centre = findcent
%==========================================================================
function centre = findcent
global st
obj    = get(st.fig,'CurrentObject');
centre = [];
cent   = [];
cp     = [];
for i=valid_handles
    for j=1:3
        if ~isempty(obj)
            if (st.vols{i}.ax{j}.ax == obj)
                cp = get(obj,'CurrentPoint');
            end
        end
        if ~isempty(cp)
            cp   = cp(1,1:2);
            is   = inv(st.Space);
            cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
            switch j
                case 1
                    cent([1 2])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,2)-1];
                case 2
                    cent([1 3])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,3)-1];
                case 3
                    if st.mode ==0
                        cent([3 2])=[cp(1)+st.bb(1,3)-1 cp(2)+st.bb(1,2)-1];
                    else
                        cent([2 3])=[st.bb(2,2)+1-cp(1) cp(2)+st.bb(1,3)-1];
                    end
            end
            break;
        end
    end
    if ~isempty(cent), break; end
end
if ~isempty(cent), centre = st.Space(1:3,1:3)*cent(:) + st.Space(1:3,4); end


%==========================================================================
% function handles = valid_handles(handles)
%==========================================================================
function handles = valid_handles(handles)
global st
if ~nargin, handles = 1:max_img; end
if isempty(st) || ~isfield(st,'vols')
    handles = [];
else
    handles = handles(:)';
    handles = handles(handles<=max_img & handles>=1 & ~rem(handles,1));
    for h=handles
        if isempty(st.vols{h}), handles(handles==h)=[]; end
    end
end


%==========================================================================
% function reset_st
%==========================================================================
function reset_st
global st
fig = spm_figure('FindWin','Graphics');
bb  = []; %[ [-78 78]' [-112 76]' [-50 85]' ];
st  = struct('n', 0, 'vols',{cell(max_img,1)}, 'bb',bb, 'Space',eye(4), ...
             'centre',[0 0 0], 'callback',';', 'xhairs',1, 'hld',1, ...
             'fig',fig, 'mode',1, 'plugins',{{}}, 'snap',[]);

xTB = spm('TBs');
if ~isempty(xTB)
    pluginbase = {spm('Dir') xTB.dir};
else
    pluginbase = {spm('Dir')};
end
for k = 1:numel(pluginbase)
    pluginpath = fullfile(pluginbase{k},'spm_orthviews');
    if exist(pluginpath,'dir') == 7
        pluginfiles = dir(fullfile(pluginpath,'spm_ov_*.m'));
        if ~isempty(pluginfiles)
            if ~isdeployed, addpath(pluginpath); end
            for l = 1:numel(pluginfiles)
                pluginname = spm_file(pluginfiles(l).name,'basename');
                st.plugins{end+1} = strrep(pluginname, 'spm_ov_','');
            end
        end
    end
end


%==========================================================================
% function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
%==========================================================================
function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1; end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(img<1)   = 1;
img(img>cml) = cml;
img(~isfinite(img)) = miscol;


%==========================================================================
% function cmap = getcmap(acmapname)
%==========================================================================
function cmap = getcmap(acmapname)
% get colormap of name acmapname
if ~isempty(acmapname)
    cmap = evalin('base',acmapname,'[]');
    if isempty(cmap) % not a matrix, is .mat file?
        acmat = spm_file(acmapname, 'ext','.mat');
        if exist(acmat, 'file')
            s    = struct2cell(load(acmat));
            cmap = s{1};
        end
    end
end
if size(cmap, 2)~=3
    warning('Colormap was not an N by 3 matrix')
    cmap = [];
end


%==========================================================================
% function item_parent = addcontext(volhandle)
%==========================================================================
function item_parent = addcontext(volhandle)
global st
% create context menu
set(0,'CurrentFigure',st.fig);
% contextmenu
item_parent = uicontextmenu('Callback',@(obj,evt) cm_pos);

% contextsubmenu 0
item00 = uimenu(item_parent, 'Label','unknown image', 'UserData','filename');
spm_orthviews('context_menu','image_info',item00,volhandle);
item0a = uimenu(item_parent, 'UserData','pos_mm', 'Separator','on', ...
    'Callback','spm_orthviews(''context_menu'',''repos_mm'');');
item0b = uimenu(item_parent, 'UserData','pos_vx', ...
    'Callback','spm_orthviews(''context_menu'',''repos_vx'');');
item0c = uimenu(item_parent, 'UserData','v_value');

% contextsubmenu 1
item1    = uimenu(item_parent,'Label','Zoom', 'Separator','on');
[zl, rl] = spm_orthviews('ZoomMenu');
for cz = numel(zl):-1:1
    if isinf(zl(cz))
        czlabel = 'Full Volume';
    elseif isnan(zl(cz))
        czlabel = 'BBox, this image > ...';
    elseif zl(cz) == 0
        czlabel = 'BBox, this image nonzero';
    else
        czlabel = sprintf('%dx%d mm', 2*zl(cz), 2*zl(cz));
    end
    item1_x = uimenu(item1, 'Label',czlabel,...
        'Callback', sprintf(...
        'spm_orthviews(''context_menu'',''zoom'',%d,%d)',zl(cz),rl(cz)));
    if isinf(zl(cz)) % default display is Full Volume
        set(item1_x, 'Checked','on');
    end
end

% contextsubmenu 2
checked   = {'off','off'};
checked{st.xhairs+1} = 'on';
item2     = uimenu(item_parent,'Label','Crosshairs','Callback','spm_orthviews(''context_menu'',''Xhair'');','Checked',checked{2});

% contextsubmenu 3

if any(any(abs(st.Space - eye(4)) > 1e-4))
    checked = {'on', 'off'};
else
    checked = {'off', 'on'};
end
item3     = uimenu(item_parent,'Label','Orientation');
item3_1   = uimenu(item3,      'Label','World space', 'Callback','spm_orthviews(''context_menu'',''orientation'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Voxel space (1st image)', 'Callback','spm_orthviews(''context_menu'',''orientation'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Voxel space (this image)', 'Callback','spm_orthviews(''context_menu'',''orientation'',1);','Checked','off');

% contextsubmenu 3
if isempty(st.snap)
    checked = {'off', 'on'};
else
    checked = {'on', 'off'};
end
item3     = uimenu(item_parent,'Label','Snap to Grid');
item3_1   = uimenu(item3,      'Label','Don''t snap', 'Callback','spm_orthviews(''context_menu'',''snap'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Snap to 1st image', 'Callback','spm_orthviews(''context_menu'',''snap'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Snap to this image', 'Callback','spm_orthviews(''context_menu'',''snap'',1);','Checked','off');

% contextsubmenu 4
if st.hld == 0
    checked = {'off', 'off', 'on'};
elseif st.hld > 0
    checked = {'off', 'on', 'off'};
else
    checked = {'on', 'off', 'off'};
end
item4     = uimenu(item_parent,'Label','Interpolation');
item4_1   = uimenu(item4,      'Label','NN',    'Callback','spm_orthviews(''context_menu'',''interpolation'',3);', 'Checked',checked{3});
item4_2   = uimenu(item4,      'Label','Trilin', 'Callback','spm_orthviews(''context_menu'',''interpolation'',2);','Checked',checked{2});
item4_3   = uimenu(item4,      'Label','Sinc',  'Callback','spm_orthviews(''context_menu'',''interpolation'',1);','Checked',checked{1});

% contextsubmenu 5
% item5     = uimenu(item_parent,'Label','Position', 'Callback','spm_orthviews(''context_menu'',''position'');');

% contextsubmenu 6
item6       = uimenu(item_parent,'Label','Intensity','Separator','on');
item6_1     = uimenu(item6,      'Label','Range');
item6_1_1   = uimenu(item6_1,    'Label','This image');
item6_1_1_1 = uimenu(item6_1_1,  'Label','Automatic', 'Callback','spm_orthviews(''context_menu'',''window'',2);');
item6_1_1_2 = uimenu(item6_1_1,  'Label','Manual', 'Callback','spm_orthviews(''context_menu'',''window'',1);');
item6_1_1_3 = uimenu(item6_1_1,  'Label','Percentiles', 'Callback','spm_orthviews(''context_menu'',''window'',3);');
item6_1_2   = uimenu(item6_1,    'Label','All images');
item6_1_2_1 = uimenu(item6_1_2,  'Label','Automatic', 'Callback','spm_orthviews(''context_menu'',''window_gl'',2);');
item6_1_2_2 = uimenu(item6_1_2,  'Label','Manual', 'Callback','spm_orthviews(''context_menu'',''window_gl'',1);');
if license('test','image_toolbox') == 1
    offon = {'off', 'on'};
    checked = offon(strcmp(st.vols{volhandle}.mapping, ...
        {'linear', 'histeq', 'loghisteq', 'quadhisteq'})+1);
    item6_2     = uimenu(item6,      'Label','Mapping');
    item6_2_1   = uimenu(item6_2,    'Label','This image');
    item6_2_1_1 = uimenu(item6_2_1,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping'',''linear'');');
    item6_2_1_2 = uimenu(item6_2_1,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping'',''histeq'');');
    item6_2_1_3 = uimenu(item6_2_1,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping'',''loghisteq'');');
    item6_2_1_4 = uimenu(item6_2_1,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping'',''quadhisteq'');');
    item6_2_2   = uimenu(item6_2,    'Label','All images');
    item6_2_2_1 = uimenu(item6_2_2,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''linear'');');
    item6_2_2_2 = uimenu(item6_2_2,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''histeq'');');
    item6_2_2_3 = uimenu(item6_2_2,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''loghisteq'');');
    item6_2_2_4 = uimenu(item6_2_2,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''quadhisteq'');');
end

% contextsubmenu 7
item7     = uimenu(item_parent,'Label','Overlay');
item7_1   = uimenu(item7,      'Label','Add blobs');
item7_1_1 = uimenu(item7_1,    'Label','This image',  'Callback','spm_orthviews(''context_menu'',''add_blobs'',2);');
item7_1_2 = uimenu(item7_1,    'Label','All images', 'Callback','spm_orthviews(''context_menu'',''add_blobs'',1);');
item7_2   = uimenu(item7,      'Label','Add image');
item7_2_1 = uimenu(item7_2,    'Label','This image',  'Callback','spm_orthviews(''context_menu'',''add_image'',2);');
item7_2_2 = uimenu(item7_2,    'Label','All images', 'Callback','spm_orthviews(''context_menu'',''add_image'',1);');
item7_3   = uimenu(item7,      'Label','Add coloured blobs','Separator','on');
item7_3_1 = uimenu(item7_3,    'Label','This image',  'Callback','spm_orthviews(''context_menu'',''add_c_blobs'',2);');
item7_3_2 = uimenu(item7_3,    'Label','All images', 'Callback','spm_orthviews(''context_menu'',''add_c_blobs'',1);');
item7_4   = uimenu(item7,      'Label','Add coloured image');
item7_4_1 = uimenu(item7_4,    'Label','This image',  'Callback','spm_orthviews(''context_menu'',''add_c_image'',2);');
item7_4_2 = uimenu(item7_4,    'Label','All images', 'Callback','spm_orthviews(''context_menu'',''add_c_image'',1);');
item7_5   = uimenu(item7,      'Label','Remove blobs',        'Visible','off','Separator','on');
item7_6   = uimenu(item7,      'Label','Remove coloured blobs','Visible','off');
item7_6_1 = uimenu(item7_6,    'Label','This image', 'Visible','on');
item7_6_2 = uimenu(item7_6,    'Label','All images','Visible','on');
item7_7   = uimenu(item7,      'Label','Set blobs max', 'Visible','off');

for i=1:3
    set(st.vols{volhandle}.ax{i}.ax,'UIcontextmenu',item_parent);
    st.vols{volhandle}.ax{i}.cm = item_parent;
end

% process any plugins
for k = 1:numel(st.plugins)
    feval(['spm_ov_', st.plugins{k}],'context_menu',volhandle,item_parent);
    if k==1
        h = get(item_parent,'Children');
        set(h(1),'Separator','on'); 
    end
end


%==========================================================================
% function addcontexts(handles)
%==========================================================================
function addcontexts(handles)
for ii = valid_handles(handles)
    addcontext(ii);
end


%==========================================================================
% function rmcontexts(handles)
%==========================================================================
function rmcontexts(handles)
global st
for ii = valid_handles(handles)
    for i=1:3
        set(st.vols{ii}.ax{i}.ax,'UIcontextmenu',[]);
        try, st.vols{ii}.ax{i} = rmfield(st.vols{ii}.ax{i},'cm'); end
    end
end


%==========================================================================
% function c_menu(varargin)
%==========================================================================
function c_menu(varargin)
global st

switch lower(varargin{1})
    case 'image_info'
        if nargin <3
            current_handle = get_current_handle;
        else
            current_handle = varargin{3};
        end
        if isfield(st.vols{current_handle},'fname')
            [p,n,e,v] = spm_fileparts(st.vols{current_handle}.fname);
            if isfield(st.vols{current_handle},'n')
                v = sprintf(',%d',st.vols{current_handle}.n);
            end
            set(varargin{2}, 'Label',[n e v]);
        end
        delete(get(varargin{2},'children'));
        if exist('p','var')
            item1 = uimenu(varargin{2}, 'Label', p);
        end
        if isfield(st.vols{current_handle},'descrip')
            item2 = uimenu(varargin{2}, 'Label',...
                st.vols{current_handle}.descrip);
        end
        dt = st.vols{current_handle}.dt(1);
        item3 = uimenu(varargin{2}, 'Label', sprintf('Data type: %s', spm_type(dt)));
        str   = 'Intensity: varied';
        if size(st.vols{current_handle}.pinfo,2) == 1
            if st.vols{current_handle}.pinfo(2)
                str = sprintf('Intensity: Y = %g X + %g',...
                    st.vols{current_handle}.pinfo(1:2)');
            else
                str = sprintf('Intensity: Y = %g X', st.vols{current_handle}.pinfo(1)');
            end
        end
        item4  = uimenu(varargin{2}, 'Label',str);
        item5  = uimenu(varargin{2}, 'Label', 'Image dimensions', 'Separator','on');
        item51 = uimenu(varargin{2}, 'Label',...
            sprintf('%dx%dx%d', st.vols{current_handle}.dim(1:3)));
        
        prms   = spm_imatrix(st.vols{current_handle}.mat);
        item6  = uimenu(varargin{2}, 'Label', 'Voxel size', 'Separator','on');
        item61 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', prms(7:9)));
        
        O      = st.vols{current_handle}.mat\[0 0 0 1]'; O=O(1:3)';
        item7  = uimenu(varargin{2}, 'Label', 'Origin', 'Separator','on');
        item71 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', O));
        
        R      = spm_matrix([0 0 0 prms(4:6)]);
        item8  = uimenu(varargin{2}, 'Label', 'Rotations', 'Separator','on');
        item81 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(1,1:3)));
        item82 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(2,1:3)));
        item83 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(3,1:3)));
        item9  = uimenu(varargin{2},...
            'Label','Specify other image...',...
            'Callback','spm_orthviews(''context_menu'',''swap_img'');',...
            'Separator','on');
        
    case 'repos_mm'
        oldpos_mm = spm_orthviews('pos');
        newpos_mm = spm_input('New Position (mm)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_mm),3);
        spm_orthviews('reposition',newpos_mm);
        
    case 'repos_vx'
        current_handle = get_current_handle;
        oldpos_vx = spm_orthviews('pos', current_handle);
        newpos_vx = spm_input('New Position (voxels)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_vx),3);
        newpos_mm = st.vols{current_handle}.mat*[newpos_vx;1];
        spm_orthviews('reposition',newpos_mm(1:3));
        
    case 'zoom'
        zoom_all(varargin{2:end});
        bbox;
        redraw_all;
        
    case 'xhair'
        spm_orthviews('Xhairs',varargin{2:end});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = findobj(cm_handles(i),'label','Crosshairs');
            if st.xhairs
                set(z_handle,'Checked','on');
            else
                set(z_handle,'Checked','off');
            end
        end
        
    case 'orientation'
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Orientation'),'Children');
            set(z_handle,'Checked','off');
        end
        if varargin{2} == 3
            spm_orthviews('Space');
            for i = 1:numel(cm_handles),
                z_handle = findobj(cm_handles(i),'label','World space');
                set(z_handle,'Checked','on');
            end
        elseif varargin{2} == 2,
            spm_orthviews('Space',1);
            for i = 1:numel(cm_handles)
                z_handle = findobj(cm_handles(i),'label',...
                    'Voxel space (1st image)');
                set(z_handle,'Checked','on');
            end
        else
            spm_orthviews('Space',get_current_handle);
            z_handle = findobj(st.vols{get_current_handle}.ax{1}.cm, ...
                'label','Voxel space (this image)');
            set(z_handle,'Checked','on');
            return;
        end
        
    case 'snap'
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle,'Checked','off');
        end
        if varargin{2} == 3
            st.snap = [];
        elseif varargin{2} == 2
            st.snap = 1;
        else
            st.snap = get_current_handle;
            z_handle = get(findobj(st.vols{get_current_handle}.ax{1}.cm,'label','Snap to Grid'),'Children');
            set(z_handle(1),'Checked','on');
            return;
        end
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle(varargin{2}),'Checked','on');
        end
        
    case 'interpolation'
        tmp        = [-4 1 0];
        st.hld     = tmp(varargin{2});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Interpolation'),'Children');
            set(z_handle,'Checked','off');
            set(z_handle(varargin{2}),'Checked','on');
        end
        redraw_all;
        
    case 'window'
        current_handle = get_current_handle;
        if varargin{2} == 2
            spm_orthviews('window',current_handle);
        elseif varargin{2} == 3
            pc = spm_input('Percentiles', '+1', 'w', '3 97', 2, 100);
            wn = spm_summarise(st.vols{current_handle}, 'all', ...
                @(X) spm_percentile(X, pc));
            spm_orthviews('window',current_handle,wn);
        else
            if isnumeric(st.vols{current_handle}.window)
                defstr = sprintf('%.2f %.2f', st.vols{current_handle}.window);
            else
                defstr = '';
            end
            [w,yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w,yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            spm_orthviews('window',current_handle,w);
        end
        
    case 'window_gl'
        if varargin{2} == 2
            for i = 1:numel(get_cm_handles)
                st.vols{i}.window = 'auto';
            end
        else
            current_handle = get_current_handle;
            if isnumeric(st.vols{current_handle}.window)
                defstr = sprintf('%d %d', st.vols{current_handle}.window);
            else
                defstr = '';
            end
            [w,yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w,yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            for i = 1:numel(get_cm_handles)
                st.vols{i}.window = w;
            end
        end
        redraw_all;
        
    case 'mapping'
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        current_handle = get_current_handle;
        cm_handles = get_cm_handles;
        st.vols{current_handle}.mapping = varargin{2};
        z_handle = get(findobj(cm_handles(current_handle), ...
            'label','Intensity mapping'),'Children');
        for k = 1:numel(z_handle)
            c_handle = get(z_handle(k), 'Children');
            set(c_handle, 'checked', 'off');
            set(c_handle(checked), 'checked', 'on');
        end
        redraw_all;
        
    case 'mapping_gl'
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        cm_handles = get_cm_handles;
        for k = valid_handles
            st.vols{k}.mapping = varargin{2};
            z_handle = get(findobj(cm_handles(k), ...
                'label','Intensity mapping'),'Children');
            for l = 1:numel(z_handle)
                c_handle = get(z_handle(l), 'Children');
                set(c_handle, 'checked', 'off');
                set(c_handle(checked), 'checked', 'on');
            end
        end
        redraw_all;
        
    case 'swap_img'
        current_handle = get_current_handle;
        newimg = spm_select(1,'image','select new image');
        if ~isempty(newimg)
            new_info = spm_vol(newimg);
            fn = fieldnames(new_info);
            for k=1:numel(fn)
                st.vols{current_handle}.(fn{k}) = new_info.(fn{k});
            end
            spm_orthviews('context_menu','image_info',get(gcbo, 'parent'));
            redraw_all;
        end
        
    case 'add_blobs'
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            for i = 1:numel(cm_handles)
                addblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M);
                % Add options for removing blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','This image','Callback','spm_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1,
                    item7_3_2 = uimenu(c_handle,'Label','All images','Callback','spm_orthviews(''context_menu'',''remove_blobs'',1);');
                end
                % Add options for setting maxima for blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Set blobs max');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                uimenu(c_handle,'Label','This image','Callback','spm_orthviews(''context_menu'',''setblobsmax'',2);');
                if varargin{2} == 1
                    uimenu(c_handle,'Label','All images','Callback','spm_orthviews(''context_menu'',''setblobsmax'',1);');
                end
            end
            redraw_all;
        end
        
    case 'remove_blobs'
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        for i = 1:numel(cm_handles)
            rmblobs(cm_handles(i));
            % Remove options for removing blobs
            c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove blobs');
            delete(get(c_handle,'Children'));
            set(c_handle,'Visible','off');
            % Remove options for setting maxima for blobs
            c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Set blobs max');
            set(c_handle,'Visible','off');
        end
        redraw_all;
        
    case 'add_image'
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        fname = spm_select(1,'image','select image');
        if ~isempty(fname)
            for i = 1:numel(cm_handles)
                addimage(cm_handles(i),fname);
                % Add options for removing blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','This image','Callback','spm_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1
                    item7_3_2 = uimenu(c_handle,'Label','All images','Callback','spm_orthviews(''context_menu'',''remove_blobs'',1);');
                end
                % Add options for setting maxima for blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Set blobs max');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                uimenu(c_handle,'Label','This image','Callback','spm_orthviews(''context_menu'',''setblobsmax'',2);');
                if varargin{2} == 1
                    uimenu(c_handle,'Label','All images','Callback','spm_orthviews(''context_menu'',''setblobsmax'',1);');
                end
            end
            redraw_all;
        end
        
    case 'add_c_blobs'
        % Add blobs to the image - in full colour
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            c = spm_input('Colour','+1','m',...
                'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',xSPM.title,c_names{c});
            for i = 1:numel(cm_handles)
                addcolouredblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M,colours(c,:),xSPM.title);
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove coloured blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1   = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',2,c);',...
                    'UserData',c);
                if varargin{2} == 1
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end
            end
            redraw_all;
        end
        
    case 'remove_c_blobs'
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
        for i = 1:numel(cm_handles)
            if isfield(st.vols{cm_handles(i)},'blobs')
                for j = 1:numel(st.vols{cm_handles(i)}.blobs)
                    if all(st.vols{cm_handles(i)}.blobs{j}.colour == colours(varargin{3},:));
                        if isfield(st.vols{cm_handles(i)}.blobs{j},'cbar')
                            delete(st.vols{cm_handles(i)}.blobs{j}.cbar);
                        end
                        st.vols{cm_handles(i)}.blobs(j) = [];
                        break;
                    end
                end
                rm_c_menu = findobj(st.vols{cm_handles(i)}.ax{1}.cm,'Label','Remove coloured blobs');
                delete(gcbo);
                if isempty(st.vols{cm_handles(i)}.blobs)
                    st.vols{cm_handles(i)} = rmfield(st.vols{cm_handles(i)},'blobs');
                    set(rm_c_menu, 'Visible', 'off');
                end
            end
        end
        redraw_all;
        
    case 'add_c_image'
        % Add truecolored image
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        fname = spm_select([1 Inf],'image','select image(s)');
        for k = 1:size(fname,1)
            c = spm_input(sprintf('Image %d: Colour',k),'+1','m',...
                'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',fname(k,:),c_names{c});
            for i = 1:numel(cm_handles)
                addcolouredimage(cm_handles(i),fname(k,:),colours(c,:));
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove coloured blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1 = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',2,c);','UserData',c);
                if varargin{2} == 1
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end
            end
            redraw_all;
        end
        
    case 'setblobsmax'
        if varargin{2} == 1
            % global
            cm_handles = valid_handles;
            mx = -inf;
            for i = 1:numel(cm_handles)
                if ~isfield(st.vols{cm_handles(i)}, 'blobs'), continue, end
                for j = 1:numel(st.vols{cm_handles(i)}.blobs)
                    mx = max(mx, st.vols{cm_handles(i)}.blobs{j}.max);
                end
            end
            mx = spm_input('Maximum value', '+1', 'r', mx, 1);
            for i = 1:numel(cm_handles)
                if ~isfield(st.vols{cm_handles(i)}, 'blobs'), continue, end
                for j = 1:numel(st.vols{cm_handles(i)}.blobs)
                    st.vols{cm_handles(i)}.blobs{j}.max = mx;
                end
            end
        else
            % local (should handle coloured blobs, but not implemented yet)
            cm_handle = get_current_handle;
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            if ~isfield(st.vols{cm_handle}, 'blobs'), return, end
            for j = 1:numel(st.vols{cm_handle}.blobs)
                if nargin < 4 || ...
                        all(st.vols{cm_handle}.blobs{j}.colour == colours(varargin{3},:))
                    mx = st.vols{cm_handle}.blobs{j}.max;
                    mx = spm_input('Maximum value', '+1', 'r', mx, 1);
                    st.vols{cm_handle}.blobs{j}.max = mx;
                end
            end
        end
        redraw_all;
end


%==========================================================================
% function current_handle = get_current_handle
%==========================================================================
function current_handle = get_current_handle
cm_handle      = get(gca,'UIContextMenu');
cm_handles     = get_cm_handles;
current_handle = find(cm_handles==cm_handle);


%==========================================================================
% function cm_pos
%==========================================================================
function cm_pos(handles)
global st
if ~nargin, handles = get_current_handle; end
for i=valid_handles(handles)
    if isfield(st.vols{i}.ax{1},'cm')
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_mm'),...
            'Label',sprintf('mm:  %.1f %.1f %.1f',spm_orthviews('pos')));
        pos = spm_orthviews('pos',i);
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_vx'),...
            'Label',sprintf('vx:  %.1f %.1f %.1f',pos));
        try
            Y = spm_sample_vol(st.vols{i},pos(1),pos(2),pos(3),st.hld);
        catch
            Y = NaN;
            fprintf('Cannot access file "%s".\n', st.vols{i}.fname);
        end
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','v_value'),...
            'Label',sprintf('Y = %g',Y));
    end
end


%==========================================================================
% function cm_handles = get_cm_handles
%==========================================================================
function cm_handles = get_cm_handles
global st
cm_handles = [];
for i = valid_handles
    try, cm_handles = [cm_handles st.vols{i}.ax{1}.cm]; end
end


%==========================================================================
% function zoom_all(zoom,res)
%==========================================================================
function zoom_all(zoom,res)
cm_handles = get_cm_handles;
zoom_op(zoom,res);
for i = 1:numel(cm_handles)
    z_handle = get(findobj(cm_handles(i),'label','Zoom'),'Children');
    set(z_handle,'Checked','off');
    if isinf(zoom)
        set(findobj(z_handle,'Label','Full Volume'),'Checked','on');
    elseif zoom > 0
        set(findobj(z_handle,'Label',sprintf('%dx%d mm', 2*zoom, 2*zoom)),'Checked','on');
    end % leave all unchecked if either bounding box option was chosen
end


%==========================================================================
% function m = max_img
%==========================================================================
function m = max_img
m = 24;
