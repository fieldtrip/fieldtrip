function [t,sts] = cfg_getfile(varargin)
% File selector
% FORMAT [t,sts] = cfg_getfile(n,typ,mesg,sel,wd,filt,frames)
%     n    - Number of files
%            A single value or a range.  e.g.
%            1       - Select one file
%            Inf     - Select any number of files
%            [1 Inf] - Select 1 to Inf files
%            [0 1]   - select 0 or 1 files
%            [10 12] - select from 10 to 12 files
%     typ  - file type
%           'any'   - all files
%           'batch' - SPM batch files (.m, .mat and XML)
%           'dir'   - select a directory
%           'image' - Image files (".img" and ".nii")
%                     Note that it gives the option to select
%                     individual volumes of the images.
%           'mat'   - Matlab .mat files or .txt files (assumed to contain
%                     ASCII representation of a 2D-numeric array)
%           'mesh'  - Mesh files (".gii" and ".mat")
%           'nifti' - NIfTI files without the option to select frames
%           'xml'   - XML files
%           Other strings act as a filter to regexp.  This means
%           that e.g. DCM*.mat files should have a typ of '^DCM.*\.mat$'
%      mesg - a prompt (default 'Select files...')
%      sel  - list of already selected files
%      wd   - Directory to start off in
%      filt - value for user-editable filter (default '.*')
%      frames - Image frame numbers to include (default '1')
%
%      t    - selected files
%      sts  - status (1 means OK, 0 means window quit)
%
% FORMAT [t,ind] = cfg_getfile('Filter',files,typ,filt,frames)
% filter the list of files (cell array) in the same way as the
% GUI would do. There is an additional typ 'extimage' which will match
% images with frame specifications, too. The 'frames' argument
% is currently ignored, i.e. image files will not be filtered out if
% their frame numbers do not match.
% When filtering directory names, the filt argument will be applied to the
% last directory in a path only.
% t returns the filtered list (cell array), ind an index array, such that 
% t = files(ind).
%
% FORMAT cpath = cfg_getfile('CPath',path,cwd)
% function to canonicalise paths: Prepends cwd to relative paths, processes
% '..' & '.' directories embedded in path.
% path     - string matrix containing path name
% cwd      - current working directory [default '.']
% cpath    - conditioned paths, in same format as input path argument
%
% FORMAT [files,dirs]=cfg_getfile('List',direc,filt)
% Returns files matching the filter (filt) and directories within dire
% direc    - directory to search
% filt     - filter to select files with (see regexp) e.g. '^w.*\.img$'
% files    - files matching 'filt' in directory 'direc'
% dirs     - subdirectories of 'direc'
% FORMAT [files,dirs]=cfg_getfile('ExtList',direc,filt,frames)
% As above, but for selecting frames of 4D NIfTI files
% frames   - vector of frames to select (defaults to 1, if not
%            specified). If the frame number is Inf, all frames for the
%            matching images are listed. 
% FORMAT [files,dirs]=cfg_getfile('FPList',direc,filt)
% FORMAT [files,dirs]=cfg_getfile('ExtFPList',direc,filt,frames)
% As above, but returns files with full paths (i.e. prefixes direc to each)
% FORMAT [files,dirs]=cfg_getfile('FPListRec',direc,filt)
% FORMAT [files,dirs]=cfg_getfile('ExtFPListRec',direc,filt,frames)
% As above, but returns files with full paths (i.e. prefixes direc to
% each) and searches through sub directories recursively.
%
% FORMAT cfg_getfile('prevdirs',dir)
% Add directory dir to list of previous directories.
% FORMAT dirs=cfg_getfile('prevdirs')
% Retrieve list of previous directories.
%
% This code is based on the file selection dialog in SPM5, with virtual
% file handling turned off.
%____________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% John Ashburner and Volkmar Glauche
% $Id: cfg_getfile.m 6251 2014-10-30 13:58:44Z volkmar $

t = {};
sts = false;
if nargin > 0 && ischar(varargin{1})
    switch lower(varargin{1})
        case {'addvfiles', 'clearvfiles', 'vfiles'}
            cfg_message('matlabbatch:deprecated:vfiles', ...
                        'Trying to use deprecated ''%s'' call.', ...
                        lower(varargin{1}));
        case 'cpath'
            cfg_message(nargchk(2,Inf,nargin,'struct'));
            t = cpath(varargin{2:end});
            sts = true;
        case 'filter'
            filt    = mk_filter(varargin{3:end});
            t       = varargin{2};
            if numel(t) == 1 && isempty(t{1})
                sts = 1;
                return;
            end;
            t1 = cell(size(t));
            if any(strcmpi(varargin{3},{'dir','extdir'}))
                % only filter last directory in path
                for k = 1:numel(t)
                    t{k} = cpath(t{k});
                    if t{k}(end) == filesep
                        [p n] = fileparts(t{k}(1:end-1));
                    else
                        [p n] = fileparts(t{k});
                    end
                    if strcmpi(varargin{3},'extdir')
                        t1{k} = [n filesep];
                    else
                        t1{k} = n;
                    end
                end
            else
                % only filter filenames, not paths
                for k = 1:numel(t)
                    [p n e] = fileparts(t{k});
                    t1{k} = [n e];
                end
            end
            [t1,sts1] = do_filter(t1,filt.ext);
            [t1,sts2] = do_filter(t1,filt.filt);
            sts = sts1(sts2);
            t = t(sts);
        case {'list', 'fplist', 'extlist', 'extfplist'}
            if nargin > 3
                frames = varargin{4};
            else
                frames = 1; % (ignored in listfiles if typ==any)
            end;
            if regexpi(varargin{1}, 'ext') % use frames descriptor
                typ = 'extimage';
            else
                typ = 'any';
            end
            filt    = mk_filter(typ, varargin{3}, frames);
            [t sts] = listfiles(varargin{2}, filt); % (sts is subdirs here)
            sts = sts(~(strcmp(sts,'.')|strcmp(sts,'..'))); % remove '.' and '..' entries
            if regexpi(varargin{1}, 'fplist') % return full pathnames
                direc = cfg_getfile('cpath', varargin{2});
                % remove trailing path separator if present
                direc = regexprep(direc, [filesep '$'], '');
                if ~isempty(t)
                    t = strcat(direc, filesep, t);
                end
                if nargout > 1
                    % subdirs too
                    sts = cellfun(@(sts1)cpath(sts1, direc), sts, 'UniformOutput',false);
                end
            end
        case {'fplistrec', 'extfplistrec'}
            % list directory
            [f1 d1] = cfg_getfile(varargin{1}(1:end-3),varargin{2:end});
            f2 = cell(size(d1));
            d2 = cell(size(d1));
            for k = 1:numel(d1)
                % recurse into sub directories
                [f2{k} d2{k}] = cfg_getfile(varargin{1}, d1{k}, ...
                                            varargin{3:end});
            end
            t = vertcat(f1, f2{:});
            if nargout > 1
                sts = vertcat(d1, d2{:});
            end
        case 'prevdirs',
            if nargin > 1
                prevdirs(varargin{2});
            end;
            if nargout > 0 || nargin == 1
                t = prevdirs;
                sts = true;
            end;
        otherwise
            cfg_message('matlabbatch:usage','Inappropriate usage.');
    end
else
    [t,sts] = selector(varargin{:});
end
%=======================================================================

%=======================================================================
function [t,ok] = selector(n,typ,mesg,already,wd,filt,frames,varargin)
if nargin<1 || ~isnumeric(n) || numel(n) > 2 
    n = [0 Inf];
else
    if numel(n)==1,   n    = [n n];    end;
    if n(1)>n(2),     n    = n([2 1]); end;
    if ~isfinite(n(1)), n(1) = 0;        end;
end

if nargin<2 || ~ischar(typ), typ     = 'any';   end;

if nargin<3 || ~(ischar(mesg) || iscellstr(mesg))
    mesg    = 'Select files...';
elseif iscellstr(mesg)
    mesg = char(mesg);
end

if nargin<4 || isempty(already) || (iscell(already) && isempty(already{1}))
    already = {};
else
    % Add folders of already selected files to prevdirs list
    pd1 = cellfun(@(a1)strcat(fileparts(a1),filesep), already, ...
                  'UniformOutput',false);
    prevdirs(pd1);
end

if nargin<5 || isempty(wd) || ~ischar(wd)
    if isempty(already)
        wd      = pwd;
    else
        wd = fileparts(already{1});
        if isempty(wd)
            wd = pwd;
        end
    end;
end

if nargin<6 || ~ischar(filt), filt    = '.*';    end;

if nargin<7 || ~(isnumeric(frames) || ischar(frames))
    frames  = '1';
elseif isnumeric(frames)
    frames = char(gencode_rvalue(frames(:)'));
elseif ischar(frames)
    try
        ev = eval(frames);
        if ~isnumeric(ev)
            frames = '1';
        end
    catch
        frames = '1';
    end
end
ok  = 0;

t = '';
sfilt = mk_filter(typ,filt,eval(frames));

[col1,col2,col3,lf,bf] = colours;

% delete old selector, if any
fg = findobj(0,'Tag',mfilename);
if ~isempty(fg)
    delete(fg);
end
% create figure
fg = figure('IntegerHandle','off',...
    'Tag',mfilename,...
    'Name',mesg,...
    'NumberTitle','off',...
    'Units','Pixels',...
    'MenuBar','none',...
    'DefaultTextInterpreter','none',...
    'DefaultUicontrolInterruptible','on',...
    'Visible','off');
    
cfg_onscreen(fg);
set(fg,'Visible','on');

sellines = min([max([n(2) numel(already)]), 4]);
[pselp pcntp pfdp pdirp] = panelpositions(fg, sellines+1);

uicontrol(fg,...
          'style','text',...
          'units','normalized',...
          'position',posinpanel([0 sellines/(sellines+1) 1 1/(sellines+1)],pselp),...
          lf,...
          'BackgroundColor',get(fg,'Color'),...
          'ForegroundColor',col3,...
          'HorizontalAlignment','left',...
          'string',mesg,...
          'tag','msg');

% Selected Files
sel = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',posinpanel([0 0 1 sellines/(sellines+1)],pselp),...
    lf,...
    'Callback',@unselect,...
    'tag','selected',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'Max',10000,...
    'Min',0,...
    'String',already,...
    'Value',1);
c0 = uicontextmenu('Parent',fg);
set(sel,'uicontextmenu',c0);
uimenu('Label','Unselect All', 'Parent',c0,'Callback',@unselect_all);

% get cwidth for buttons
tmp=uicontrol('style','text','string',repmat('X',[1,50]),bf,...
    'units','normalized','visible','off');
fnp = get(tmp,'extent');
delete(tmp);
cw = 3*fnp(3)/50;

if strcmpi(typ,'image'),
    uicontrol(fg,...
        'style','edit',...
        'units','normalized',...
        'Position',posinpanel([0.61 0 0.37 .45],pcntp),...
        'Callback',@update_frames,...
        'tag','frame',...
        lf,...
        'BackgroundColor',col1,...
        'String',frames,'UserData',eval(frames));
% 'ForegroundGolor',col3,...
end;

% Help
uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',posinpanel([0.02 .5 cw .45],pcntp),...
    bf,...
    'Callback',@heelp,...
    'tag','?',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','?',...
    'ToolTipString','Show Help');

uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',posinpanel([0.03+cw .5 cw .45],pcntp),...
    bf,...
    'Callback',@editwin,...
    'tag','Ed',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Ed',...
    'ToolTipString','Edit Selected Files');

uicontrol(fg,...
   'Style','pushbutton',...
   'units','normalized',...
   'Position',posinpanel([0.04+2*cw .5 cw .45],pcntp),...
   bf,...
   'Callback',@select_rec,...
   'tag','Rec',...
   'ForegroundColor',col3,...
   'BackgroundColor',col1,...
   'String','Rec',...
   'ToolTipString','Recursively Select Files with Current Filter');

% Done
dne = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',posinpanel([0.05+3*cw .5 0.45-3*cw .45],pcntp),...
    bf,...
    'Callback',@(h,e)delete(h),...
    'tag','D',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Done',...
    'Enable','off',...
    'DeleteFcn',@null);

if numel(already)>=n(1) && numel(already)<=n(2),
    set(dne,'Enable','on');
end;

% Filter Button
uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',posinpanel([0.51 .5 0.1 .45],pcntp),...
    bf,...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'Callback',@clearfilt,...
    'String','Filt');

% Filter
uicontrol(fg,...
    'style','edit',...
    'units','normalized',...
    'Position',posinpanel([0.61 .5 0.37 .45],pcntp),...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    lf,...
    'Callback',@(ob,ev)update(ob),...
    'tag','regexp',...
    'String',filt,...
    'UserData',sfilt);

% Directories
db = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',posinpanel([0.02 0 0.47 1],pfdp),...
    lf,...
    'Callback',@click_dir_box,...
    'tag','dirs',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'Max',1,...
    'Min',0,...
    'String','',...
    'UserData',wd,...
    'Value',1);

% Files
tmp = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',posinpanel([0.51 0 0.47 1],pfdp),...
    lf,...
    'Callback',@click_file_box,...
    'tag','files',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'UserData',n,...
    'Max',10240,...
    'Min',0,...
    'String','',...
    'Value',1);
c0 = uicontextmenu('Parent',fg);
set(tmp,'uicontextmenu',c0);
uimenu('Label','Select All', 'Parent',c0,'Callback',@select_all);

% Drives
if strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64'),
    % get fh for lists
    tmp=uicontrol('style','text','string','X',lf,...
        'units','normalized','visible','off');
    fnp = get(tmp,'extent');
    delete(tmp);
    fh = 2*fnp(4); % Heuristics: why do we need 2*
    sz = get(db,'Position');
    sz(4) = sz(4)-fh-2*0.01;
    set(db,'Position',sz);
    uicontrol(fg,...
        'style','text',...
        'units','normalized',...
        'Position',posinpanel([0.02 1-fh-0.01 0.10 fh],pfdp),...
              'HorizontalAlignment','left',...
        lf,...
        'BackgroundColor',get(fg,'Color'),...
        'ForegroundColor',col3,...
        'String','Drive');
    uicontrol(fg,...
        'style','popupmenu',...
        'units','normalized',...
        'Position',posinpanel([0.12 1-fh-0.01 0.37 fh],pfdp),...
        lf,...
        'Callback',@setdrive,...
        'tag','drive',...
        'BackgroundColor',col1,...
        'ForegroundColor',col3,...
        'String',listdrives(false),...
        'Value',1);
end;

[pd,vl] = prevdirs([wd filesep]);

% Previous dirs
uicontrol(fg,...
    'style','popupmenu',...
    'units','normalized',...
    'Position',posinpanel([0.12 .05 0.86 .95*1/3],pdirp),...
    lf,...
    'Callback',@click_dir_list,...
    'tag','previous',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String',pd,...
    'Value',vl);
uicontrol(fg,...
    'style','text',...
    'units','normalized',...
    'Position',posinpanel([0.02 0 0.10 .95*1/3],pdirp),...
          'HorizontalAlignment','left',...
    lf,...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Prev');

% Parent dirs
uicontrol(fg,...
    'style','popupmenu',...
    'units','normalized',...
    'Position',posinpanel([0.12 1/3+.05 0.86 .95*1/3],pdirp),...
    lf,...
    'Callback',@click_dir_list,...
    'tag','pardirs',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String',pardirs(wd));
uicontrol(fg,...
    'style','text',...
    'units','normalized',...
    'Position',posinpanel([0.02 1/3 0.10 .95*1/3],pdirp),...
          'HorizontalAlignment','left',...
    lf,...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Up');

% Directory
uicontrol(fg,...
    'style','edit',...
    'units','normalized',...
    'Position',posinpanel([0.12 2/3 0.86 .95*1/3],pdirp),...
    lf,...
    'Callback',@edit_dir,...
    'tag','edit',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String','');
uicontrol(fg,...
    'style','text',...
    'units','normalized',...
    'Position',posinpanel([0.02 2/3 0.10 .95*1/3],pdirp),...
          'HorizontalAlignment','left',...
    lf,...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Dir');

resize_fun(fg);
set(fg, 'ResizeFcn',@resize_fun);
update(sel,wd)
set(fg,'windowstyle', 'modal');

waitfor(dne);
drawnow;
if ishandle(sel),
    t  = get(sel,'String');
    if isempty(t)
        t = {''};
    elseif sfilt.code == -1
        % canonicalise non-empty folder selection
        t = cellfun(@(t1)cpath(t1, pwd), t, 'UniformOutput',false);
    end;
    ok = 1;
end;
if ishandle(fg),  delete(fg); end;
drawnow;
return;
%=======================================================================

%=======================================================================
function apos = posinpanel(rpos,ppos)
% Compute absolute positions based on panel position and relative
% position
apos = [ppos(1:2)+ppos(3:4).*rpos(1:2) ppos(3:4).*rpos(3:4)];
%=======================================================================

%=======================================================================
function [pselp, pcntp, pfdp, pdirp] = panelpositions(fg, sellines)
if nargin == 1
    na = numel(get(findobj(fg,'Tag','selected'),'String'));
    n  = get(findobj(fg,'Tag','files'),'Userdata');
    sellines = min([max([n(2) na]), 4]);
end
lf = cfg_get_defaults('cfg_ui.lfont');
bf = cfg_get_defaults('cfg_ui.bfont');
% Create dummy text to estimate character height
t=uicontrol('style','text','string','Xg','units','normalized','visible','off',lf);
lfh = 1.05*get(t,'extent');
delete(t)
t=uicontrol('style','text','string','Xg','units','normalized','visible','off',bf);
bfh = 1.05*get(t,'extent');
delete(t)
% panel heights
% 3 lines for directory, parent and prev directory list
% variable height for dir/file navigation
% 2 lines for buttons, filter etc
% sellines plus scrollbar for selected files
pselh = sellines*lfh(4) + 1.2*lfh(4);
pselp = [0 0 1 pselh];
pcnth = 2*bfh(4);
pcntp = [0 pselh 1 pcnth];
pdirh = 3*lfh(4);
pdirp = [0 1-pdirh 1 pdirh];
pfdh  = 1-(pselh+pcnth+pdirh);
pfdp  = [0 pselh+pcnth 1 pfdh];

%=======================================================================

%=======================================================================
function null(varargin)
%=======================================================================

%=======================================================================
function omsg = msg(ob,str)
ob = sib(ob,'msg');
omsg = get(ob,'String');
set(ob,'String',str);
if nargin>=3,
    set(ob,'ForegroundColor',[1 0 0],'FontWeight','bold');
else
    set(ob,'ForegroundColor',[0 0 0],'FontWeight','normal');
end;
drawnow;
return;
%=======================================================================

%=======================================================================
function setdrive(ob,varargin)
st = get(ob,'String');
vl = get(ob,'Value');
update(ob,st{vl});
return;
%=======================================================================

%=======================================================================
function resize_fun(fg,varargin)
% do nothing
return;
[pselp pcntp pfdp pdirp] = panelpositions(fg);
set(findobj(fg,'Tag','msg'), 'Position',pselp);
set(findobj(fg,'Tag','pcnt'), 'Position',pcntp);
set(findobj(fg,'Tag','pfd'), 'Position',pfdp);
set(findobj(fg,'Tag','pdir'), 'Position',pdirp);
return;
%=======================================================================

%=======================================================================
function [d,mch] = prevdirs(d)
persistent pd
if ~iscell(pd), pd = {}; end;
if nargin == 0
    d = pd;
else
    if ~iscell(d)
        d = cellstr(d);
    end
    d   = unique(d(:));
    mch = cellfun(@(d1)find(strcmp(d1,pd)), d, 'UniformOutput',false);
    sel = cellfun(@isempty, mch);
    npd = numel(pd);
    pd  = [pd(:);d(sel)];
    mch = [mch{~sel} npd+(1:nnz(sel))];
    d = pd;
end
return;
%=======================================================================

%=======================================================================
function pd = pardirs(wd)
if ispc
    fs = '\\';
else
    fs = filesep;
end
pd1 = textscan(wd,'%s','delimiter',fs,'MultipleDelimsAsOne',1);
if ispc
    pd = cell(size(pd1{1}));
    pd{end} = pd1{1}{1};
    for k = 2:numel(pd1{1})
        pd{end-k+1} = fullfile(pd1{1}{1:k},filesep);
    end
else
    pd = cell(numel(pd1{1})+1,1);
    pd{end} = filesep;
    for k = 1:numel(pd1{1})
        pd{end-k} = fullfile(filesep,pd1{1}{1:k},filesep);
    end
end
%=======================================================================

%=======================================================================
function clearfilt(ob,varargin)
set(sib(ob,'regexp'),'String','.*');
update(ob);
return;
%=======================================================================

%=======================================================================
function click_dir_list(ob,varargin)
vl = get(ob,'Value');
ls = get(ob,'String');
update(ob,deblank(ls{vl}));
return;
%=======================================================================

%=======================================================================
function edit_dir(ob,varargin)
update(ob,get(ob,'String'));
return;
%=======================================================================

%=======================================================================
function c = get_current_char(lb)
fg = sib(lb, mfilename);
c = get(fg, 'CurrentCharacter');
if ~isempty(c)
    % reset CurrentCharacter
    set(fg, 'CurrentCharacter', char(13));
end
%=======================================================================

%=======================================================================
function click_dir_box(lb,varargin)
c = get_current_char(lb);
if isempty(c) || isequal(c,char(13))
    vl  = get(lb,'Value');
    str = get(lb,'String');
    pd  = get(sib(lb,'edit'),'String');
    while ~isempty(pd) && strcmp(pd(end),filesep)
        pd=pd(1:end-1);      % Remove any trailing fileseps
    end
    sel = str{vl};
    if strcmp(sel,'..'),     % Parent directory
        [dr odr] = fileparts(pd);
    elseif strcmp(sel,'.'),  % Current directory
        dr = pd;
        odr = '';
    else
        dr = fullfile(pd,sel);
        odr = '';
    end;
    update(lb,dr);
    if ~isempty(odr)
        % If moving up one level, try to set focus on previously visited
        % directory
        cdrs = get(lb, 'String');
        dind = find(strcmp(odr, cdrs));
        if ~isempty(dind)
            set(lb, 'Value',dind(1));
        end
    end
end
return;
%=======================================================================

%=======================================================================
function re = getfilt(ob)
ob  = sib(ob,'regexp');
ud  = get(ob,'UserData');
re  = struct('code',ud.code,...
             'frames',get(sib(ob,'frame'),'UserData'),...
             'ext',{ud.ext},...
             'filt',{{get(sib(ob,'regexp'),'String')}});
return;
%=======================================================================

%=======================================================================
function update(lb,dr)
lb = sib(lb,'dirs');
if nargin<2 || isempty(dr) || ~ischar(dr)
    dr = get(lb,'UserData');
end;
if ~(strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64'))
    dr    = [filesep dr filesep];
else
    dr    = [dr filesep];
end;
dr(strfind(dr,[filesep filesep])) = [];
[f,d] = listfiles(dr,getfilt(lb));
if isempty(d),
    dr    = get(lb,'UserData');
    [f,d] = listfiles(dr,getfilt(lb));
else
    set(lb,'UserData',dr);
end;
set(lb,'Value',1,'String',d);
set(sib(lb,'files'),'Value',1,'String',f);
set(sib(lb,'pardirs'),'String',pardirs(dr),'Value',1);
[ls,mch] = prevdirs(dr);
set(sib(lb,'previous'),'String',ls,'Value',mch);
set(sib(lb,'edit'),'String',dr);

if numel(dr)>1 && dr(2)==':',
    str = char(get(sib(lb,'drive'),'String'));
    mch = find(lower(str(:,1))==lower(dr(1)));
    if ~isempty(mch),
        set(sib(lb,'drive'),'Value',mch);
    end;
end;
return;
%=======================================================================

%=======================================================================
function update_frames(lb,varargin)
str = get(lb,'String');
%r   = get(lb,'UserData');
try
    r = eval(['[',str,']']);
catch
    msg(lb,['Failed to evaluate "' str '".'],'r');
    beep;
    return;
end;
if ~isnumeric(r),
    msg(lb,['Expression non-numeric "' str '".'],'r');
    beep;
else
    set(lb,'UserData',r);
    msg(lb,'');
    update(lb);
end;
%=======================================================================

%=======================================================================
function select_all(ob,varargin)
lb = sib(ob,'files');
set(lb,'Value',1:numel(get(lb,'String')));
drawnow;
click_file_box(lb);
return;
%=======================================================================

%=======================================================================
function click_file_box(lb,varargin)
c = get_current_char(lb);
if isempty(c) || isequal(c, char(13))
    vlo  = get(lb,'Value');
    if isempty(vlo),
        msg(lb,'Nothing selected');
        return;
    end;
    lim  = get(lb,'UserData');
    ob   = sib(lb,'selected');
    str3 = get(ob,'String');
    str  = get(lb,'String');
    lim1  = min([max([lim(2)-numel(str3),0]),numel(vlo)]);
    if lim1==0,
        msg(lb,['Selected ' num2str(size(str3,1)) '/' num2str(lim(2)) ' already.']);
        beep;
        set(sib(lb,'D'),'Enable','on');
        return;
    end;

    vl   = vlo(1:lim1);
    msk  = false(size(str,1),1);
    if vl>0, msk(vl) = true; else msk = []; end;
    str1 = str( msk);
    str2 = str(~msk);
    dr   = get(sib(lb,'edit'), 'String');
    str1 = strcat(dr, str1);

    set(lb,'Value',min([vl(1),numel(str2)]),'String',str2);
    r    = (1:numel(str1))+numel(str3);
    str3 = [str3(:);str1(:)];
    set(ob,'String',str3,'Value',r);
    if numel(vlo)>lim1,
        msg(lb,['Retained ' num2str(lim1) '/' num2str(numel(vlo))...
            ' of selection.']);
        beep;
    elseif isfinite(lim(2))
        if lim(1)==lim(2),
            msg(lb,['Selected ' num2str(numel(str3)) '/' num2str(lim(2)) ' files.']);
        else
            msg(lb,['Selected ' num2str(numel(str3)) '/' num2str(lim(1)) '-' num2str(lim(2)) ' files.']);
        end;
    else
        if size(str3,1) == 1, ss = ''; else ss = 's'; end;
        msg(lb,['Selected ' num2str(numel(str3)) ' file' ss '.']);
    end;
    if ~isfinite(lim(1)) || numel(str3)>=lim(1),
        set(sib(lb,'D'),'Enable','on');
    end;
end
return;
%=======================================================================

%=======================================================================
function obj = sib(ob,tag)
persistent fg;
if isempty(fg) || ~ishandle(fg)
    fg = findobj(0,'Tag',mfilename);
end
obj = findobj(fg,'Tag',tag);
return;
%if isempty(obj),
%    cfg_message('matlabbatch:usage',['Can''t find object with tag "' tag '".']);
%elseif length(obj)>1,
%    cfg_message('matlabbatch:usage',['Found ' num2str(length(obj)) ' objects with tag "' tag '".']);
%end;
%return;
%=======================================================================

%=======================================================================
function unselect(lb,varargin)
vl      = get(lb,'Value');
if isempty(vl), return; end;
str     = get(lb,'String');
msk     = true(numel(str),1);
if vl~=0, msk(vl) = false; end;
str2    = str(msk);
set(lb,'Value',min(vl(1),numel(str2)),'String',str2);
lim = get(sib(lb,'files'),'UserData');
if numel(str2)>= lim(1) && numel(str2)<= lim(2),
    set(sib(lb,'D'),'Enable','on');
else 
    set(sib(lb,'D'),'Enable','off');
end;

if numel(str2) == 1, ss1 = ''; else ss1 = 's'; end;
%msg(lb,[num2str(size(str2,1)) ' file' ss ' remaining.']);
if numel(vl) == 1, ss = ''; else ss = 's'; end;
msg(lb,['Unselected ' num2str(numel(vl)) ' file' ss '. ' ...
        num2str(numel(str2)) ' file' ss1 ' remaining.']);
return;
%=======================================================================

%=======================================================================
function unselect_all(ob,varargin)
lb = sib(ob,'selected');
set(lb,'Value',[],'String',{},'ListBoxTop',1);
msg(lb,'Unselected all files.');
lim = get(sib(lb,'files'),'UserData');
if lim(1)>0, set(sib(lb,'D'),'Enable','off'); end;
return;
%=======================================================================

%=======================================================================
function [f,d] = listfiles(dr,filt)
try
    ob = sib(gco,'msg');
    domsg = ~isempty(ob);
catch
    domsg = false;
end
if domsg
    omsg = msg(ob,'Listing directory...');
end
if nargin<2, filt = '';  end;
if nargin<1, dr   = '.'; end;
de      = dir(dr);
if ~isempty(de),
    d     = {de([de.isdir]).name};
    if ~any(strcmp(d, '.'))
        d = [{'.'}, d(:)'];
    end;
    if filt.code~=-1,
        f = {de(~[de.isdir]).name};
    else
        % f = d(3:end);
        f = d;
    end;
else
    d = {'.','..'};
    f = {};
end;

if domsg
    msg(ob,['Filtering ' num2str(numel(f)) ' files...']);
end
f  = do_filter(f,filt.ext);
f  = do_filter(f,filt.filt);
ii = cell(1,numel(f));
if filt.code==1 && (numel(filt.frames)~=1 || filt.frames(1)~=1),
    if domsg
        msg(ob,['Reading headers of ' num2str(numel(f)) ' images...']);
    end
    for i=1:numel(f),
        try
            ni = nifti(fullfile(dr,f{i}));
            dm = [ni.dat.dim 1 1 1 1 1];
            d4 = (1:dm(4))';
        catch
            d4 = 1;
        end;
        if all(isfinite(filt.frames))
            msk = false(size(filt.frames));
            for j=1:numel(msk), msk(j) = any(d4==filt.frames(j)); end;
            ii{i} = filt.frames(msk);
        else
            ii{i} = d4;
        end;
    end
elseif filt.code==1 && (numel(filt.frames)==1 && filt.frames(1)==1),
    for i=1:numel(f),
        ii{i} = 1;
    end;
end;

if domsg
    msg(ob,['Listing ' num2str(numel(f)) ' files...']);
end

[f,ind] = sortrows(f(:));
ii      = ii(ind);
msk     = true(1,numel(f));
for i=2:numel(f),
    if strcmp(f{i-1},f{i}),
        if filt.code==1,
            tmp      = sort([ii{i}(:) ; ii{i-1}(:)]);
            tmp(~diff(tmp,1)) = [];
            ii{i}    = tmp;
        end;
        msk(i-1) = false;
    end;
end;
f        = f(msk);
if filt.code==1,
    % Combine filename and frame number(s)
    ii       = ii(msk);
    nii      = cellfun(@numel, ii);
    c        = cell(sum(nii),1);
    fi       = cell(numel(f),1);
    for k = 1:numel(fi)
        fi{k} = k*ones(1,nii(k));
    end
    ii = [ii{:}];
    fi = [fi{:}];
    for i=1:numel(c),
        c{i} = sprintf('%s,%d', f{fi(i)}, ii(i));
    end;
    f        = c;
elseif filt.code==-1,
    fs = filesep;
    for i=1:numel(f),
        f{i} = [f{i} fs];
    end;
end;
f        = f(:);
d        = unique(d(:));
if domsg
    msg(ob,omsg);
end
return;
%=======================================================================

%=======================================================================
function [f,ind] = do_filter(f,filt)
t2 = false(numel(f),1);
filt_or = sprintf('(%s)|',filt{:});
t1 = regexp(f,filt_or(1:end-1));
if numel(f)==1 && ~iscell(t1), t1 = {t1}; end;
for i=1:numel(t1),
    t2(i) = ~isempty(t1{i});
end;
ind = find(t2);
f   = f(t2);
return;
%=======================================================================

%=======================================================================
function heelp(ob,varargin)
[col1,col2,col3,fn] = colours;
fg = sib(ob,mfilename);
t  = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.01 0.01 0.98 0.98],...
    fn,...
    'BackgroundColor',col2,...
    'ForegroundColor',col3,...
    'Max',0,...
    'Min',0,...
    'tag','HelpWin',...
    'String','                   ');
c0 = uicontextmenu('Parent',fg);
set(t,'uicontextmenu',c0);
uimenu('Label','Done', 'Parent',c0,'Callback',@helpclear);

str  = cfg_justify(t, {[...
'File Selection help. You can return to selecting files via the right mouse button (the "Done" option). '],...
'',[...
'The panel at the bottom shows files that are already selected. ',...
'Clicking a selected file will un-select it. To un-select several, you can ',...
'drag the cursor over the files, and they will be gone on release. ',...
'You can use the right mouse button to un-select everything.'],...
'',[...
'Directories are navigated by editing the name of the current directory (where it says "Dir"), ',...
'by going to one of the previously entered directories ("Prev"), ',...
'by going to one of the parent directories ("Up") or by navigating around ',...
'the parent or subdirectories listed in the left side panel.'],...
'',[...
'Files matching the filter ("Filt") are shown in the panel on the right. ',...
'These can be selected by clicking or dragging.  Use the right mouse button if ',...
'you would like to select all files.  Note that when selected, the files disappear ',...
'from this panel.  They can be made to reappear by re-specifying the directory ',...
'or the filter. '],...
'',[...
'Both directory and file lists can also be browsed by typing the leading ',...
'character(s) of a directory or file name.'],...
'',[...
'Note that the syntax of the filter differs from that used by most other file selectors. ',...
'The filter works using so called ''regular expressions''. Details can be found in the ',...
'MATLAB help text on ''regexp''. ',...
'The following is a list of symbols with special meaning for filtering the filenames:'],...
'    ^     start of string',...
'    $     end of string',...
'    .     any character',...
'    \     quote next character',...
'    *     match zero or more',...
'    +     match one or more',...
'    ?     match zero or one, or match minimally',...
'    {}    match a range of occurrances',...
'    []    set of characters',...
'    [^]   exclude a set of characters',...
'    ()    group subexpression',...
'    \w    match word [a-z_A-Z0-9]',...
'    \W    not a word [^a-z_A-Z0-9]',...
'    \d    match digit [0-9]',...
'    \D    not a digit [^0-9]',...
'    \s    match white space [ \t\r\n\f]',...
'    \S    not a white space [^ \t\r\n\f]',...
'    \<WORD\>    exact word match',...
'',[...
'Individual time frames of image files can also be selected.  The frame filter ',...
'allows specified frames to be shown, which is useful for image files that ',...
'contain multiple time points.  If your images are only single time point, then ',...
'reading all the image headers can be avoided by specifying a frame filter of "1". ',...
'The filter should contain a list of integers indicating the frames to be used. ',...
'This can be generated by e.g. "1:100", or "1:2:100".'],...
'',[...
'The recursive selection button (Rec) allows files matching the regular expression to ',...
'be recursively selected.  If there are many directories to search, then this can take ',...
'a while to run.'],...
'',[...
'There is also an edit button (Ed), which allows you to edit your selection of files. ',...
'When you are done, then use the menu-button of your mouse to either cancel or accept your changes.'],''});
set(t,'String',str);
return;
%=======================================================================

%=======================================================================
function helpclear(ob,varargin)
ob = get(ob,'Parent');
ob = get(ob,'Parent');
ob = findobj(ob,'Tag','HelpWin');
delete(ob);
%=======================================================================

%=======================================================================
function t = cpath(t,d)
switch filesep,
case '/',
    mch = '^/';
    fs  = '/';
    fs1 = '/';
case '\',
    mch = '^.:\\';
    fs  = '\';
    fs1 = '\\';
otherwise;
    cfg_message('matlabbatch:usage','What is this filesystem?');
end

if isempty(regexp(t,mch,'once')),
    if (nargin<2)||isempty(d), d = pwd; end;
    t = [d fs t];
end;

% Replace occurences of '/./' by '/' (problems with e.g. /././././././')
re = [fs1 '\.' fs1];
while ~isempty(regexp(t,re, 'once' )),
    t  = regexprep(t,re,fs);
end;
t  = regexprep(t,[fs1 '\.' '$'], fs);

% Replace occurences of '/abc/../' by '/'
re = [fs1 '[^' fs1 ']+' fs1 '\.\.' fs1];
while ~isempty(regexp(t,re, 'once' )),
    t  = regexprep(t,re,fs,'once');
end;
t  = regexprep(t,[fs1 '[^' fs1 ']+' fs1 '\.\.' '$'],fs,'once');

% Replace '//'
t  = regexprep(t,[fs1 '+'], fs);
%=======================================================================

%=======================================================================
function editwin(ob,varargin)
[col1,col2,col3,lf,bf] = colours;
fg   = gcbf;
lb   = sib(ob,'selected');
str  = get(lb,'String');
ac   = allchild(fg);
acv  = get(ac,'Visible');
h    = uicontrol(fg,...
    'Style','Edit',...
    'units','normalized',...
    'String',str,...
    lf,...
    'Max',2,...
    'Tag','EditWindow',...
    'HorizontalAlignment','Left',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'Position',[0.01 0.08 0.98 0.9],...
    'Userdata',struct('ac',{ac},'acv',{acv}));
ea   = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[.01 .01,.32,.07],...
    bf,...
    'Callback',@editdone,...
    'tag','EditWindowAccept',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Accept');
ee   = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[.34 .01,.32,.07],...
    bf,...
    'Callback',@editeval,...
    'tag','EditWindowEval',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Eval');
ec   = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[.67 .01,.32,.07],...
    bf,...
    'Callback',@editclear,...
    'tag','EditWindowCancel',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Cancel');
set(ac,'visible','off');
%=======================================================================

%=======================================================================
function editeval(ob,varargin)
ob  = get(ob, 'Parent');
ob  = sib(ob, 'EditWindow');
str = get(ob, 'String');
if ~isempty(str)
    [out,sts] = cfg_eval_valedit(char(str));
    if sts && (iscellstr(out) || ischar(out))
        set(ob, 'String', cellstr(out));
    else
        fgc = get(ob, 'ForegroundColor');
        set(ob, 'ForegroundColor', 'red');
        pause(1);
        set(ob, 'ForegroundColor', fgc);
    end
end
%=======================================================================

%=======================================================================
function editdone(ob,varargin)
ob  = get(ob,'Parent');
ob  = sib(ob,'EditWindow');
str = cellstr(get(ob,'String'));
if isempty(str) || isempty(str{1})
    str = {};
else
    dstr = deblank(str);
    if ~isequal(str, dstr)
        c = questdlg(['Some of the filenames contain trailing blanks. This may ' ...
            'be due to copy/paste of strings between MATLAB and the ' ...
            'edit window. Do you want to remove any trailing blanks?'], ...
            'Trailing Blanks in Filenames', ...
            'Remove', 'Keep', 'Remove');
        switch lower(c)
            case 'remove'
                str = dstr;
        end
    end
    filt = getfilt(ob);
    if filt.code >= 0 % filter files, but not dirs
        [p,n,e] = cellfun(@fileparts, str, 'uniformoutput',false);
        fstr = strcat(n, e);
        [fstr1,fsel] = do_filter(fstr, filt.ext);
        str = str(fsel);
    end
end
lim = get(sib(ob,'files'),'UserData');
if numel(str)>lim(2),
    msg(ob,['Retained ' num2str(lim(2)) ' of the ' num2str(numel(str)) ' files.']);
    beep;
    str = str(1:lim(2));
elseif isfinite(lim(2)),
    if lim(1)==lim(2),
        msg(ob,['Selected ' num2str(numel(str)) '/' num2str(lim(2)) ' files.']);
    else
        msg(ob,['Selected ' num2str(numel(str)) '/' num2str(lim(1)) '-' num2str(lim(2)) ' files.']);
    end;
else
    if numel(str) == 1, ss = ''; else ss = 's'; end;
    msg(ob,['Specified ' num2str(numel(str)) ' file' ss '.']);
end;
if ~isfinite(lim(1)) || numel(str)>=lim(1),
    set(sib(ob,'D'),'Enable','on');
else
    set(sib(ob,'D'),'Enable','off');
end;
set(sib(ob,'selected'),'String',str,'Value',[]);
acs = get(ob,'Userdata');
fg = gcbf;
delete(findobj(fg,'-regexp','Tag','^EditWindow.*'));
set(acs.ac,{'Visible'},acs.acv);
%=======================================================================

%=======================================================================
function editclear(ob,varargin)
fg = gcbf;
acs = get(findobj(fg,'Tag','EditWindow'),'Userdata');
delete(findobj(fg,'-regexp','Tag','^EditWindow.*'));
set(acs.ac,{'Visible'},acs.acv);
%=======================================================================

%=======================================================================
function [c1,c2,c3,lf,bf] = colours
c1 = [1 1 1];
c2 = [1 1 1];
c3 = [0 0 0];
lf = cfg_get_defaults('cfg_ui.lfont');
bf = cfg_get_defaults('cfg_ui.bfont');
if isempty(lf)
    lf = struct('FontName',get(0,'FixedWidthFontName'), ...
                'FontWeight','normal', ...
                'FontAngle','normal', ...
                'FontSize',14, ...
                'FontUnits','points');
end
if isempty(bf)
    bf = struct('FontName',get(0,'FixedWidthFontName'), ...
                'FontWeight','normal', ...
                'FontAngle','normal', ...
                'FontSize',14, ...
                'FontUnits','points');
end
%=======================================================================

%=======================================================================
function select_rec(ob, varargin)
start = get(sib(ob,'edit'),'String');
filt  = get(sib(ob,'regexp'),'Userdata');
filt.filt = {get(sib(ob,'regexp'), 'String')};
fob = sib(ob,'frame');
if ~isempty(fob)
    filt.frames = get(fob,'Userdata');
else
    filt.frames = [];
end;
ptr    = get(gcbf,'Pointer');
try
    set(gcbf,'Pointer','watch');
    sel    = select_rec1(start,filt);
catch
    set(gcbf,'Pointer',ptr);
    sel = {};
end;
set(gcbf,'Pointer',ptr);
already= get(sib(ob,'selected'),'String');
fb     = sib(ob,'files');
lim    = get(fb,'Userdata');
limsel = min(lim(2)-size(already,1),size(sel,1));
set(sib(ob,'selected'),'String',[already(:);sel(1:limsel)],'Value',[]);
msg(ob,sprintf('Added %d/%d matching files to selection.', limsel, size(sel,1)));
if ~isfinite(lim(1)) || size(sel,1)>=lim(1),
    set(sib(ob,'D'),'Enable','on');
else
    set(sib(ob,'D'),'Enable','off');
end;
%=======================================================================

%=======================================================================
function sel=select_rec1(cdir,filt)
sel={};
[t,d] = listfiles(cdir,filt);
if ~isempty(t)
    sel = strcat([cdir,filesep],t);
end;
for k = 1:numel(d)
    if ~any(strcmp(d{k},{'.','..'}))
        sel1 = select_rec1(fullfile(cdir,d{k}),filt);
        sel  = [sel(:); sel1(:)];
    end;
end;
%=======================================================================

%=======================================================================
function sfilt=mk_filter(typ,filt,frames)
if nargin<3, frames  = 1;     end;
if nargin<2, filt    = '.*';    end;
if nargin<1, typ     = 'any';   end;
switch lower(typ),
case {'any','*'}, code = 0; ext = {'.*'};
case {'image'},   code = 1; ext = {'.*\.nii(,\d+){0,2}$','.*\.img(,\d+){0,2}$','.*\.NII(,\d+){0,2}$','.*\.IMG(,\d+){0,2}$'};
case {'mesh'},    code = 0; ext = {'.*\.gii$','.*\.GII$','.*\.mat$','.*\.MAT$'};
case {'nifti'},   code = 0; ext = {'.*\.nii$','.*\.img$','.*\.NII$','.*\.IMG$'};
case {'gifti'},   code = 0; ext = {'.*\.gii$','.*\.GII$'};
case {'extimage'},   code = 1; ext = {'.*\.nii(,[0-9]*){0,2}$',...
                            '.*\.img(,[0-9]*){0,2}$',...
                            '.*\.NII(,[0-9]*){0,2}$',...
                            '.*\.IMG(,[0-9]*){0,2}$'};
case {'xml'},     code = 0; ext = {'.*\.xml$','.*\.XML$'};
case {'mat'},     code = 0; ext = {'.*\.mat$','.*\.MAT$','.*\.txt','.*\.TXT'};
case {'batch'},   code = 0; ext = {'.*\.mat$','.*\.MAT$','.*\.m$','.*\.M$','.*\.xml$','.*\.XML$'};
case {'dir'},     code =-1; ext = {'.*'};
case {'extdir'},     code =-1; ext = {['.*' filesep '$']};
otherwise,        code = 0; ext = {typ};
end;
sfilt = struct('code',code,'frames',frames,'ext',{ext},...
                             'filt',{{filt}});
%=======================================================================

%=======================================================================
function drivestr = listdrives(reread)
persistent mydrivestr;
if isempty(mydrivestr) || reread
    driveLett = strcat(cellstr(char(('C':'Z')')), ':');
    dsel = false(size(driveLett));
    for i=1:numel(driveLett)
        dsel(i) = exist([driveLett{i} '\'],'dir')~=0;
    end
    mydrivestr = driveLett(dsel);
end;
drivestr = mydrivestr;
%=======================================================================

%=======================================================================
