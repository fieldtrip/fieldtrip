function [t,sts] = cfg_getfile(varargin)
% File selector
% FORMAT [t,sts] = cfg_getfile(n,typ,mesg,sel,wd,filt,prms)
%     n    - Number of files
%            A single value or a range.  e.g.
%            1       - Select one file
%            Inf     - Select any number of files
%            [1 Inf] - Select 1 to Inf files
%            [0 1]   - select 0 or 1 files
%            [10 12] - select from 10 to 12 files
%     typ  - file type
%            'any'   - all files
%            'batch' - matlabbatch batch files (.m, .mat and XML)
%            'dir'   - select a directory. By default, hidden ('.xyz') and
%                      MATLAB class directories are not shown.
%            'mat'   - Matlab .mat files or .txt files (assumed to contain
%                      ASCII representation of a 2D-numeric array)
%            'xml'   - XML files
%            Other strings act as a filter to regexp.  This means
%            that e.g. DCM*.mat files should have a typ of '^DCM.*\.mat$'.
%            A combination of types can be specified as a cellstr list of
%            types. A file must match at least one of the specified types.
%     mesg - a prompt (default 'Select files...')
%     sel  - list of already selected files
%     wd   - Directory to start off in
%     filt - value for user-editable filter (default '.*'). Can be a
%            string or cellstr. If it is a cellstr, filter expressions
%            will be concatenated by '|' (regexp OR operator).
%     prms - Type specific parameters
%
%     t    - selected files
%     sts  - status (1 means OK, 0 means window quit)
%
% FORMAT [t,ind] = cfg_getfile('Filter',files,typ,filt,prms,...)
% filter the list of files (column cell array) in the same way as the
% GUI would do. The 'prms...' argument(s) will be passed to a typ
% specific filtering function, if available.
% When filtering directory names, the filt argument will be applied to the
% last directory in a path only.
% t returns the filtered list (cell array), ind an index array, such that 
% t = files(ind).
%
% FORMAT cpath = cfg_getfile('CPath',path,cwd)
% function to canonicalise paths: Prepends cwd to relative paths, processes
% '..' & '.' directories embedded in path.
% path     - cellstr containing path names
% cwd      - cellstr containing current working directory [default '.']
% cpath    - conditioned paths, in same format as input path argument
%
% FORMAT [files,dirs]=cfg_getfile('List'[,direc[,typ[,filt[,prms]]]])
% Returns files matching the filter (filt) and directories within direc
% direc    - directory to search. Defaults to pwd.
% typ      - file type
% filt     - additional filter to select files with (see regexp)
%            e.g. '^w.*\.img$' 
% files    - files matching 'typ' and 'filt' in directory 'direc'
% dirs     - subdirectories of 'direc'
% FORMAT [files,dirs]=cfg_getfile('FPList'[,direc[,typ[,filt[,prms]]]])
% As above, but returns files with full paths (i.e. prefixes direc to
% each).
% FORMAT [files,dirs]=cfg_getfile('FPListRec'[,direc[,typ[,filt[,prms]]]])
% As above, but returns files with full paths (i.e. prefixes direc to
% each) and searches through sub directories recursively.
%
% FORMAT cfg_getfile('PrevDirs',dir)
% Add directory dir to list of previous directories.
% FORMAT dirs=cfg_getfile('prevdirs')
% Retrieve list of previous directories.
%
% FORMAT cfg_getfile('DirFilters', filter_list)
% Specify a list of regular expressions to filter directory names. To show
% all directories, use {'.*'}. Default is {'^[^.@]'}, i.e. directory names
% starting with '.' or '@' will not be shown.
%
% FORMAT cfg_getfile('ListDrives'[, reread])
% On PCWIN(64) machines, list all available drive letters. If reread is
% true, refresh internally cached list of drive letters.
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
% $Id: cfg_getfile.m 6467 2015-06-03 14:47:41Z guillaume $

t = {};
sts = false;
if nargin > 0 && ischar(varargin{1})
    switch lower(varargin{1})
        case 'cpath'
            if nargin >= 2 && nargin <= 3
                if all(cellfun(@iscellstr,varargin(2:end)))
                    t = cpath(varargin{2:end});
                    sts = true;
                else
                    cfg_message('cfg_getfile:notcellstr','Inputs to %s(''%s'',...) must be cellstr.', mfilename, varargin{1});
                end
            else
                cfg_message('cfg_getfile:nargin','%s(''%s'',...) Wrong number of inputs.', mfilename, varargin{1});
            end
        case 'filter'
            t = varargin{2};
            if isempty(t) || (numel(t) == 1 && isempty(t{1}))
                sts = 1;
                return;
            end
            if ndims(t) > 2 || size(t,2) > 1
                cfg_message('cfg_getfile:notcolumncell','Input file lists to %s(''%s'',...) must be a column cellstr.', mfilename, varargin{1});
            end
            t = cpath(t);
            [unused,n,e] = cellfun(@fileparts,t,'UniformOutput',false);
            t1           = strcat(n,e);
            filt         = mk_filter(varargin{3:end});
            [unused,sts] = do_filter_exp(t1,filt);
            t            = t(sts);
        case {'list', 'fplist', 'fplistrec'}
            if nargin < 2
                direc = pwd;
            else
                direc = varargin{2};
            end
            if nargin < 3
                typ = 'any';
            else
                typ = varargin{3};
            end
            if nargin < 4
                filt = '.*';
            else
                filt = varargin{4};
            end
            if nargin < 5
                prms = {};
            else
                prms = varargin(5:end);
            end
            filt    = mk_filter(typ, filt, prms{:});
            direc   = cpath(cellstr(direc));
            t       = cell(numel(direc),1);
            sts     = cell(numel(direc),1);
            for k = 1:numel(t)
                if regexpi(varargin{1},'rec')
                    [t{k}, sts{k}] = select_rec1(direc{k}, filt, false);
                else
                    [t{k}, sts{k}] = listfiles(direc{k}, filt, false); % (sts is subdirs here)
                    sts{k} = sts{k}(~(strcmp(sts{k},'.')|strcmp(sts{k},'..'))); % remove '.' and '..' entries
                    if regexpi(varargin{1}, 'fplist') % return full pathnames
                        if ~isempty(t{k})
                            % starting with r2012b, this could be written as
                            % t{k} = fullfile(direc{k}, t{k});
                            t{k} = strcat(direc{k}, filesep, t{k});
                        end
                        if (nargout > 1) && ~isempty(sts{k})
                            sts{k} = cpath(sts{k}, direc(k));
                        end
                    end
                end
            end
            t   = cat(1,t{:});
            sts = cat(1,sts{:});
        case 'prevdirs',
            if nargin > 1
                prevdirs(varargin{2});
            end
            if nargout > 0 || nargin == 1
                t = prevdirs;
                sts = true;
            end
        case 'regfilter'
            t = reg_filter(varargin{2:end});
        case 'dirfilters'
            df       = reg_filter('dir');
            df.regex = varargin{2};
            df       = struct2cell(df);
            reg_filter(df{:});
        case 'listdrives'
            if ispc
                t = listdrives(varargin{2:end});
            else
                t = {};
            end
        otherwise
            cfg_message('matlabbatch:usage','Inappropriate usage.');
    end
else
    [t,sts] = selector(varargin{:});
end
%=======================================================================

%=======================================================================
function [t,ok] = selector(n,typ,mesg,already,wd,filt,varargin)
% Limits
if nargin<1 || ~isnumeric(n) || numel(n) > 2 
    n = [0 Inf];
else
    n(~isfinite(n)) = Inf;
    if numel(n)==1,     n    = [n n];    end
    if n(1)>n(2),       n    = n([2 1]); end
    if ~isfinite(n(1)), n(1) = 0;        end
end

% Filter
if nargin<2 || ~(ischar(typ) || iscellstr(typ)), typ     = 'any';   end

if nargin<6 || ~(ischar(filt) || iscellstr(filt)), filt    = '.*';    end
if iscellstr(filt) && numel(filt) > 1
    filt_or = sprintf('(%s)|',filt{:});
    filt_or = filt_or(1:end-1);
else
    filt_or = char(filt);
end
ok  = 0;

t = '';
sfilt = mk_filter(typ,filt,varargin{:});

% Message
if nargin<3 || ~(ischar(mesg) || iscellstr(mesg))
    mesg    = 'Select files...';
elseif iscellstr(mesg)
    mesg = char(mesg);
end

% Already selected files
if nargin<4 || isempty(already) || (iscell(already) && isempty(already{1}))
    already = {};
else
    % Canonicalise already selected paths
    already = cpath(already);
    % Filter already selected files by type, but not by user defined filter
    [pd, nam, ext] = cellfun(@(a1)fileparts(a1), already, ...
        'UniformOutput',false);
    sfilt1      = sfilt;
    sfilt1.filt = '.*';
    [unused, ind] = do_filter_exp(strcat(nam,ext),sfilt1);
    already = already(ind);
    % Add folders of already selected files to prevdirs list
    prevdirs(pd(ind));
end

% Working directory
if nargin<5 || isempty(wd) || ~(iscellstr(wd) || ischar(wd))
    if isempty(already)
        wd      = pwd;
    else
        wd = fileparts(already{1});
        if isempty(wd)
            wd = pwd;
        end
    end
end
wd = char(cpath(cellstr(wd)));

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
filtlines = sum(arrayfun(@(f)numel(f.prms.val), sfilt.tfilt));
[pselp, pfiltp, pcntp, pfdp, pdirp] = panelpositions(fg, sellines+1, filtlines);

% Messages are displayed as title of Selected Files uipanel
psel = uipanel(fg,...
    'units','normalized',...
    'Position',pselp,...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'BorderType','none',...
    'Tag','msg');

% Selected Files
sel = uicontrol(psel,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0 0 1 1],...
    lf{:},...
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

% Filters uipanel
if filtlines > 0
    pfilt = uipanel(fg,...
        'units','normalized',...
        'Position',pfiltp,...
        lf{:},...
        'BackgroundColor',get(fg,'Color'),...
        'ForegroundColor',col3,...
        'BorderType','none',...
        'Tag','pfilt');
    filth = 1/filtlines;
    cfiltp = 0;
    for cf = 1:numel(sfilt.tfilt)
        [id, stop, vals] = list(sfilt.tfilt(cf).prms,cfg_findspec({{'class','cfg_entry'}}),cfg_tropts(cfg_findspec,0,inf,0,inf,false),{'name','val'});
        for ci = 1:numel(id)
            % Label
            uicontrol(pfilt,...
                'style','text',...
                'units','normalized',...
                'Position',[0.02 cfiltp 0.47 filth],...
                'HorizontalAlignment','left',...
                lf{:},...
                'BackgroundColor',get(fg,'Color'),...
                'ForegroundColor',col3,...
                'String',vals{1}{ci});
            % Value
            uid.cf = cf;
            uid.id = id{ci};
            uicontrol(pfilt,...
                'style','edit',...
                'units','normalized',...
                'Position',[0.51 cfiltp 0.47 filth],...
                'ForegroundColor',col3,...
                'BackgroundColor',col1,...
                lf{:},...
                'Callback',@(ob,ev)update_filt(get(sib('dirs'),'Userdata')),...
                'String',char(gencode_rvalue(vals{2}{ci}{1})),...
                'UserData',uid);
            cfiltp = cfiltp+filth;
        end
    end
end

% Buttons uipanel
pcnt = uipanel(fg,...
    'units','normalized',...
    'Position',pcntp,...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'BorderType','none',...
    'Tag','pcnt');

% get cwidth for buttons
tmp=uicontrol('style','text','string',repmat('X',[1,50]),bf{:},...
    'units','normalized','visible','off');
fnp = get(tmp,'extent');
delete(tmp);
cw = 3*fnp(3)/50;

% Help
uicontrol(pcnt,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.02 0 cw 1],...
    bf{:},...
    'Callback',@heelp,...
    'tag','?',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','?',...
    'ToolTipString','Show Help');

uicontrol(pcnt,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.03+cw 0 cw 1],...
    bf{:},...
    'Callback',@editwin,...
    'tag','Ed',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Ed',...
    'ToolTipString','Edit Selected Files');

uicontrol(pcnt,...
   'Style','pushbutton',...
   'units','normalized',...
   'Position',[0.04+2*cw 0 cw 1],...
   bf{:},...
   'Callback',@select_rec,...
   'tag','Rec',...
   'ForegroundColor',col3,...
   'BackgroundColor',col1,...
   'String','Rec',...
   'ToolTipString','Recursively Select Files with Current Filter');

% Done
dne = uicontrol(pcnt,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.05+3*cw 0 0.44-3*cw 1],...
    bf{:},...
    'Callback',@(h,e)delete(h),...
    'tag','D',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Done',...
    'Enable','off',...
    'DeleteFcn',@null);

% Filter label
uicontrol(pcnt,...
    'style','text',...
    'units','normalized',...
    'Position',[0.51 0 0.08 1],...
    'HorizontalAlignment','left',...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Filter');

% Filter Button
uicontrol(pcnt,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.59 0 0.08 1],...
    bf{:},...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'Callback',@clearfilt,...
    'String','Reset');

% Filter
uicontrol(pcnt,...
    'style','edit',...
    'units','normalized',...
    'Position',[0.68 0 0.3 1],...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    lf{:},...
    'Callback',@(ob,ev)update(get(sib('dirs'),'Userdata')),...
    'tag','regexp',...
    'String',filt_or,...
    'UserData',sfilt);

% File/Dir uipanel
pfd = uipanel(fg,...
    'units','normalized',...
    'Position',pfdp,...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'BorderType','none',...
    'Tag','pfd');

% Directories
db = uicontrol(pfd,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.02 0 0.47 1],...
    lf{:},...
    'Callback',@click_dir_box,...
    'tag','dirs',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'Max',1,...
    'Min',0,...
    'String','',...
    'ToolTipString','Navigate Directories',...
    'UserData',wd,...
    'Value',1);

% Files
tmp = uicontrol(pfd,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.51 0 0.47 1],...
    lf{:},...
    'Callback',@click_file_box,...
    'tag','files',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'UserData',n,...
    'Max',10240,...
    'Min',0,...
    'String','',...
    'ToolTipString','Select Items',...
    'Value',1);
c0 = uicontextmenu('Parent',fg);
set(tmp,'uicontextmenu',c0);
uimenu('Label','Select All', 'Parent',c0,'Callback',@select_all);

updatedir_fun = @(ob,ev)(update(char(cpath(subsref(get(ob,'String'),substruct('()',{get(ob,'Value')}))))));

% Drives
if strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64'),
    % get fh for lists
    tmp=uicontrol('style','text','string','X',lf{:},...
        'units','normalized','visible','off');
    fnp = get(tmp,'extent');
    delete(tmp);
    fh = 2*fnp(4); % Heuristics: why do we need 2*
    sz = get(db,'Position');
    sz(4) = sz(4)-fh-0.05;
    set(db,'Position',sz);
    uicontrol(pfd,...
        'style','text',...
        'units','normalized',...
        'Position',[0.02 1-fh-0.01 0.10 fh],...
              'HorizontalAlignment','left',...
        lf{:},...
        'BackgroundColor',get(fg,'Color'),...
        'ForegroundColor',col3,...
        'String','Drive');
    uicontrol(pfd,...
        'style','popupmenu',...
        'units','normalized',...
        'Position',[0.12 1-fh-0.01 0.37 fh],...
        lf{:},...
        'Callback',updatedir_fun,...
        'tag','drive',...
        'BackgroundColor',col1,...
        'ForegroundColor',col3,...
        'String',listdrives(false),...
        'Value',1);
end

% Directories uipanel
pdir = uipanel(fg,...
    'units','normalized',...
    'Position',pdirp,...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'BorderType','none',...
    'Tag','pdir');

[pd,vl] = prevdirs(wd);

% Previous dirs
uicontrol(pdir,...
    'style','popupmenu',...
    'units','normalized',...
    'Position',[0.12 0.05 0.86 .95*1/3],...
    lf{:},...
    'Callback',updatedir_fun,...
    'tag','previous',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String',pd,...
    'Value',vl);
uicontrol(pdir,...
    'style','text',...
    'units','normalized',...
    'Position',[0.02 0 0.10 .95*1/3],...
          'HorizontalAlignment','left',...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Prev');

% Parent dirs
uicontrol(pdir,...
    'style','popupmenu',...
    'units','normalized',...
    'Position',[0.12 1/3+.05 0.86 .95*1/3],...
    lf{:},...
    'Callback',updatedir_fun,...
    'tag','pardirs',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String',pardirs(wd));
uicontrol(pdir,...
    'style','text',...
    'units','normalized',...
    'Position',[0.02 1/3 0.10 .95*1/3],...
          'HorizontalAlignment','left',...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Up');

% Directory
uicontrol(pdir,...
    'style','edit',...
    'units','normalized',...
    'Position',[0.12 2/3 0.86 .95*1/3],...
    lf{:},...
    'Callback',@(ob,ev)(update(char(cpath({get(ob,'String')})))),...
    'tag','edit',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String','');
uicontrol(pdir,...
    'style','text',...
    'units','normalized',...
    'Position',[0.02 2/3 0.10 .95*1/3],...
          'HorizontalAlignment','left',...
    lf{:},...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Dir');

resize_fun(fg);
set(fg, 'ResizeFcn',@resize_fun);
update(wd)
checkdone('Initial selection.');
set(fg,'windowstyle', 'modal');

waitfor(dne);
drawnow;
if ishandle(sel),
    t  = get(sel,'String');
    if isempty(t)
        t = {''};
    elseif any(strcmp({sfilt.tfilt.typ},'dir'))
        % canonicalise non-empty folder selection
        t = cpath(t, {pwd});
    end
    ok = 1;
end
if ishandle(fg),  delete(fg); end
drawnow;
return;
%=======================================================================

%=======================================================================
function drivestr = listdrives(reread)
persistent mydrivestr;
if isempty(mydrivestr) || (nargin > 0 && reread)
    driveLett  = strcat(cellstr(char(('C':'Z')')), ':');
    dsel       = cellfun(@(dl)exist([dl '\'],'dir'),driveLett)~=0;
    mydrivestr = driveLett(dsel);
end
drivestr = mydrivestr;
%=======================================================================

%=======================================================================
function [d,mch] = prevdirs(d)
persistent pd
if ~iscell(pd), pd = {}; end
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
pd1 = pathparts(cellstr(wd));
if ispc
    pd = cell(numel(pd1{1}),1);
    pd{end} = pd1{1}{1};
    for k = 2:numel(pd1{1})
        pd{end-k+1} = fullfile(pd1{1}{1:k});
    end
else
    pd = cell(numel(pd1{1})+1,1);
    pd{end} = filesep;
    for k = 1:numel(pd1{1})
        pd{end-k} = fullfile(filesep,pd1{1}{1:k});
    end
end
%=======================================================================

%=======================================================================
function click_dir_box(lb,varargin)
c = get_current_char;
if isempty(c) || isequal(c,char(13))
    vl  = get(lb,'Value');
    str = get(lb,'String');
    pd  = get(sib('edit'),'String');
    sel = str{vl};
    if strcmp(sel,'..'),     % Parent directory
        [dr, odr] = fileparts(pd);
    elseif strcmp(sel,'.'),  % Current directory
        dr = pd;
        odr = '';
    else
        dr = fullfile(pd,sel);
        odr = '';
    end
    update(dr);
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
function update(dr)
lb = sib('dirs');
if nargin<1
    dr = get(lb,'UserData');
end
[f,d] = listfiles(dr,getfilt);
if isempty(d),
    dr    = get(lb,'UserData');
    [f,d] = listfiles(dr,getfilt);
else
    set(lb,'UserData',dr);
end
set(lb,'Value',1,'String',d);
set(sib('files'),'Value',1,'String',f);
set(sib('pardirs'),'String',pardirs(dr),'Value',1);
[ls,mch] = prevdirs(dr);
set(sib('previous'),'String',ls,'Value',mch);
set(sib('edit'),'String',dr);

if ispc && numel(dr)>1 && dr(2)==':',
    str = char(get(sib('drive'),'String'));
    mch = find(lower(str(:,1))==lower(dr(1)));
    if isempty(mch),
        str = listdrives(true);
        cstr = char(str);
        mch = find(lower(cstr(:,1))==lower(dr(1)));
        if ~isempty(mch)
            set(sib('drive'),'String',str, 'Value',mch);
        end
    else
        set(sib('drive'),'Value',mch);
    end
end
return;
%=======================================================================

%=======================================================================
function update_filt(dr)
uid  = get(gcbo, 'Userdata');
filt = get(sib('regexp'), 'Userdata');
valstr = get(gcbo, 'String');
[val, sts] = cfg_eval_valedit(valstr);
citem = subsref(filt.tfilt(uid.cf).prms, uid.id);
if ~sts
    % Try with quotes/brackets
    if citem.strtype == 's'
        valstr1 = sprintf('''%s''', valstr);
    else
        valstr1 = sprintf('[%s]', valstr);
    end
    [val, sts] = cfg_eval_valedit(valstr1);
end
if sts
    citem.val{1} = val;
    filt.tfilt(uid.cf).prms = subsasgn(filt.tfilt(uid.cf).prms, uid.id, citem);
    set(sib('regexp'), 'Userdata',filt);
end
set(gcbo, 'String',char(gencode_rvalue(citem.val{1})));
update(dr);
%=======================================================================

%=======================================================================
function checkdone(selmsg)
nsel = numel(get(sib('selected'),'String'));
fb   = sib('files');
lim  = get(fb,'Userdata');
if nsel == 1, s = ''; else s = 's'; end
if isequal(lim,[0 inf])
    msg('Selected %d file%s. (%s)',nsel,s,selmsg);
elseif lim(2) == inf
    msg('Selected %d/[%d-...] file%s. (%s)',nsel,lim(1),s,selmsg);
else
    msg('Selected %d/[%d-%d] file%s. (%s)',nsel,lim,s,selmsg);
end
if nsel>=lim(1) && (~isfinite(lim(2)) || nsel<=lim(2))
    set(sib('D'),'Enable','on');
else
    set(sib('D'),'Enable','off');
end
%=======================================================================

%=======================================================================
function select_all(varargin)
lb = sib('files');
set(lb,'Value',1:numel(get(lb,'String')));
drawnow;
click_file_box(lb);
return;
%=======================================================================

%=======================================================================
function click_file_box(lb,varargin)
c = get_current_char;
if isempty(c) || isequal(c, char(13))
    vlo  = get(lb,'Value');
    if isempty(vlo),
        msg('Nothing selected');
        return;
    end
    str  = get(lb,'String');
    dr   = get(sib('edit'), 'String');
    nsel = select(cellfun(@(str1)fullfile(dr,str1),str(vlo),'UniformOutput',false),false);
    vrem = setdiff(1:numel(str), vlo(1:nsel));
    set(lb,'Value',min([vlo(1),numel(str)-nsel]),'String',str(vrem));
end
return;
%=======================================================================

%=======================================================================
function nsel = select(new,repl)
if repl
    already = {};
else
    already = get(sib('selected'),'String');
end
fb       = sib('files');
lim      = get(fb,'Userdata');
nnew     = numel(new);
nalready = numel(already);
nsel     = min(lim(2)-nalready,nnew);
set(sib('selected'),'String',[already(:);new(1:nsel)],'Value',[]);
if nsel == 1, s = ''; else s = 's'; end
checkdone(sprintf('Added %d/%d file%s.',nsel,nnew,s));
return;
%=======================================================================

%=======================================================================
function select_rec(varargin)
start = get(sib('edit'),'String');
filt  = get(sib('regexp'),'Userdata');
filt.filt = {get(sib('regexp'), 'String')};
ptr    = get(gcbf,'Pointer');
try
    set(gcbf,'Pointer','watch');
    sel    = select_rec1(start,filt,true);
    select(sel,false);
end
set(gcbf,'Pointer',ptr);
%=======================================================================

%=======================================================================
function [f,d]=select_rec1(cdir,filt,cdflag)
[f,d] = listfiles(cdir,filt,cdflag);
if isempty(f)
    f = {};
else
    f = strcat(cdir,filesep,f);
end
dsel = cellfun(@(d1)any(strcmp(d1,{'.','..'})),d);
if any(~dsel)
    d       = strcat(cdir,filesep,d(~dsel));
    [f1,d1] = cellfun(@(d1)select_rec1(d1,filt,cdflag),d,'UniformOutput',false);
    f       = [f(:);cat(1,f1{:})];
    d       = [d(:);cat(1,d1{:})];
else
    d       = {};
end
%=======================================================================

%=======================================================================
function unselect(lb,varargin)
vl      = get(lb,'Value');
if isempty(vl)||(numel(vl)==1 && vl==0), return; end
str     = get(lb,'String');
msk     = true(numel(str),1);
msk(vl) = false;
str2    = str(msk);
set(lb,'Value',min(vl(1),numel(str2)),'String',str2);
if nnz(~msk) == 1, s = ''; else s = 's'; end
checkdone(sprintf('Unselected %d file%s.',numel(vl),s));
return;
%=======================================================================

%=======================================================================
function unselect_all(varargin)
lb = sib('selected');
set(lb,'Value',[],'String',{},'ListBoxTop',1);
checkdone('Unselected all files.');
return;
%=======================================================================

%=======================================================================
function [f,d] = listfiles(dr,filt,cdflag)
try
    ob = sib('msg');
    domsg = ~isempty(ob);
catch
    domsg = false;
end
if domsg
    omsg = msg('Listing directory...');
end
if nargin<3, cdflag = true; end
if nargin<2, filt = '';  end
if nargin<1, dr   = '.'; end
de      = dir(dr);
if isempty(de),
    d = {'..'};
    f = {};
    if domsg
        msg(omsg);
    end
    return
end

d = sort({de([de.isdir]).name})';
d = d(~strcmp(d,'.'));
dfilt = mk_filter('dir');
d = do_filter(d,dfilt.tfilt.regex);
% add back '..'
if ~any(strcmp('..',d))
    d = [{'..'}; d(:)];
end
dsel = strcmp({filt.tfilt.typ},'dir');
if any(dsel)
    fd = d(~strcmp('..',d));
    if cdflag && ~any(strcmp('.',fd))
        fd = [{'.'}; fd(:)];
    end
    % reset regex filter
    filt.tfilt(dsel).regex = {'.*'};
else
    fd = {};
end
if any(~dsel)
    ff = {de(~[de.isdir]).name}';
else
    ff = {};
end
f = [sort(fd(:)); sort(ff(:))];

if ~isempty(f)
    if domsg
        msg('Filtering %d files...',numel(f));
    end
    f  = sort(do_filter(f,filt.filt));
    if domsg
        msg('Listing %d files...',numel(f));
    end
    f1   = cell(size(filt.tfilt));
    i1   = cell(size(filt.tfilt));
    for k = 1:numel(filt.tfilt)
        [f11,i11] = do_filter(f,filt.tfilt(k).regex);
        if isempty(f11)||isempty(filt.tfilt(k).fun)
            f1{k} = f11;
            i1{k} = i11;
        else
            [unused,prms] = harvest(filt.tfilt(k).prms, filt.tfilt(k).prms, false, false);
            [f1{k},i12] = filt.tfilt(k).fun('list',dr,f11,prms);
            i1{k}       = i11(i12);
        end
    end
    % files might have been matched by multiple filters. Sort into roughly
    % alphabetical order before removing duplicates with 'stable' option.
    f = cat(1,f1{:});
    o = cat(1,i1{:});
    % [un,so] = sort(o);
    % f       = unique(f(so), 'stable');
    [un,fi,fj] = unique(f);
    ufi        = unique(fi(fj));
    [un,so]    = sort(o(ufi));
    f = f(ufi(so));
end
d = unique(d(:));
if domsg
    msg(omsg);
end
return;
%=======================================================================

%=======================================================================
function pp = pathparts(p)
% parse paths in cellstr p
% returns cell array of path component cellstr arrays
% For PC (WIN) targets, both '\' and '/' are accepted as filesep, similar
% to MATLAB fileparts
if ispc
    fs = '\\/';
else
    fs = filesep;
end
pp = cellfun(@(p1)textscan(p1,'%s','delimiter',fs,'MultipleDelimsAsOne',1),p);
if ispc
    for k = 1:numel(pp)
        if ~isempty(regexp(pp{k}{1}, '^[a-zA-Z]:$', 'once'))
            pp{k}{1} = strcat(pp{k}{1}, filesep);
        elseif ~isempty(regexp(p{k}, '^\\\\', 'once'))
            pp{k}{1} = strcat(filesep, filesep, pp{k}{1});
        end
    end
end
%=======================================================================

%=======================================================================
function t = cpath(t,d)
% canonicalise paths to full path names, removing xxx/./yyy and xxx/../yyy
% constructs
% t must be a cell array of (relative or absolute) paths, d must be a
% single cell containing the base path of relative paths in t
if ispc % valid absolute paths
    % Allow drive letter or UNC path
    mch = '^([a-zA-Z]:)|(\\\\[^\\]*)';
else
    mch = '^/';
end
if (nargin<2)||isempty(d), d = {pwd}; end
% Find partial paths, prepend them with d
ppsel    = cellfun(@isempty, regexp(t,mch,'once'));
t(ppsel) = cellfun(@(t1)fullfile(d{1},t1),t(ppsel),'UniformOutput',false);
% Break paths into cell lists of folder names
pt = pathparts(t);
% Remove single '.' folder names
sd = cellfun(@(pt1)strcmp(pt1,'.'),pt,'UniformOutput',false);
for cp = 1:numel(pt)
    pt{cp} = pt{cp}(~sd{cp});
end
% Go up one level for '..' folders, don't remove drive letter/server name
% from PC path
if ispc
    ptstart = 2;
else
    ptstart = 1;
end
for cp = 1:numel(pt)
    tmppt = {};
    for cdir = ptstart:numel(pt{cp})
        if strcmp(pt{cp}{cdir},'..')
            tmppt = tmppt(1:end-1);
        else
            tmppt{end+1} = pt{cp}{cdir};
        end
    end
    if ispc
        pt{cp} = [pt{cp}(1) tmppt];
    else
        pt{cp} = tmppt;
    end
end
% Assemble paths
if ispc
    t         = cellfun(@(pt1)fullfile(pt1{:}),pt,'UniformOutput',false);
else
    t         = cellfun(@(pt1)fullfile(filesep,pt1{:}),pt,'UniformOutput',false);
end
%=======================================================================

%=======================================================================
function clearfilt(varargin)
set(sib('regexp'),'String','.*');
update;
return;
%=======================================================================

%=======================================================================
function re = getfilt
% Get filter parameters from GUI widgets.
ob  = sib('regexp');
re  = get(ob,'UserData');
re.filt = {get(ob,'String')};
return;
%=======================================================================

%=======================================================================
function [f,ind] = do_filter_exp(f,filt)
% Filter a cellstr list of filenames by typ/regexp and filter handlers.
% FORMAT [f,ind] = do_filter_exp(f,filt)
% f    - cellstr list of filenames
% filt - a filter struct as returned by mk_filter
% ind  - index into f_input such that f_output = f_input(ind)
if (isempty(filt.filt) || any(strcmp(filt.filt,'.*'))) 
    ind0     = (1:numel(f))';
else
    [f,ind0] = do_filter(f,filt.filt);
end
ind1 = cell(size(filt.tfilt));
for k = 1:numel(filt.tfilt)
    [f1, ind1{k}] = do_filter(f,filt.tfilt(k).regex);
    if ~(isempty(f1) || isempty(filt.tfilt(k).fun))
        [unused,prms]         = harvest(filt.tfilt(k).prms, filt.tfilt(k).prms, false, false);
        [unused,ind2]         = filt.tfilt(k).fun('filter',f1,prms);
        ind1{k}               = ind1{k}(ind2);
    end
end
sel = unique(cat(1,ind1{:}));
ind = ind0(sel);
f   = f(sel);
return
%=======================================================================

%=======================================================================
function [f,ind] = do_filter(f,filt)
% Filter a cellstr list of filenames with a regular expression filter.
% FORMAT [f,ind] = do_filter(f,filt)
% f    - cellstr list of filenames
% filt - cellstr list of filter expressions (these will be ORed together).
% ind  - index into f_input such that f_output = f_input(ind)
filt_or = sprintf('(%s)|',filt{:});
t1      = regexp(f,filt_or(1:end-1));
if numel(f)==1 && ~iscell(t1), t1 = {t1}; end
ind = find(~cellfun(@isempty,t1));
f   = f(ind);
return;
%=======================================================================

%=======================================================================
function sfilt=mk_filter(typ,filt,varargin)
% Create a filter structure for a specific typ and set of parameters. For
% each typ, parameters should be a struct of parameter values with
% fieldnames as defined by the tags of the registered cfg_entry items.
if nargin<2, filt    = '.*';    end
if nargin<1, typ     = 'any';   end
sfilt.tfilt = reg_filter(typ);
for k = 1:numel(sfilt.tfilt)
    if nargin >= k+2
        sfilt.tfilt(k).prms = initialise(sfilt.tfilt(k).prms, varargin{k}, false);
    end
end
sfilt.filt  = cellstr(filt);
return
%=======================================================================

%=======================================================================
function filt = reg_filter(typ,filt,cflag,fun,prms)
% Store and retrieve filter definitions
% FORMAT reg_filter(typ,filt,[fun[,prms]])
% typ       - char identifying the filter
% filt      - filter expression (cellstr containing regular expressions).
%             This filter must cover both real filenames to be matched by
%             this type as well as expanded filenames as returned by the
%             associated handler function (if any).
% cflag     - case sensitivity - true for case sensitive matching, false
%             for case insensitive matching (default).
% fun       - optional function to expand typ specific file lists. Calling
%             syntax is 
%             files       = fun('list',dir,files,prms) 
%             [files ind] = fun('filter',files,val1,prms)
%             dir   - directory
%             files - cellstr of filenames (in 'list' mode: relative to
%                     dir)
%             prms  - optional parameters - a struct with fields with
%                     tagnames according to registered prms cfg_entry
%                     objects
%             ind   - 'filter' mode: index vector into the output list of
%                     files such that ind(k) = l if files_out{k} is
%                     expanded from files_in{l}. 
%             ind   - 'list' mode: index vector into the input list of
%                     files such that files_out = files_in(ind). 
%             In 'list' mode, the handler may modify the list of
%             files (delete non-matching files, add indices to filenames
%             etc). In 'filter' mode, the handler may delete files from the
%             list, but must not modify filenames.
% prms      - optional cell array of cfg_entry items, describing parameter
%             inputs for handler function
% FORMAT rfilt = reg_filter(typ)
% typ     - char or cellstr list of registered types
% filt    - list of matching filters.
persistent lfilt;
if isempty(lfilt)
    dtypes = {'any'  ,'xml'     ,'mat'              ,'batch'                   ,'dir'};
    dregex = {{'.*'} ,{'\.xml$'},{'\.mat$','\.txt$'},{'\.mat$','\.m$','\.xml$'},{'^[^.@]'}};
    dcflag = false;
    dfun   = {''     ,''        ,''                 ,''                        ,@local_isdir};
    dprms  = {cfg_branch};
    lfilt = struct('typ',dtypes, 'regex',dregex, 'cflag',dcflag, 'fun',dfun, 'prms',dprms);
end
if nargin == 0
    filt = lfilt;
elseif nargin == 1
    filt  = struct('typ',{}, 'regex',{}, 'cflag',{}, 'fun',{}, 'prms',{});
    ctyp  = cellstr(typ);
    fsel  = any(cell2mat(cellfun(@(t)strcmpi({lfilt.typ},t),ctyp,'UniformOutput',false)'),1);
    tsel  = any(cell2mat(cellfun(@(t)strcmpi(ctyp,t),{lfilt.typ},'UniformOutput',false)'),1);
    if any(fsel)
        filt = lfilt(fsel);
        % sort filters into typ order
        [unused, st] = sort(ctyp(tsel));
        [unused, sf] = sort({filt.typ});
        filt = filt(st(sf));
    end
    if any(~tsel)
        filt(end+1) = struct('typ','ad-hoc', 'regex',{ctyp(~tsel)}, 'cflag',true, 'fun','', 'prms',cfg_branch);
    end
else
    nfilt = struct('typ',typ, 'regex',{cellstr(filt)}, 'cflag',false, 'fun','', 'prms',cfg_branch);
    if nargin > 2
        nfilt.cflag = cflag;
    end
    if nargin > 3
        nfilt.fun   = fun;
    end
    if nargin > 4
        if isa(prms, 'cfg_branch')
            nfilt.prms = prms;
        else
            nfilt.prms.val  = prms;
        end
    end
    sel  = strcmpi({lfilt.typ},typ);
    if any(sel)
        lfilt(sel) = nfilt;
    else
        lfilt = [lfilt nfilt];
    end
end
return
%=======================================================================

%=======================================================================
function editwin(varargin)
[col1,col2,col3,lf,bf] = colours;
fg   = sib(mfilename);
lb   = sib('selected');
str  = get(lb,'String');
ac   = allchild(fg);
acv  = get(ac,'Visible');
h    = uicontrol(fg,...
    'Style','Edit',...
    'units','normalized',...
    'String',str,...
    lf{:},...
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
    bf{:},...
    'Callback',@editdone,...
    'tag','EditWindowAccept',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Accept');
ee   = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[.34 .01,.32,.07],...
    bf{:},...
    'Callback',@editeval,...
    'tag','EditWindowEval',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Eval');
ec   = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[.67 .01,.32,.07],...
    bf{:},...
    'Callback',@editclear,...
    'tag','EditWindowCancel',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Cancel');
set(ac,'visible','off');
%=======================================================================

%=======================================================================
function editeval(varargin)
ob  = sib('EditWindow');
str = get(ob, 'String');
if isempty(str), return, end
if ~isempty(str)
    [out, sts] = cfg_eval_valedit(char(str));
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
function editdone(varargin)
ob  = sib('EditWindow');
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
    filt = getfilt;
    filt.filt = {};
    [p, n, e] = cellfun(@fileparts, str, 'uniformoutput',false);
    fstr = strcat(n, e);
    [fstr1, fsel] = do_filter_exp(fstr, filt);
    str = str(fsel);
end
select(str,true);
editclear;
%=======================================================================

%=======================================================================
function editclear(varargin)
acs = get(sib('EditWindow'),'Userdata');
delete(findobj(sib(mfilename),'-regexp','Tag','^EditWindow.*'));
set(acs.ac,{'Visible'},acs.acv);
%=======================================================================

%=======================================================================
function heelp(varargin)
[col1,col2,col3,lf,bf] = colours;
fg = sib(mfilename);
ac   = allchild(fg);
acv  = get(ac,'Visible');
set(ac,'Visible','off');
t  = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.01 0.1 0.98 0.9],...
    lf{:},...
    'BackgroundColor',col2,...
    'ForegroundColor',col3,...
    'Max',0,...
    'Min',0,...
    'tag','HelpWin',...
    'String','                   ',...
    'Userdata',struct('ac',{ac},'acv',{acv}));
ea   = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[.01 .01,.99,.07],...
    bf{:},...
    'Callback',@helpclear,...
    'tag','HelpWinClear',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Close');

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
acs = get(sib('HelpWin'),'Userdata');
delete(findobj(sib(mfilename),'-regexp','Tag','^HelpWin.*'));
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
    lf = {'FontName',get(0,'FixedWidthFontName'), ...
                'FontWeight','normal', ...
                'FontAngle','normal', ...
                'FontSize',14, ...
                'FontUnits','points'};
end
if isempty(bf)
    bf = {'FontName',get(0,'FixedWidthFontName'), ...
                'FontWeight','normal', ...
                'FontAngle','normal', ...
                'FontSize',14, ...
                'FontUnits','points'};
end
%=======================================================================

%=======================================================================
function [pselp, pfiltp, pcntp, pfdp, pdirp] = panelpositions(fg, sellines, filtlines)
if nargin == 1
    na = numel(get(findobj(fg,'Tag','selected'),'String'));
    n  = get(findobj(fg,'Tag','files'),'Userdata');
    sellines = min([max([n(2) na]), 4]);
    sfilt = get(findobj(fg,'Tag','regexp'),'Userdata');
    filtlines = sum(arrayfun(@(f)numel(f.prms.val), sfilt.tfilt));
end
lf = cfg_get_defaults('cfg_ui.lfont');
bf = cfg_get_defaults('cfg_ui.bfont');
% Create dummy text to estimate character height
t=uicontrol('style','text','string','Xg','units','normalized','visible','off',lf{:});
lfh = 1.05*get(t,'extent');
delete(t)
t=uicontrol('style','text','string','Xg','units','normalized','visible','off',bf{:});
bfh = 1.05*get(t,'extent');
delete(t)
% panel heights
% 3 lines for directory, parent and prev directory list
% variable height for dir/file navigation
% 1 line for buttons and regexp filter
% filtlines for extra filter input
% sellines plus scrollbar for selected files
% build up from bottom to top
pselh  = sellines*lfh(4) + 1.2*lfh(4);
pselp  = [0 0 1 pselh];
pfilth = filtlines*bfh(4);
pfiltp = [0 pselh 1 pfilth];
pcnth  = 1*bfh(4);
pcntp  = [0 pselh+pfilth 1 pcnth];
pdirh  = 3*lfh(4);
pdirp  = [0 1-pdirh 1 pdirh];
pfdh   = 1-(pselh+pfilth+pcnth+pdirh);
pfdp   = [0 pselh+pfilth+pcnth 1 pfdh];
%=======================================================================

%=======================================================================
function resize_fun(fg,varargin)
[pselp, pfiltp, pcntp, pfdp, pdirp] = panelpositions(fg);
if pfdp(4) <= 0
    pfdp(4)=eps;
end
set(findobj(fg,'Tag','msg'), 'Position',pselp);
if pfiltp(4) > 0
    set(findobj(fg,'Tag','pfilt'), 'Position',pfiltp);
end
set(findobj(fg,'Tag','pcnt'), 'Position',pcntp);
set(findobj(fg,'Tag','pfd'), 'Position',pfdp);
set(findobj(fg,'Tag','pdir'), 'Position',pdirp);
return;
%=======================================================================

%=======================================================================
function null(varargin)
%=======================================================================

%=======================================================================
function omsg = msg(varargin)
ob = sib('msg');
omsg = get(ob,'Title');
set(ob,'Title',sprintf(varargin{:}));
drawnow;
return;
%=======================================================================

%=======================================================================
function c = get_current_char
fg = sib(mfilename);
c = get(fg, 'CurrentCharacter');
if ~isempty(c)
    % reset CurrentCharacter
    set(fg, 'CurrentCharacter', char(13));
end
%=======================================================================

%=======================================================================
function obj = sib(tag)
persistent fg;
if isempty(fg) || ~ishandle(fg)
    fg = findobj(0,'Tag',mfilename);
end
obj = findobj(fg,'Tag',tag);
return;
%=======================================================================

%=======================================================================
function varargout = local_isdir(cmd,varargin)
switch lower(cmd)
    case 'list'
        d  = dir(varargin{1});
        dn = {d([d.isdir]).name}';
        [varargout{1:2}] = intersect(varargin{2}, dn);
    case 'filter'
        varargout{1} = varargin{1};
        varargout{2} = 1:numel(varargout{1});
end
%=======================================================================

%=======================================================================
