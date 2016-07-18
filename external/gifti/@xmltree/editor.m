function editor(tree)
%XMLTREE/EDITOR A Graphical User Interface for an XML tree
%  EDITOR(TREE) opens a new figure displaying the xmltree object TREE.
%  H = EDITOR(TREE) also returns the figure handle H.
%
%  This is a beta version of <xmltree/view> successor
%
%  See also XMLTREE
%__________________________________________________________________________
% Copyright (C) 2002-2015  http://www.artefact.tk/

% Guillaume Flandin
% $Id: editor.m 6480 2015-06-13 01:08:30Z guillaume $


%error(nargchk(1,1,nargin));

if ~isempty(getfilename(tree))
    title = getfilename(tree);
elseif ~isempty(inputname(1))
    title = ['Variable ''' inputname(1) ''''];
else
    title = 'Untitled';
end
h = initWindow(title);

setappdata(h,'handles',guihandles(h));
setappdata(h,'save',0);

tree = set(tree,root(tree),'show',1);
setappdata(h,'tree',tree);

doUpdate([],[],h);

set(h,'HandleVisibility','callback');

%==========================================================================

function h = initWindow(title)

wincolor = struct('bg',    [0.8 0.8 0.8], ...
                  'fg',    [1.0 1.0 1.0], ...
                  'title', [0.9 0.9 0.9]);

title = [':: XMLTree Editor :: ' title];

h = figure('Name', title, ...
    'Units', 'Points', ...
    'NumberTitle', 'off', ...
    'Resize', 'on', ...
    'Color', wincolor.bg,...
    'Position', [200 200 440 330], ...
    'MenuBar', 'none', ... 
    'Tag', mfilename);
           
set(h, 'CloseRequestFcn', {@doClose,h});

%- Left box
uicontrol('Style', 'listbox', ...
    'HorizontalAlignment','left', ...
    'Units','Normalized', ...
    'Visible','on',...
    'BackgroundColor', wincolor.fg, ...
    'Max', 1, ...
    'Value', 1 , ...
    'Enable', 'on', ...
    'Position', [0.04 0.12 0.3 0.84], ...
    'Callback', {@doList,h}, ...
    'String', ' ', ...
    'Tag', 'xmllistbox');

%- Right box
uicontrol('Style', 'listbox', ...
    'HorizontalAlignment','left', ...
    'Units','Normalized', ...
    'Visible','on',...
    'BackgroundColor', wincolor.bg, ...
    'Min', 0, ...
    'Max', 2, ...
    'Value', [], ...
    'Enable', 'inactive', ...
    'Position', [0.38 0.50 0.58 0.46], ...
    'Callback', '', ...
    'String', ' ', ...
    'Tag', 'xmllist');

%- The Add button
uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.04 0.03 0.11 0.06], ...
    'String', 'Add', ...
    'Enable','on',...
    'Tag', 'addbutton', ...
    'Callback', {@doAdd,h});
   
%- The Modify button
uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.175 0.03 0.11 0.06], ...
    'String', 'Modify', ...
    'Enable','on',...
    'Tag', 'modifybutton', ...
    'Callback', {@doModify,h});

%- The Copy button
uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.310 0.03 0.11 0.06], ...
    'String', 'Copy', ...
    'Enable','on',...
    'Tag', 'copybutton', ...
    'Callback', {@doCopy,h});

%- The Delete button
uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.445 0.03 0.11 0.06], ...
    'String', 'Delete', ...
    'Enable','on',...
    'Tag', 'delbutton', ...
    'Callback', {@doDelete,h});

%- The Save button
uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.580 0.03 0.11 0.06], ...
    'String', 'Save', ...
    'Tag', 'savebutton', ...
    'Callback', {@doSave,h});

%- The Run button  
%uicontrol('Style', 'pushbutton', ...
%   'Units', 'Normalized', ...
%   'Position', [0.715 0.03 0.11 0.06], ...
%   'String', 'Run', ...
%   'Enable', 'off', ...
%   'Tag', 'runbutton', ...
%   'Callback', {@doRun,h});

%- The Attributes button
uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.715 0.03 0.11 0.06], ...
    'String', 'Attr.', ...
    'Enable', 'on', ...
    'Tag', 'runbutton', ...
    'Callback', {@doAttributes,h});

%- The Close button
uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.850 0.03 0.11 0.06], ...
    'String', 'Close', ...
    'Tag', 'closebutton', ...
    'Callback', {@doClose,h});

%==========================================================================

function doClose(fig,evd,h)
    s = getappdata(h, 'save');
    status = 1;
    if s
        button = questdlg(sprintf('Save changes to the XML tree?'),...
            'XMLTree Editor','Yes','No','Cancel','Yes');
        if strcmp(button,'Yes')
            status = doSave([],[],h);
        elseif strcmp(button,'Cancel')
            status = 0;
        end
    end
    if status
        delete(h);
    end
    
function doAdd(fig,evd,h)
    tree    = getappdata(h, 'tree');
    uidList = getappdata(h, 'uidlist');
    handles = getappdata(h, 'handles');
    pos     = get(handles.xmllistbox,'value');
    uid     = uidList(pos);
    answer = questdlg('Which kind of item to add?', ...
                      'XMLTree Editor :: Add','Node','Leaf','Node');
    switch answer
        case 'Node'
            tree = add(tree,uid,'element','New_Node');
        case 'Leaf'
            tree = add(tree,uid,'element','New_Leaf');
            l = length(tree);
            tree = add(tree,l,'chardata','default');
    end
    tree = set(tree,uid,'show',1);
    setappdata(h, 'tree', tree);
    setappdata(h, 'save', 1);
    doUpdate([],[],h);
    doList([],[],h);

function doModify(fig,evd,h)
    tree    = getappdata(h, 'tree');
    uidList = getappdata(h, 'uidlist');
    handles = getappdata(h, 'handles');
    pos     = get(handles.xmllistbox,'value');
    uid     = uidList(pos);
    contents = children(tree,uid);
    if length(contents) > 0 && ...
       strcmp(get(tree,contents(1),'type'),'chardata')
        str = get(tree,contents(1),'value');
        prompt = {'Name :','New value:'};
        def = {get(tree,uid,'name'),str};
        title = sprintf('Modify %s',get(tree,uid,'name'));
        lineNo = 1;
        answer = inputdlg(prompt,title,lineNo,def);
        if ~isempty(answer)
            tree = set(tree,uid,'name',answer{1});
            str = answer{2};
            tree = set(tree,contents(1),'value',str);
            setappdata(h, 'tree', tree);
            setappdata(h, 'save', 1);
        end
    else
        str = ['Tag ' get(tree,uid,'name')];
        prompt = {'Name :'};
        def = {get(tree,uid,'name'),str};
        title = sprintf('Modify %s tag',get(tree,uid,'name'));
        lineNo = 1;
        answer = inputdlg(prompt,title,lineNo,def);
        if ~isempty(answer)
            tree = set(tree,uid,'name',answer{1});
            str = ['Tag ' get(tree,uid,'name')];
            setappdata(h, 'tree', tree);
            setappdata(h, 'save', 1);
        end
    end
    %- Trying to keep the slider active
    set(handles.xmllist, 'Enable', 'on');
    set(handles.xmllist, 'String', ' ');
    set(handles.xmllist, 'Enable', 'Inactive');
    set(handles.xmllist, 'String', str);
    doUpdate([],[],h);
    doList([],[],h);

function doCopy(fig,evd,h)
    tree    = getappdata(h, 'tree');
    uidList = getappdata(h, 'uidlist');
    handles = getappdata(h, 'handles');
    pos     = get(handles.xmllistbox,'value');
    uid     = uidList(pos);
    if pos ~= 1
        tree = copy(tree,uid);
        setappdata(h, 'tree', tree);
        setappdata(h, 'save', 1);
        doUpdate([],[],h);
        doList([],[],h);
    end
    

function doDelete(fig,evd,h)
    tree    = getappdata(h, 'tree');
    uidList = getappdata(h, 'uidlist');
    handles = getappdata(h, 'handles');
    pos     = get(handles.xmllistbox,'value');
    uid     = uidList(pos);
    if pos > 1
        tree = delete(tree,uid);
        set(handles.xmllistbox,'value',pos-1);
        setappdata(h, 'save', 1);
    end
    setappdata(h, 'tree', tree);
    doUpdate([],[],h);
    doList([],[],h);

function status = doSave(fig,evd,h)
    tree = getappdata(h, 'tree');
    [filename, pathname] = uiputfile({'*.xml' 'XML file (*.xml)'}, ...
                                     'Save XML file as');
    status = ~(isequal(filename,0) | isequal(pathname,0));
    if status
        save(tree,fullfile(pathname, filename));
        set(h,'Name',[':: XMLTree Editor :: ' filename]);
        setappdata(h, 'save', 0);
    end

function doRun(fig,evd,h)
    warndlg('Not implemented','XMLTree :: Run');

function doAttributes(fig,evd,h)
    tree    = getappdata(h, 'tree');
    uidList = getappdata(h, 'uidlist');
    handles = getappdata(h, 'handles');
    pos     = get(handles.xmllistbox,'value');
    uid     = uidList(pos);
    attr = attributes(tree,'get',uid);
    if ~isempty(attr)
        fprintf('This element has %d attributes.\n',length(attr));
        %%%
        %%% Do what you want with 'attr'
        %%% to modify them, use:
        %%%   tree = attributes(tree,'set',uid,n,key,val)
        %%% to add one, use:
        %%%   tree = attributes(tree,'add',uid,key,val)
        %%% to delete one, use:
        %%%   tree = attributes(tree,'del',uid[,n]);
        %%%
        setappdata(h, 'tree', tree);
        setappdata(h, 'save', 1); %- only if attributes have been modified
    end
    
%==========================================================================

function doList(fig,evd,h)
    tree =    getappdata(h, 'tree');
    uidList = getappdata(h, 'uidlist');
    handles = getappdata(h, 'handles');
    uid     = uidList(get(handles.xmllistbox, 'value'));
    
    %- Single mouse click
    if strcmp(get(h,'SelectionType'),'normal')
        contents = children(tree, uid);
        if length(contents) > 0 && ...
           strcmp(get(tree,contents(1),'type'),'chardata')
            str = get(tree,contents(1),'value');
            set(handles.addbutton,'Enable','off');
        elseif length(contents) == 0
            str = '';
            tree = add(tree,uid,'chardata',str);
            setappdata(h, 'tree', tree);
            set(handles.addbutton,'Enable','off');
        else
            str = ['Tag ' get(tree,uid,'name')]; 
            set(handles.addbutton,'Enable','on');
        end
        if get(handles.xmllistbox,'value') == 1
            set(handles.copybutton,'Enable','off');
            set(handles.delbutton,'Enable','off');
        else
            set(handles.copybutton,'Enable','on');
            set(handles.delbutton,'Enable','on');
        end
        %- Trying to keep the slider active
        set(handles.xmllist, 'Enable', 'on');
        set(handles.xmllist, 'String', ' ');
        set(handles.xmllist, 'Enable', 'Inactive');
        set(handles.xmllist, 'String', str);
    %- Double mouse click
    else
        tree = doFlip(tree, uid);
        setappdata(h, 'tree', tree);
        doUpdate([],[],h);
    end

function doUpdate(fig,evd,h)
    tree    = getappdata(h, 'tree');
    handles = getappdata(h, 'handles');
    
    [batchString, uidList] = doUpdateR(tree);
    set(handles.xmllistbox, 'String', batchString);
    setappdata(h, 'uidlist', uidList);
    
%==========================================================================

function [batchString, uidList] = doUpdateR(tree, uid, o)
    if nargin < 2, uid = root(tree); end
    if nargin < 3 || o == 0
        o = 0;
        sep = ' ';
    else
        sep = blanks(4*o);
    end
    if attributes(tree,'length',uid) > 0
        batchString     = {[sep, get(tree, uid, 'name') ' *']};
    else
        batchString     = {[sep, get(tree, uid, 'name')]};
    end
    uidList         = [get(tree,uid,'uid')];
    haselementchild = 0;
    contents        = get(tree, uid, 'contents');
    if isfield(tree, uid, 'show') && get(tree, uid, 'show') == 1
        for i=1:length(contents)
            if strcmp(get(tree,contents(i),'type'),'element')
                [subbatchString, subuidList] = doUpdateR(tree,contents(i),o+1);
                batchString = {batchString{:} subbatchString{:}};
                uidList = [uidList subuidList];
                haselementchild = 1;
            end
        end
        if haselementchild == 1, batchString{1}(length(sep)) = '-'; end
    else
        for i=1:length(contents)
            if strcmp(get(tree,contents(i),'type'),'element')
                haselementchild = 1;
            end
        end
        if haselementchild == 1, batchString{1}(length(sep)) = '+'; end
    end

function tree = doFlip(tree, uid)
    if isfield(tree,uid,'show')
        show = get(tree,uid,'show');
    else
        show = 0;
    end
    tree = set(tree,uid,'show',~show);
