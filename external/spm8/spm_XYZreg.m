function varargout=spm_XYZreg(varargin)
% Registry for GUI XYZ locations, and point list utility functions
%
%                           ----------------
%
% PointList & voxel centre utilities...
%
% FORMAT [xyz,d] = spm_XYZreg('RoundCoords',xyz,M,D)
% FORMAT [xyz,d] = spm_XYZreg('RoundCoords',xyz,V)
% Rounds specified xyz location to nearest voxel centre
% xyz - (Input) 3-vector of X, Y & Z locations, in "real" co-ordinates
% M   - 4x4 transformation matrix relating voxel to "real" co-ordinates
% D   - 3 vector of image X, Y & Z dimensions (DIM)
% V   - 9-vector of image and voxel sizes, and origin [DIM,VOX,ORIGIN]'
%       M derived as [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]]
%       DIM    - D
%       VOX    - Voxel dimensions in units of "real" co-ordinates
%       ORIGIN - Origin of "real" co-ordinates in voxel co-ordinates
% xyz - (Output) co-ordinates of nearest voxel centre in "real" co-ordinates
% d   - Euclidean distance between requested xyz & nearest voxel centre
%
% FORMAT i = spm_XYZreg('FindXYZ',xyz,XYZ)
% finds position of specified voxel in XYZ pointlist
% xyz - 3-vector of co-ordinates
% XYZ - Pointlist: 3xn matrix of co-ordinates
% i   - Column(s) of XYZ equal to xyz
%
% FORMAT [xyz,i,d] = spm_XYZreg('NearestXYZ',xyz,XYZ)
% find nearest voxel in pointlist to specified location
% xyz - (Input) 3-vector of co-ordinates
% XYZ - Pointlist: 3xn matrix of co-ordinates
% xyz - (Output) co-ordinates of nearest voxel in XYZ pointlist
%       (ties are broken in favour of the first location in the pointlist)
% i   - Column of XYZ containing co-ordinates of nearest pointlist location
% d   - Euclidean distance between requested xyz & nearest pointlist location
%
% FORMAT d = spm_XYZreg('Edist',xyz,XYZ)
% Euclidean distances between co-ordinates xyz & points in XYZ pointlist
% xyz - 3-vector of co-ordinates
% XYZ - Pointlist: 3xn matrix of co-ordinates
% d   - n row-vector of Euclidean distances between xyz & points of XYZ
%
%                           ----------------
% Registry functions
%
% FORMAT [hReg,xyz] = spm_XYZreg('InitReg',hReg,M,D,xyz)
% Initialise registry in graphics object
% hReg - Handle of HandleGraphics object to build registry in. Object must
%        be un'Tag'ged and have empty 'UserData'
% M    - 4x4 transformation matrix relating voxel to "real" co-ordinates, used
%        and stored for checking validity of co-ordinates
% D    - 3 vector of image X, Y & Z dimensions (DIM), used
%        and stored for checking validity of co-ordinates
% xyz  - (Input) Initial co-ordinates [Default [0;0;0]]
%        These are rounded to the nearest voxel centre
% hReg - (Output) confirmation of registry handle
% xyz  - (Output) Current registry co-ordinates, after rounding
%
% FORMAT spm_XYZreg('UnInitReg',hReg)
% Clear registry information from graphics object
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object.
%        Object's 'Tag' & 'UserData' are cleared
%
% FORMAT xyz = spm_XYZreg('GetCoords',hReg)
% Get current registry co-ordinates
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% 
% FORMAT [xyz,d] = spm_XYZreg('SetCoords',xyz,hReg,hC,Reg)
% Set co-ordinates in registry & update registered HGobjects/functions
% xyz  - (Input) desired co-ordinates
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
%        If hReg doesn't contain a registry, a warning is printed.
% hC   - Handle of caller object (to prevent circularities) [Default 0]
%        If caller object passes invalid registry handle, then spm_XYZreg
%        attempts to blank the 'hReg' fiend of hC's 'UserData', printing
%        a warning notification.
% Reg  - Alternative nx2 cell array Registry of handles / functions
%        If specified, overrides use of registry held in hReg
%        [Default getfield(get(hReg,'UserData'),'Reg')]
% xyz  - (Output) Desired co-ordinates are rounded to nearest voxel if hC
%        is not specified, or is zero. Otherwise, caller is assummed to
%        have checked verity of desired xyz co-ordinates. Output xyz returns
%        co-ordinates actually set.
% d    - Euclidean distance between desired and set co-ordinates.
%
% FORMAT nReg = spm_XYZreg('XReg',hReg,{h,Fcn}pairs)
% Cross registration object/function pairs with the registry, push xyz co-ords
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% h    - Handle of HandleGraphics object to be registered
%        The 'UserData' of h must be a structure with an 'Reg' field, which
%        is set to hReg, the handle of the registry (back registration)
% Fcn  - Handling function for HandleGraphics object h
%        This function *must* accept XYZ updates via the call:
%                feval(Fcn,'SetCoords',xyz,h,hReg)
%        and should *not* call back the registry with the update!
%        {h,Fcn} are appended to the registry (forward registration)
% nReg - New registry cell array: Handles are checked for validity before
%        entry. Invalid handles are omitted, generating a warning.
%
% FORMAT nReg = spm_XYZreg('Add2Reg',hReg,{h,Fcn}pairs)
% Add object/function pairs for XYZ updates to registry (forward registration)
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% h    - Handle of HandleGraphics object to be registered
% Fcn  - Handling function for HandleGraphics object h
%        This function *must* accept XYZ updates via the call:
%                feval(Fcn,'SetCoords',xyz,h,hReg)
%        and should *not* call back the registry with the update!
%        {h,Fcn} are appended to the registry (forward registration)
% nReg - New registry cell array: Handles are checked for validity before
%        entry. Invalid handles are omitted, generating a warning.
%
% FORMAT spm_XYZreg('SetReg',h,hReg)
% Set registry field of object's UserData (back registration)
% h    - Handle of HandleGraphics object to be registered
%        The 'UserData' of h must be a structure with an 'Reg' field, which
%        is set to hReg, the handle of the registry (back registration)
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
%
% FORMAT nReg = spm_XYZreg('unXReg',hReg,hD1,hD2,hD3,...)
% Un-cross registration of HandleGraphics object hD
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% hD?  - Handles of HandleGraphics object to be unregistered
%        The 'UserData' of hD must be a structure with a 'Reg' field, which
%        is set to empty (back un-registration)
% nReg - New registry cell array: Registry entries with handle entry hD are 
%        removed from the registry (forward un-registration)
%        Handles not in the registry generate a warning
%
% FORMAT nReg = spm_XYZreg('Del2Reg',hReg,hD)
% Delete HandleGraphics object hD from registry (forward un-registration)
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% hD?  - Handles of HandleGraphics object to be unregistered
% nReg - New registry cell array: Registry entries with handle entry hD are 
%        removed from the registry. Handles not in registry generate a warning
%
% FORMAT spm_XYZreg('UnSetReg',h)
% Unset registry field of object's UserData (back un-registration)
% h - Handle of HandleGraphics object to be unregistered
%     The 'UserData' of hD must be a structure with a 'Reg' field, which
%     is set to empty (back un-registration)
%
% FORMAT spm_XYZreg('CleanReg',hReg)
% Clean invalid handles from registry
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
%
% FORMAT Reg = spm_XYZreg('VReg',Reg,Warn)
% Prune invalid handles from Registry cell array
% Reg  - (Input) nx2 cell array of {handle,function} pairs
% Warn - If specified, print warning if find invalid handles
% Reg  - (Output) mx2 cell array of valid {handle,function} pairs
%
% FORMAT hReg = spm_XYZreg('FindReg',h)
% Find/check registry object
% h    - handle of Registry, or figure containing Registry (default gcf)
%        If ischar(h), then uses spm_figure('FindWin',h) to locate named figures
% hReg - handle of confirmed registry object
%        Errors if h is not a registry or a figure containing a unique registry
%        Registry object is identified by 'hReg' 'Tag'
%_______________________________________________________________________
%
% spm_XYZreg provides a framework for modular inter-GUI communication of
% XYZ co-orginates, and various utility functions for pointlist handling
% and rounding in voxel co-ordinates.
%
%-----------------------------------------------------------------------
%                                                           THE REGISTRY
%
% The concept of the registry is of a central entity which "knows"
% about other GUI objects holding XYZ co-ordinates, and keeps them all
% in sync. Changes to the registry's XYZ co-ordinates are passed on to
% registered functions by the registry (forward registration).
% Individual objects which can change the XYZ co-ordinates should
% therefore update the registry with the new co-ordinates (back
% registration), so that the registry can tell all registered objects
% about the new location, and a framework is provided for this.
%
% The registry is held as the 'UserData of a HandleGraphics object,
% whose handle therefore identifies the registry. The registry object
% is 'Tag'ged 'hReg' for identification (though this 'Tag' is not used
% for locating the registry, so multiple registry incarnations are
% possible). The registry object's 'UserData' is a structure containing
% the current XYZ co-ordinates, the voxel-to-co-ordinates matrix M, the
% image dimensions D, and the Registry itself. The registry is a nx2
% cell array containing n handle/function pairs.
%
% The model is that all GUI objects requiring linking to a common XYZ
% location via the registry each be identified by a HandleGraphics
% handle. This handle can be the handle of the particular instantiation
% of the GUI control itself (as is the case with the MIP-GUI of
% spm_mip_ui where the axis handle is used to identify the MIP to use);
% the handle of another HandleGraphics object associated with the GUI
% control (as is the case with the XYZ editable widgets of
% spm_results_ui where the handle of the bounding frame uicontrol is
% used); or may be 0, the handle of the root object, which allows non
% GUI functions (such as a function that just prints information) to be
% added to the registry. The registry itself thus conforms to this
% model. Each object has an associated "handling function" (so this
% function is the registry's handling function). The registry itself
% consists of object-handle/handling-function pairs.
%
% If an object and it's handling function are entered in the registry,
% then the object is said to be "forward registered", because the
% registry will now forward all location updates to that object, via
% it's handling function. The assummed syntax is:
% feval(Fcn,'SetCoords',xyz,h,hReg), where Fcn is the handling function
% for the GUI control identified by handle h, xyz are the new
% co-ordinates, and hReg is the handle of the registry.
%
% An optional extension is "back registration", whereby the GUI
% controls inform the registry of the new location when they are
% updated. All that's required is that the objects call the registry's
% 'SetCoords' function: spm_XYZreg('SetCoords',xyz,hReg,hC), where hReg
% is the registry object's handle, and hC is the handle associated with
% the calling GUI control. The specification of the caller GUI control
% allows the registry to avoid circularities: If the object is "forward
% registered" for updates, then the registry function doesn't try to
% update the object which just updated the registry! (Similarly, the
% handle of the registry object, hReg, is passed to the handling
% function during forward XYZ updating, so that the handling function's
% 'SetCoords' facility can be constructed to accept XYZ updates from
% various sources, and only inform the registry if not called by the
% registry, and hence avoid circularities.)
%
% A framework is provided for "back" registration. Really all that is
% required is that the GUI controls know of the registry object (via
% it's handle hReg), and call the registry's 'SetCoords' facility when
% necessary. This can be done in many ways, but a simple structure is
% provided, mirroring that of the registry's operation. This framework
% assummes that the GUI controls identification object's 'UserData' is
% a structure with a field named 'hReg', which stores the handle of the
% registry (if back registered), or is empty (if not back registered,
% i.e. standalone). spm_XYZreg provides utility functions for
% setting/unsetting this field, and for "cross registering" - that is
% both forward and back registration in one command. Cross registering
% involves adding the handle/function pair to the registry, and setting
% the registry handle in the GUI control object's 'UserData' 'hReg'
% field. It's up to the handling function to read the registry handle
% from it's objects 'UserData' and act accordingly. A simple example of
% such a function is provided in spm_XYZreg_Ex2.m, illustrated below.
%
% SubFunctions are provided for getting and setting the current
% co-ordinates; adding and deleting handle/function pairs from the
% registry (forward registration and un-registration), setting and
% removing registry handle information from the 'hReg' field of the
% 'UserData' of a HG object (backward registration & un-registration);
% cross registration and unregistration (including pushing of current
% co-ordinates); setting and getting the current XYZ location. See the
% FORMAT statements and the example below...
%
%                           ----------------
% Example
% %-Create a window:
% F = figure;
% %-Create an object to hold the registry
% hReg = uicontrol(F,'Style','Text','String','hReg',...
%   'Position',[100 200 100 025],...
%   'FontName','Times','FontSize',14,'FontWeight','Bold',...
%   'HorizontalAlignment','Center');
% %-Setup M & D
% V = [65;87;26;02;02;04;33;53;08];
% M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
% D = V(1:3);
% %-Initialise a registry in this object, with initial co-ordinates [0;0;0]
% spm_XYZreg('InitReg',hReg,M,D,[0;0;0])
% % (ans returns [0;0;0] confirming current co-ordinates
% %-Set co-ordinates to [10;10;10]
% spm_XYZreg('SetCoords',[10,10,10],hReg)
% % (warns of co-ordinate rounding to [10,10,12], & returns ans as [10;10;12])
%
% %-Forward register a command window xyz reporting function: spm_XYZreg_Ex1.m
% spm_XYZreg('Add2Reg',hReg,0,'spm_XYZreg_Ex1')
% % (ans returns new registry, containing just this handle/function pair
% %-Set co-ordinates to [0;10;12]
% [xyz,d] = spm_XYZreg('SetCoords',[0,10,12],hReg);
% % (spm_XYZreg_Ex1 called, and prints co-ordinates and handles)
% %-Have a peek at the registry information
% RD = get(hReg,'UserData')
% RD.xyz    %-The current point according to the registry
% RD.Reg    %-The nx2 cell array of handle/function pairs
%
% %-Create an example GUI XYZ control, using spm_XYZreg_Ex2.m
% hB = spm_XYZreg_Ex2('Create',M,D,xyz);
% % (A figure window with a button appears, whose label shows the current xyz
% %-Press the button, and enter new co-ordinates [0;0;0] in the Cmd window...
% % (...the button's internal notion of the current location is changed, but
% % (the registry isn't informed:
% spm_XYZreg('GetCoords',hReg)
% (...returns [0;10;12])
% %-"Back" register the button
% spm_XYZreg('SetReg',hB,hReg)
% %-Check the back registration
% if ( hReg == getfield(get(hB,'UserData'),'hReg') ), disp('yes!'), end
% %-Now press the button, and enter [0;0;0] again...
% % (...this time the registry is told, and the registry tells spm_XYZreg_Ex1,
% % (which prints out the new co-ordinates!
% %-Forward register the button to receive updates from the registry
% nReg = spm_XYZreg('Add2Reg',hReg,hB,'spm_XYZreg_Ex2')
% % (The new registry is returned as nReg, showing two entries
% %-Set new registry co-ordinates to [10;10;12]
% [xyz,d] = spm_XYZreg('SetCoords',[10;10;12],hReg);
% % (...the button updates too!
%
% %-Illustration of robustness: Delete the button & use the registry
% delete(hB)
% [xyz,d] = spm_XYZreg('SetCoords',[10;10;12],hReg);
% % (...the invalid handle hB in the registry is ignored)
% %-Peek at the registry
% getfield(get(hReg,'UserData'),'Reg')
% %-Delete hB from the registry by "cleaning"
% spm_XYZreg('CleanReg',hReg)
% % (...it's gone
%
% %-Make a new button and cross register
% hB = spm_XYZreg_Ex2('Create',M,D)
% % (button created with default co-ordinates of [0;0;0]
% nReg = spm_XYZreg('XReg',hReg,hB,'spm_XYZreg_Ex2')
% % (Note that the registry pushes the current co-ordinates to the button
% %-Use the button & spm_XYZreg('SetCoords'... at will!
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes, Chloe Hutton
% $Id$



%=======================================================================
switch lower(varargin{1}), case 'roundcoords'
%=======================================================================
% [xyz,d] = spm_XYZreg('RoundCoords',xyz,M,D)
% [xyz,d] = spm_XYZreg('RoundCoords',xyz,V)
if nargin<3, error('Insufficient arguments'), end
if nargin<4
    V = varargin{3};
    M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
    D = V(1:3);
else
    M = varargin{3};
    D = varargin{4};
end
    
%-Round xyz to coordinates of actual voxel centre
%-Do rounding in voxel coordinates & ensure within image size
%-Watch out for infinities!
%-----------------------------------------------------------------------
xyz  = [varargin{2}(:); 1];
xyz(isinf(xyz)) = 1e10*sign(xyz(isinf(xyz)));
rcp  = round(inv(M)*xyz);
rcp  = max([min([rcp';[D',1]]);[1,1,1,1]])';
rxyz = M*rcp;

%-Work out Euclidean distance between points xyz & rounded xyz
d = sqrt(sum((xyz-rxyz).^2));

varargout = {rxyz(1:3),d};



%=======================================================================
case 'findxyz'
%=======================================================================
% i = spm_XYZreg('FindXYZ',xyz,XYZ)
if nargin<3, error('Insufficient arguments'), end
XYZ = varargin{3};
xyz = varargin{2};
    
%-Find XYZ = xyz
%-----------------------------------------------------------------------
i = find(all([XYZ(1,:)==xyz(1);XYZ(2,:)==xyz(2);XYZ(3,:)==xyz(3)],1));

varargout = {i};



%=======================================================================
case 'nearestxyz'
%=======================================================================
% [xyz,i,d] = spm_XYZreg('NearestXYZ',xyz,XYZ)
if nargin<3, error('Insufficient arguments'), end
    
%-Find in XYZ nearest point to coordinates xyz (Euclidean distance) 
%-----------------------------------------------------------------------
[d,i] = min(spm_XYZreg('Edist',varargin{2},varargin{3}));
varargout = {varargin{3}(:,i),i,d};



%=======================================================================
case 'edist'
%=======================================================================
% d = spm_XYZreg('Edist',xyz,XYZ)
if nargin<3, error('Insufficient arguments'), end
    
%-Calculate (Euclidean) distances from pointlist co-ords to xyz
%-----------------------------------------------------------------------
varargout = {sqrt(sum([ (varargin{3}(1,:) - varargin{2}(1));...
            (varargin{3}(2,:) - varargin{2}(2));...
            (varargin{3}(3,:) - varargin{2}(3)) ].^2))};



%=======================================================================
case 'initreg'      % Initialise registry in handle h
%=======================================================================
% [hReg,xyz] = spm_XYZreg('InitReg',hReg,M,D,xyz)
if nargin<5, xyz=[0;0;0]; else, xyz=varargin{5}; end
if nargin<4, error('Insufficient arguments'), end
D    = varargin{4};
M    = varargin{3};
hReg = varargin{2};

%-Check availability of hReg object for building a registry in
%-----------------------------------------------------------------------
if ~isempty(get(hReg,'UserData')), error('Object already has UserData...'), end
if ~isempty(get(hReg,'Tag')), error('Object already ''Tag''ed...'), end

%-Check co-ordinates are in range
%-----------------------------------------------------------------------
[xyz,d] = spm_XYZreg('RoundCoords',xyz,M,D);
if d>0 & nargout<2, warning(sprintf('%s: Co-ords rounded to neatest voxel center: Discrepancy %.2f',mfilename,d)), end

%-Set up registry
%-----------------------------------------------------------------------
RD = struct('xyz',xyz,'M',M,'D',D,'Reg',[]);
RD.Reg = {};
set(hReg,'Tag','hReg','UserData',RD)

%-Return current co-ordinates
%-----------------------------------------------------------------------
varargout = {hReg,xyz};



%=======================================================================
case 'uninitreg'    % UnInitialise registry in handle hReg
%=======================================================================
% spm_XYZreg('UnInitReg',hReg)
hReg = varargin{2};
if ~strcmp(get(hReg,'Tag'),'hReg'), warning('Not an XYZ registry'), return, end
set(hReg,'Tag','','UserData',[])



%=======================================================================
case 'getcoords'    % Get current co-ordinates
%=======================================================================
% xyz = spm_XYZreg('GetCoords',hReg)
if nargin<2, hReg=spm_XYZreg('FindReg'); else, hReg=varargin{2}; end
if ~ishandle(hReg), error('Invalid object handle'), end
if ~strcmp(get(hReg,'Tag'),'hReg'), error('Not a registry'), end
varargout = {getfield(get(hReg,'UserData'),'xyz')};



%=======================================================================
case 'setcoords'    % Set co-ordinates & update registered functions
%=======================================================================
% [xyz,d] = spm_XYZreg('SetCoords',xyz,hReg,hC,Reg)
% d returned empty if didn't check, warning printed if d not asked for & round
% Don't check if callerhandle specified (speed)
% If Registry cell array Reg is specified, then only these handles are updated
hC=0; mfn=''; if nargin>=4
    if ~ischar(varargin{4}), hC=varargin{4}; else mfn=varargin{4}; end
end
hReg = varargin{3};

%-Check validity of hReg registry handle
%-----------------------------------------------------------------------
%-Return if hReg empty, in case calling objects functions don't check isempty
if isempty(hReg), return, end
%-Check validity of hReg registry handle, correct calling objects if necc.
if ~ishandle(hReg)
    str = sprintf('%s: Invalid registry handle (%.4f)',mfilename,hReg);
    if hC>0
        %-Remove hReg from caller
        spm_XYZreg('SetReg',hC,[])
        str = [str,sprintf('\n\t\t\t...removed from caller (%.4f)',hC)];
    end
    warning(str)
    return
end
xyz  = varargin{2};

RD      = get(hReg,'UserData');

%-Check validity of coords only when called without a caller handle
%-----------------------------------------------------------------------
if hC<=0
    [xyz,d] = spm_XYZreg('RoundCoords',xyz,RD.M,RD.D);
    if d>0 & nargout<2, warning(sprintf(...
        '%s: Co-ords rounded to neatest voxel center: Discrepancy %.2f',...
        mfilename,d)), end
else
    d = 0;
end

%-Sort out valid handles, eliminate caller handle, update co-ords with
% registered handles via their functions
%-----------------------------------------------------------------------
if nargin<5
    RD.Reg = spm_XYZreg('VReg',RD.Reg);
    Reg    = RD.Reg;
else
    Reg = spm_XYZreg('VReg',varargin{5});
end
if hC>0 & length(Reg), Reg(find([Reg{:,1}]==varargin{4}),:) = []; end
for i = 1:size(Reg,1)
    feval(Reg{i,2},'SetCoords',xyz,Reg{i,1},hReg);
end

%-Update registry (if using hReg) with location & cleaned Reg cellarray
%-----------------------------------------------------------------------
if nargin<5
    RD.xyz  = xyz;
    set(hReg,'UserData',RD)
end

varargout = {xyz,d};

if ~strcmp(mfn,'spm_graph')
    sHdl=findobj(0,'Tag','SPMGraphSatelliteFig');
    axHdl=findobj(sHdl,'Type','axes','Tag','SPMGraphSatelliteAxes');
    %tag for true axis, as legend is of type axis, too
    for j=1:length(axHdl)
        autoinp=get(axHdl(j),'UserData');
        if ~isempty(autoinp), spm_graph([],[],hReg,axHdl(j)); end
    end
end


%=======================================================================
case 'xreg'     % Cross register object handles & functions
%=======================================================================
% nReg = spm_XYZreg('XReg',hReg,{h,Fcn}pairs)
if nargin<4, error('Insufficient arguments'), end
hReg = varargin{2};

%-Quick check of registry handle
%-----------------------------------------------------------------------
if isempty(hReg),   warning('Empty registry handle'), return, end
if ~ishandle(hReg), warning('Invalid registry handle'), return, end

%-Condition nReg cell array & check validity of handles to be registered
%-----------------------------------------------------------------------
nReg = varargin(3:end);
if mod(length(nReg),2), error('Registry items must be in pairs'), end
if length(nReg)>2, nReg = reshape(nReg,length(nReg)/2,2)'; end
nReg = spm_XYZreg('VReg',nReg,'Warn');

%-Set hReg registry link for registry candidates (Back registration)
%-----------------------------------------------------------------------
for i = 1:size(nReg,1)
    spm_XYZreg('SetReg',nReg{i,1},hReg);
end

%-Append registry candidates to existing registry & write back to hReg
%-----------------------------------------------------------------------
RD     = get(hReg,'UserData');
Reg    = RD.Reg;
Reg    = cat(1,Reg,nReg);
RD.Reg = Reg;
set(hReg,'UserData',RD)

%-Synch co-ordinates of newly registered objects
%-----------------------------------------------------------------------
spm_XYZreg('SetCoords',RD.xyz,hReg,hReg,nReg);

varargout = {Reg};



%=======================================================================
case 'add2reg'      % Add handle(s) & function(s) to registry
%=======================================================================
% nReg = spm_XYZreg('Add2Reg',hReg,{h,Fcn}pairs)
if nargin<4, error('Insufficient arguments'), end
hReg = varargin{2};

%-Quick check of registry handle
%-----------------------------------------------------------------------
if isempty(hReg),   warning('Empty registry handle'), return, end
if ~ishandle(hReg), warning('Invalid registry handle'), return, end

%-Condition nReg cell array & check validity of handles to be registered
%-----------------------------------------------------------------------
nReg = varargin(3:end);
if mod(length(nReg),2), error('Registry items must be in pairs'), end
if length(nReg)>2, nReg = reshape(nReg,length(nReg)/2,2)'; end
nReg = spm_XYZreg('VReg',nReg,'Warn');

%-Append to existing registry & put back in registry object
%-----------------------------------------------------------------------
RD     = get(hReg,'UserData');
Reg    = RD.Reg;
Reg    = cat(1,Reg,nReg);
RD.Reg = Reg;
set(hReg,'UserData',RD)

varargout = {Reg};



%=======================================================================
case 'setreg'           %-Set registry field of object's UserData
%=======================================================================
% spm_XYZreg('SetReg',h,hReg)
if nargin<3, error('Insufficient arguments'), end
h    = varargin{2};
hReg = varargin{3};
if ( ~ishandle(h) | h==0 ), return, end
UD = get(h,'UserData');
if ~isstruct(UD) | ~any(strcmp(fieldnames(UD),'hReg'))
    error('No UserData structure with hReg field for this object')
end
UD.hReg = hReg;
set(h,'UserData',UD)



%=======================================================================
case 'unxreg'       % Un-cross register object handles & functions
%=======================================================================
% nReg = spm_XYZreg('unXReg',hReg,hD1,hD2,hD3,...)
if nargin<3, error('Insufficient arguments'), end
hD   = [varargin{3:end}];
hReg = varargin{2};

%-Get Registry information
%-----------------------------------------------------------------------
RD         = get(hReg,'UserData');
Reg        = RD.Reg;

%-Find registry entires to delete
%-----------------------------------------------------------------------
[null,i,e] = intersect([Reg{:,1}],hD);
hD(e)      = [];
dReg       = spm_XYZreg('VReg',Reg(i,:));
Reg(i,:)   = [];
if length(hD), warning('Not all handles were in registry'), end

%-Write back new registry
%-----------------------------------------------------------------------
RD.Reg = Reg;
set(hReg,'UserData',RD)

%-UnSet hReg registry link for hD's still existing (Back un-registration)
%-----------------------------------------------------------------------
for i = 1:size(dReg,1)
    spm_XYZreg('SetReg',dReg{i,1},[]);
end

varargout = {Reg};



%=======================================================================
case 'del2reg'      % Delete handle(s) & function(s) from registry
%=======================================================================
% nReg = spm_XYZreg('Del2Reg',hReg,hD)
if nargin<3, error('Insufficient arguments'), end
hD   = [varargin{3:end}];
hReg = varargin{2};

%-Get Registry information
%-----------------------------------------------------------------------
RD         = get(hReg,'UserData');
Reg        = RD.Reg;

%-Find registry entires to delete
%-----------------------------------------------------------------------
[null,i,e] = intersect([Reg{:,1}],hD);
Reg(i,:)   = [];
hD(e)      = [];
if length(hD), warning('Not all handles were in registry'), end

%-Write back new registry
%-----------------------------------------------------------------------
RD.Reg = Reg;
set(hReg,'UserData',RD)

varargout = {Reg};



%=======================================================================
case 'unsetreg'         %-Unset registry field of object's UserData
%=======================================================================
% spm_XYZreg('UnSetReg',h)
if nargin<2, error('Insufficient arguments'), end
spm_XYZreg('SetReg',varargin{2},[])



%=======================================================================
case 'cleanreg'     % Clean invalid handles from registry
%=======================================================================
% spm_XYZreg('CleanReg',hReg)
%if ~strcmp(get(hReg,'Tag'),'hReg'), error('Not a registry'), end
hReg = varargin{2};
RD = get(hReg,'UserData');
RD.Reg = spm_XYZreg('VReg',RD.Reg,'Warn');
set(hReg,'UserData',RD)


%=======================================================================
case 'vreg'     % Prune invalid handles from registry cell array
%=======================================================================
% Reg = spm_XYZreg('VReg',Reg,Warn)
if nargin<3, Warn=0; else, Warn=1; end
Reg = varargin{2};
if isempty(Reg), varargout={Reg}; return, end
i = find(~ishandle([Reg{:,1}]));
%-***check existance of handling functions : exist('','file')?
if Warn & length(i), warning([...
    sprintf('%s: Disregarding invalid registry handles:\n\t',...
        mfilename),sprintf('%.4f',Reg{i,1})]), end
Reg(i,:)  = [];
varargout = {Reg};



%=======================================================================
case 'findreg'          % Find/check registry object
%=======================================================================
% hReg = spm_XYZreg('FindReg',h)
if nargin<2, h=get(0,'CurrentFigure'); else, h=varargin{2}; end
if ischar(h), h=spm_figure('FindWin',h); end
if ~ishandle(h), error('invalid handle'), end
if ~strcmp(get(h,'Tag'),'hReg'), h=findobj(h,'Tag','hReg'); end
if isempty(h), error('Registry object not found'), end
if length(h)>1, error('Multiple registry objects found'), end
varargout = {h};



%=======================================================================
otherwise
%=======================================================================
warning('Unknown action string')

%=======================================================================
end
