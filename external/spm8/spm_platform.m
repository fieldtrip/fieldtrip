function varargout=spm_platform(varargin)
% Platform specific configuration parameters for SPM
%
% FORMAT ans = spm_platform(arg)
% arg  - optional string argument, can be
%        - 'bigend'  - return whether this architecture is bigendian
%                      - 0   - is little endian
%                      - 1   - is big endian
%        - 'filesys' - type of filesystem
%                      - 'unx' - UNIX
%                      - 'win' - DOS
%        - 'sepchar' - returns directory separator
%        - 'user'    - returns username
%        - 'host'    - returns system's host name
%        - 'tempdir' - returns name of temp directory
%        - 'drives'  - returns string containing valid drive letters
%
% FORMAT PlatFontNames = spm_platform('fonts')
% Returns structure with fields named after the generic (UNIX) fonts, the
% field containing the name of the platform specific font.
%
% FORMAT PlatFontName = spm_platform('font',GenFontName)
% Maps generic (UNIX) FontNames to platform specific FontNames
%
% FORMAT PLATFORM = spm_platform('init',comp)
% Initialises platform specific parameters in persistent PLATFORM
% (External gateway to init_platform(comp) subfunction)
% comp         - computer to use [defaults to MATLAB's `computer`]
% PLATFORM - copy of persistent PLATFORM
%
% FORMAT spm_platform
% Initialises platform specific parameters in persistent PLATFORM
% (External gateway to init_platform(computer) subfunction)
%
%                              ----------------
% SUBFUNCTIONS:
%
% FORMAT init_platform(comp)
% Initialise platform specific parameters in persistent PLATFORM
% comp         - computer to use [defaults to MATLAB's `computer`]
%
%--------------------------------------------------------------------------
%
% Since calls to spm_platform will be made frequently, most platform
% specific parameters are stored as a structure in the persistent variable
% PLATFORM. Subsequent calls use the information from this persistent
% variable, if it exists.
%
% Platform specific definitions are contained in the data structures at
% the beginning of the init_platform subfunction at the end of this
% file.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Matthew Brett
% $Id$


%-Initialise
%--------------------------------------------------------------------------
persistent PLATFORM
if isempty(PLATFORM), PLATFORM = init_platform; end

if nargin==0, return, end


switch lower(varargin{1}), case 'init'                     %-(re)initialise
%==========================================================================
init_platform(varargin{2:end});
varargout = {PLATFORM};
   
case 'bigend'                         %-Return endian for this architecture
%==========================================================================
varargout = {PLATFORM.bigend};

case 'filesys'                                         %-Return file system
%==========================================================================
varargout = {PLATFORM.filesys};

case 'sepchar'                            %-Return file separator character
%==========================================================================
warning('use filesep instead (supported by MathWorks)')
varargout = {PLATFORM.sepchar};

case 'rootlen'           %-Return length in chars of root directory name 
%=======================================================================
varargout = {PLATFORM.rootlen};

case 'user'                                            %-Return user string
%==========================================================================
varargout = {PLATFORM.user};

case 'host'                                               %-Return hostname
%==========================================================================
varargout = {PLATFORM.host};

case 'drives'                                               %-Return drives
%==========================================================================
varargout = {PLATFORM.drives};

case 'tempdir'                                 %-Return temporary directory
%==========================================================================
twd = getenv('SPMTMP');
if isempty(twd)
    twd = tempdir;
end 
varargout = {twd};

case {'font','fonts'}       %-Map default font names to platform font names
%==========================================================================
if nargin<2, varargout={PLATFORM.font}; return, end
switch lower(varargin{2})
    case 'times'
        varargout = {PLATFORM.font.times};
    case 'courier'
        varargout = {PLATFORM.font.courier};
    case 'helvetica'
        varargout = {PLATFORM.font.helvetica};
    case 'symbol'
        varargout = {PLATFORM.font.symbol};
    otherwise
        warning(['Unknown font ',varargin{2},', using default'])
        varargout = {PLATFORM.font.helvetica};
end

    otherwise                                       %-Unknown Action string
%==========================================================================
error('Unknown Action string')

%==========================================================================
end



%==========================================================================
%- S U B - F U N C T I O N S
%==========================================================================


function PLATFORM = init_platform(comp)     %-Initialise platform variables
%==========================================================================
if nargin<1, comp=computer; end

%-Platform definitions
%--------------------------------------------------------------------------
PDefs = {'PCWIN',     'win',   0;...
         'PCWIN64',   'win',   0;...
         'MAC',       'unx',   1;...
         'MACI',      'unx',   0;...
         'MACI64',    'unx',   0;...
         'SOL2',      'unx',   1;...
         'SOL64',     'unx',   1;...
         'GLNX86',    'unx',   0;...
         'GLNXA64',   'unx',   0};

PDefs = cell2struct(PDefs,{'computer','filesys','endian'},2);


%-Which computer?
%--------------------------------------------------------------------------
[issup, ci] = ismember(comp,{PDefs.computer});
if ~issup
    error([comp ' not supported architecture for ' spm('Ver')]);
end


%-Set byte ordering
%--------------------------------------------------------------------------
PLATFORM.bigend = PDefs(ci).endian;


%-Set filesystem type
%--------------------------------------------------------------------------
PLATFORM.filesys = PDefs(ci).filesys;


%-File separator character
%--------------------------------------------------------------------------
PLATFORM.sepchar = filesep;


%-Username
%--------------------------------------------------------------------------
switch PLATFORM.filesys
    case 'unx'
        PLATFORM.user = getenv('USER');
    case 'win'
        PLATFORM.user = getenv('USERNAME');
    otherwise
        error(['Don''t know filesystem ',PLATFORM.filesys])
end
if isempty(PLATFORM.user), PLATFORM.user = 'anonymous'; end


%-Hostname
%--------------------------------------------------------------------------
[sts, Host]  = system('hostname');
if sts
    if strcmp(PLATFORM.filesys,'win')
        Host = getenv('COMPUTERNAME');
    else
        Host = getenv('HOSTNAME'); 
    end
    Host = regexp(Host,'(.*?)\.','tokens','once');
else
    Host = Host(1:end-1);
end
PLATFORM.host = Host;


%-Drives
%--------------------------------------------------------------------------
PLATFORM.drives = '';
if strcmp(PLATFORM.filesys,'win')
    driveLett = cellstr(char(('C':'Z')'));
    for i=1:numel(driveLett)
        if exist([driveLett{i} ':\'],'dir')
            PLATFORM.drives = [PLATFORM.drives driveLett{i}];
        end
    end
end


%-Fonts
%--------------------------------------------------------------------------
switch comp
    case {'MAC','MACI','MACI64'}
        PLATFORM.font.helvetica = 'TrebuchetMS';
        PLATFORM.font.times     = 'Times';
        PLATFORM.font.courier   = 'Courier';
        PLATFORM.font.symbol    = 'Symbol';
    case {'SOL2','SOL64','GLNX86','GLNXA64'}
        PLATFORM.font.helvetica = 'Helvetica';
        PLATFORM.font.times     = 'Times';
        PLATFORM.font.courier   = 'Courier';
        PLATFORM.font.symbol    = 'Symbol';
    case {'PCWIN','PCWIN64'}
        PLATFORM.font.helvetica = 'Arial Narrow';
        PLATFORM.font.times     = 'Times New Roman';
        PLATFORM.font.courier   = 'Courier New';
        PLATFORM.font.symbol    = 'Symbol';
end
