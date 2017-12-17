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
%        - 'desktop' - returns whether or not the Desktop is in use
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
% the beginning of the init_platform subfunction at the end of this file.
%__________________________________________________________________________
% Copyright (C) 1999-2014 Wellcome Trust Centre for Neuroimaging

% Matthew Brett
% $Id: spm_platform.m 6245 2014-10-15 11:22:15Z guillaume $


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
warning('Use FILESEP instead.')
varargout = {PLATFORM.sepchar};

case 'user'                                            %-Return user string
%==========================================================================
varargout = {PLATFORM.user};

case 'host'                                               %-Return hostname
%==========================================================================
varargout = {PLATFORM.host};

case 'drives'                                               %-Return drives
%==========================================================================
warning('Use spm_select(''ListDrives'') instead.');
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

case 'desktop'                                       %-Return desktop usage
%==========================================================================
varargout = {PLATFORM.desktop};

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
if nargin<1
    if ~strcmpi(spm_check_version,'octave')
        comp = computer;
    else
        if isunix
            switch uname.machine
                case {'x86_64'}
                    comp = 'GLNXA64';
                case {'i586','i686'}
                    comp = 'GLNX86';
                case {'armv6l'}
                    comp = 'ARM';
                otherwise
                    error('%s is not supported.',comp);
            end
        elseif ispc
            comp = 'PCWIN64';
        elseif ismac
            comp = 'MACI64';
        end
    end
end

%-Platform definitions
%--------------------------------------------------------------------------
PDefs = {'PCWIN',     'win',   0;...
         'PCWIN64',   'win',   0;...
         'MAC',       'unx',   1;...
         'MACI',      'unx',   0;...
         'MACI64',    'unx',   0;...
         'GLNX86',    'unx',   0;...
         'GLNXA64',   'unx',   0;...
         'ARM',       'unx',   0};

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
end
if isempty(PLATFORM.user), PLATFORM.user = 'anonymous'; end


%-Hostname
%--------------------------------------------------------------------------
switch PLATFORM.filesys
    case 'unx'
        [sts, PLATFORM.host] = system('hostname');
        if sts
            PLATFORM.host = getenv('HOSTNAME');
        else
            PLATFORM.host = PLATFORM.host(1:end-1);
        end
    case 'win'
        PLATFORM.host = getenv('COMPUTERNAME');
end
PLATFORM.host = strtok(PLATFORM.host,'.');


%-Drives
%--------------------------------------------------------------------------
PLATFORM.drives = '';
if strcmp(PLATFORM.filesys,'win')
    driveLett = spm_select('ListDrives');
    PLATFORM.drives = strrep(strcat(driveLett{:}),':','');
end


%-Fonts
%--------------------------------------------------------------------------
switch comp
    case {'MAC','MACI','MACI64'}
        PLATFORM.font.helvetica = 'TrebuchetMS';
        PLATFORM.font.times     = 'Times';
        PLATFORM.font.courier   = 'Courier';
        PLATFORM.font.symbol    = 'Symbol';
    case {'GLNX86','GLNXA64'}
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


%-Desktop
%--------------------------------------------------------------------------
try
    PLATFORM.desktop = desktop('-inuse');
catch
    PLATFORM.desktop = false;
end
