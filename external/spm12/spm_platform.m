function varargout = spm_platform(varargin)
% Platform specific configuration parameters
%
% FORMAT ans = spm_platform(param)
% param - optional string argument, can be
%         - 'bigend'  - return whether this architecture is big endian
%                       - false  - is little endian
%                       - true   - is big endian
%         - 'mexext'  - return MEX filename extension
%         - 'soext'   - return shared library filename extension
%         - 'user'    - return username
%         - 'host'    - return system's host name
%         - 'tempdir' - return name of temp directory
%         - 'desktop' - return whether or not the Desktop is in use
%
% FORMAT PlatFontNames = spm_platform('fonts')
% Return structure with fields named after the generic (UNIX) fonts, the
% field containing the name of the platform specific font.
%
% FORMAT PlatFontName = spm_platform('font',GenFontName)
% Map generic (UNIX) FontNames to platform specific FontNames
%
% FORMAT meminfo = spm_platform('memory',['available','total'])
% Return memory information concerning the amount of available physical
% memory or the total amount of physical memory.
%
% FORMAT PLATFORM = spm_platform
% Initialise platform specific parameters in persistent variable.
% PLATFORM - copy of persistent variable containing platform specific
% parameters.
%
% FORMAT PLATFORM = spm_platform('init')
% (Re)initialise platform specific parameters in persistent variable.
%
%--------------------------------------------------------------------------
% Since calls to spm_platform will be made frequently, most platform
% specific parameters are stored in a persistent variable.
% Subsequent calls use the information from this persistent variable, if
% it exists.
%__________________________________________________________________________

% Matthew Brett, Guillaume Flandin
% Copyright (C) 1999-2023 Wellcome Centre for Human Neuroimaging


%-Initialise
%--------------------------------------------------------------------------
persistent PLATFORM
if isempty(PLATFORM), PLATFORM = init_platform; end

if ~nargin, varargout = {PLATFORM}; return, end


switch lower(varargin{1}), case 'init'                     %-(re)initialise
%==========================================================================
PLATFORM = init_platform;
varargout = {PLATFORM};
   
case 'bigend'                     %-Return endianness for this architecture
%==========================================================================
varargout = {PLATFORM.bigend};

case 'filesys'                            %-Return file system (deprecated)
%==========================================================================
varargout = {PLATFORM.filesys};

case 'mexext'                               %-Return MEX filename extension
%==========================================================================
varargout = {PLATFORM.mexext};

case 'soext'                     %-Return shared library filename extension
%==========================================================================
varargout = {PLATFORM.soext};

case 'user'                                              %-Return user name
%==========================================================================
varargout = {PLATFORM.user};

case 'host'                                              %-Return host name
%==========================================================================
varargout = {PLATFORM.host};

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

case 'memory'                                   %-Return memory information
%==========================================================================
varargout = {meminfo(varargin{2:end})};
    
otherwise                                           %-Unknown Action string
%==========================================================================
error('Unknown Action string')

%==========================================================================
end



%==========================================================================
%- S U B - F U N C T I O N S
%==========================================================================


function PLATFORM = init_platform           %-Initialise platform variables
%==========================================================================
if strcmpi(spm_check_version,'matlab')
    comp = computer;
else
    if ismac
        comp = uname.machine;
        switch comp
            case {'x86_64'}
                comp = 'MACI64';
            case {'arm64'}
                comp = 'ARM';
            otherwise
                error('%s is not supported.',comp);
        end
    elseif isunix
        comp = uname.machine;
        switch comp
            case {'x86_64'}
                comp = 'GLNXA64';
            case {'armv6l','armv7l','armv8l','armv9l','aarch64','arm64'}
                comp = 'ARM';
            otherwise
                error('%s is not supported.',comp);
        end
    elseif ispc
        comp = 'PCWIN64';
    end
end

%-Platform definitions
%--------------------------------------------------------------------------
PDefs = {'PCWIN64',   'win',   false;...
         'MACI64',    'unx',   false;...
         'MACA64',    'unx',   false;...
         'GLNXA64',   'unx',   false;...
         'ARM',       'unx',   false};

PDefs = cell2struct(PDefs,{'computer','filesys','endian'},2);

%-Which computer?
%--------------------------------------------------------------------------
[issup, ci] = ismember(comp,{PDefs.computer});
if ~issup
    error([comp ' not supported architecture for ' spm('Ver')]);
end


%-Byte ordering
%--------------------------------------------------------------------------
PLATFORM.bigend = PDefs(ci).endian;


%-Filesystem type
%--------------------------------------------------------------------------
PLATFORM.filesys = PDefs(ci).filesys;


%-File separator character
%--------------------------------------------------------------------------
PLATFORM.sepchar = filesep;


%-MEX filename extension
%--------------------------------------------------------------------------
PLATFORM.mexext = mexext;


%-Shared library filename extension
%--------------------------------------------------------------------------
switch comp
    case {'MACI64','MACA64'}
        PLATFORM.soext = 'dylib';
    case {'GLNXA64','ARM'}
        PLATFORM.soext = 'so';
    case 'PCWIN64'
        PLATFORM.soext = 'dll';
    otherwise
        error(['Unknown platform "' comp '"']);
end


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


%-Fonts
%--------------------------------------------------------------------------
switch comp
    case {'MACI64','MACA64'}
        PLATFORM.font.helvetica = 'TrebuchetMS';
        PLATFORM.font.times     = 'Times';
        PLATFORM.font.courier   = 'Courier';
        PLATFORM.font.symbol    = 'Symbol';
    case {'GLNXA64','ARM'}
        PLATFORM.font.helvetica = 'Helvetica';
        PLATFORM.font.times     = 'Times';
        PLATFORM.font.courier   = 'Courier';
        PLATFORM.font.symbol    = 'Symbol';
    case 'PCWIN64'
        PLATFORM.font.helvetica = 'Arial Narrow';
        PLATFORM.font.times     = 'Times New Roman';
        PLATFORM.font.courier   = 'Courier New';
        PLATFORM.font.symbol    = 'Symbol';
    otherwise
        error(['Unknown platform "' comp '"']);
end


%-Desktop
%--------------------------------------------------------------------------
try
    PLATFORM.desktop = usejava('desktop');
catch
    PLATFORM.desktop = false;
end


function mem = meminfo(opt)                            %-Memory information
%==========================================================================
try
    if ispc
        % https://www.mathworks.com/help/matlab/ref/memory.html
        [uv,sv]   = memory;
        mem.avail = sv.PhysicalMemory.Available; % or uv.MemAvailableAllArrays
        mem.total = sv.PhysicalMemory.Total;
    elseif ismac
        % https://www.unix.com/man-page/osx/1/vm_stat/
        [sts,m]   = system('vm_stat'); % (page size of 4096 bytes)
        m         = strsplit(m,{':',sprintf('\n')});
        mem.avail = str2double(m{find(ismember(m,'Pages free'))+1}) * 4096;
        mem.avail = mem.avail + str2double(m{find(ismember(m,'Pages inactive'))+1}) * 4096;
        mem.total = mem.avail + str2double(m{find(ismember(m,'Pages active'))+1}) * 4096;
        mem.total = mem.total + str2double(m{find(ismember(m,'Pages speculative'))+1}) * 4096;
        mem.total = mem.total + str2double(m{find(ismember(m,'Pages wired down'))+1}) * 4096;
        mem.total = mem.total + str2double(m{find(ismember(m,'Pages occupied by compressor'))+1}) * 4096;
    else
        % http://man7.org/linux/man-pages/man5/proc.5.html
        m         = strsplit(fileread('/proc/meminfo')); % (in kB)
        mem.avail = str2double(m{find(ismember(m,'MemAvailable:'))+1}) * 1024;
        mem.total = str2double(m{find(ismember(m,'MemTotal:'))+1}) * 1024;
    end
catch
    mem = struct('avail',NaN,'total',NaN);
end

if ~nargin, return, end

switch lower(opt)
    case 'total'
        mem = mem.total;
    case 'available'
        mem = mem.avail;
    otherwise
        error('Unknown memory option.');
end
