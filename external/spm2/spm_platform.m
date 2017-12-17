function varargout=spm_platform(varargin)
% Platform specific configuration parameters for SPM
%
% FORMAT ans = spm_platform(arg)
% arg  - optional string argument, can be
%        - 'bigend'  - return whether this architecture is bigendian
%                      - Inf - is not IEEE floating point
%                      - 0   - is little end
%                      - 1   - big end
%        - 'filesys' - type of filesystem
%                      - 'unx' - UNIX
%                      - 'win' - DOS
%                      - 'mac' - Macintosh
%                      - 'vms' - VMS
%        - 'sepchar' - returns directory separator
%        - 'rootlen' - returns number of chars in root directory name
%        - 'user'    - returns username
%        - 'tempdir' - returns name of temp directory
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
% comp         - computer to use [defaults to MatLab's `computer`]
% PLATFORM - copy of persistent PLATFORM
%
% FORMAT spm_platform
% Initialises platform specific parameters in persistent PLATFORM
% (External gateway to init_platform(computer) subfunction)
%
%                           ----------------
% SUBFUNCTIONS:
%
% FORMAT init_platform(comp)
% Initialise platform specific parameters in persistent PLATFORM
% comp         - computer to use [defaults to MatLab's `computer`]
%
%-----------------------------------------------------------------------
%
% Since calls to spm_platform will be made frequently, most platform
% specific parameters are stored as a structure in the persistent variable
% PLATFORM. Subsequent calls use the information from this persistent
% variable, if it exists.
%
% Platform specific difinitions are contained in the data structures at
% the beginning of the init_platform subfunction at the end of this
% file.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Matthew Brett
% $Id$



%-Initialise
%-----------------------------------------------------------------------
persistent PLATFORM
if isempty(PLATFORM), PLATFORM = init_platform; end

if nargin==0, return, end


switch lower(varargin{1}), case 'init'                  %-(re)initialise
%=======================================================================
init_platform(varargin{2:end})
varargout = {PLATFORM};
   
case 'bigend'                      %-Return endian for this architecture
%=======================================================================
varargout = {PLATFORM.bigend};
if ~isfinite(PLATFORM.bigend),
	if isnan(PLATFORM.bigend)
		error(['I don''t know if "',computer,'" is big-endian.'])
	else
		error(['I don''t think that "',computer,...
			'" uses IEEE floating point ops.'])
	end
end

case 'filesys'                                      %-Return file system
%=======================================================================
varargout = {PLATFORM.filesys};

case 'sepchar'                         %-Return file separator character
%=======================================================================
warning('use filesep instead (supported by MathWorks)')
varargout = {PLATFORM.sepchar};

case 'rootlen'           %-Return length in chars of root directory name 
%=======================================================================
varargout = {PLATFORM.rootlen};

case 'user'                                         %-Return user string
%=======================================================================
varargout = {PLATFORM.user};

case 'tempdir'                              %-Return temporary directory
%=======================================================================
twd = getenv('SPMTMP');
if isempty(twd)
	twd = tempdir;
end 
varargout = {twd};


case {'font','fonts'}    %-Map default font names to platform font names
%=======================================================================
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

otherwise                                        %-Unknown Action string
%=======================================================================
error('Unknown Action string')

%=======================================================================
end



%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================


function PLATFORM = init_platform(comp)             %-Initialise platform variables
%=======================================================================
if nargin<1, comp=computer; end

%-Platform definitions
%-----------------------------------------------------------------------
PDefs = {	'PCWIN',	'win',	0;...
		'MAC',		'unx',	1;...
		'MACI64',	'unx',	0;...
		'SUN4',		'unx',	1;...
		'SOL2',		'unx',	1;...
		'HP700',	'unx',	1;...
		'SGI',		'unx',	1;...
		'SGI64',	'unx',	1;...
		'IBM_RS',	'unx',	1;...
		'ALPHA',	'unx',	0;...
		'AXP_VMSG',	'vms',	Inf;...
		'AXP_VMSIEEE',	'vms',	0;...
		'LNX86',	'unx',	0;...
		'GLNX86',	'unx',  0;...
		'GLNXA64',      'unx',  0;...
		'VAX_VMSG',	'vms',	Inf;...
		'VAX_VMSD',	'vms',	Inf	};

PDefs = cell2struct(PDefs,{'computer','filesys','endian'},2);


%-Which computer?
%-----------------------------------------------------------------------
ci = find(strcmp({PDefs.computer},comp));
if isempty(ci), error([comp,' not supported architecture for SPM']), end


%-Set bigend
%-----------------------------------------------------------------------
PLATFORM.bigend = PDefs(ci).endian;


%-Set filesys
%-----------------------------------------------------------------------
PLATFORM.filesys = PDefs(ci).filesys;


%-Set filesystem dependent stuff
%-----------------------------------------------------------------------
%-File separators character
%-Length of root directory strings
%-User name finding
%-(mouse button labels?)
switch (PLATFORM.filesys)
case 'unx'
	PLATFORM.sepchar = '/';
	PLATFORM.rootlen = 1;
	PLATFORM.user    = getenv('USER');
case 'win'
	PLATFORM.sepchar = '\';
	PLATFORM.rootlen = 3;
	PLATFORM.user    = getenv('USERNAME');
	if isempty(PLATFORM.user)
		PLATFORM.user = spm_win32utils('username'); end
otherwise
	error(['Don''t know filesystem ',PLATFORM.filesys])
end

%-Fonts
%-----------------------------------------------------------------------
switch comp
case {'SOL2'}	%-Some Sol2 platforms give segmentation violations with Helvetica
	PLATFORM.font.helvetica = 'Lucida';
	PLATFORM.font.times     = 'Times';
	PLATFORM.font.courier   = 'Courier';
	PLATFORM.font.symbol    = 'Symbol';
case {'SUN4','SOL2','HP700','SGI','SGI64','IBM_RS','ALPHA','LNX86','GLNX86','GLNXA64','MAC'}
	PLATFORM.font.helvetica = 'Helvetica';
	PLATFORM.font.times     = 'Times';
	PLATFORM.font.courier   = 'Courier';
	PLATFORM.font.symbol    = 'Symbol';
case {'PCWIN'}
	PLATFORM.font.helvetica = 'Arial Narrow';
	PLATFORM.font.times     = 'Times New Roman';
	PLATFORM.font.courier   = 'Courier New';
	PLATFORM.font.symbol    = 'Symbol';
end
