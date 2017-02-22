function varargout = spm_check_version(tbx,chk)
% Check a version number against a Toolbox version
%
% FORMAT tbx = spm_check_version
% tbx    - toolbox name {'matlab','octave',...}
%
% FORMAT v = spm_check_version(tbx)
% tbx    - toolbox name {'matlab','octave','signal',...}
%
% v      - version number {string}
%
% FORMAT status = spm_check_version(tbx,chk)
% tbx    - toolbox name {'matlab','octave','signal',...}
% chk    - version number to be checked {string}
%
% status - outcome of the comparison:
%          -1: Toolbox version is earlier than the user supplied version
%           0: Toolbox and user versions are the same
%           1: Toolbox version is later than the user supplied version
%              Think of it this way, the sign of status is determined by
%              MATLAB_TOOLBOX_VERSION - USER_VERSION (i.e., THE VERSION YOU
%              INPUT).
%__________________________________________________________________________
%
% This function checks if a user supplied version number is less than,
% equal to or greater than the version number of specified toolbox. If no
% toolbox is specified the function checks the version of MATLAB. User
% defined toolboxes can be checked but the Contents.m file must conform
% to the specification defined in ver.m
%
% This function assumes that the version number is really a text string
% with fields major.minor.revision.build. Other formats are not supported.
% Checking is done to the level specified by the input version. Thus an
% input of '7' will be rated as the same version as 7.1, but 7.0 would be
% rated as earlier than 7.1.
%__________________________________________________________________________
%
% EXAMPLES:
%
% If the MATLAB version is 7.1.0.83, and the user supplied version is '7':
% status = spm_check_version('matlab','7');
% returns status == 0   : major revision numbers are the same.
%
% If the MATLAB version is 7.1.0.0, and the user supplied version is '7.1':
% status = spm_check_version('matlab','7');
% returns status == 0   : major and minor revision numbers are the same.
%
% If the MATLAB version is 7.1.0.83, and the user supplied version is '7.2':
% status = spm_check_version('matlab','7.2');
% returns status == -1   : major + minor revision is earlier for MATLAB.
%
% If the MATLAB version is 6.5.1, and the user supplied version is '6.5.0'.
% status = spm_check_version('matlab','6.5.0');
% returns status == 1     : major + minor + release revision is later
%                           for MATLAB
% The statement ( spm_check_version('matlab','6.5.0') > 0 ) is true for
% all MATLAB Toolbox versions after 6.5.0.
%__________________________________________________________________________
%
% See also VERSION, VER, VERLESSTHAN.
%__________________________________________________________________________
% Copyright (C) 2006-2014 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman
% $Id: spm_check_version.m 6156 2014-09-05 17:34:53Z guillaume $


%-Detect software used
%==========================================================================
if ~nargin || isempty(tbx)
    if exist('OCTAVE_VERSION','builtin')
        tbx = 'octave';
    else
        tbx = 'matlab';
    end
    if ~nargin, varargout = {tbx}; return; end
end


%-Get the requested toolbox version
%==========================================================================
if strcmpi(tbx,'matlab')
    % MATLAB is a special case. The ver command does not report the
    % entire version string.
    tbxVer = strtok(builtin('version'));
elseif strcmpi(tbx,'octave')
    tbxVer = version;
    tbxVer = strrep(tbxVer,'+','');
else
    tbxStruct = ver(tbx);
    if isempty(tbxStruct)
        error('Cannot find given toolbox.');
    elseif numel(tbxStruct) > 1
        error('Too many toolboxes found for given toolbox name.')
    end
    tbxVer = tbxStruct.Version;
end

if nargin == 1, varargout = {tbxVer}; return; end


%-Parse user supplied version number
%==========================================================================

if strcmpi(tbx,'matlab') && strcmpi(spm_check_version,'octave')
    varargout = {1}; % hack
    return;
end

% If a number is supplied then convert to text
%--------------------------------------------------------------------------
if isnumeric(chk)
    chk = num2str(chk);
end

% If too many fields in input then error
%--------------------------------------------------------------------------
if numel(strfind(chk,'.')) > 3
    error(['Input string has too many fields. ',...
        'Only major.minor.release.build fields are supported.']);
end

% parse the input string into a cell array of strings
%--------------------------------------------------------------------------
vCHK = textscan(chk,'%s %s %s %s','delimiter','.');

% find empty cells. These will be set to 0 in both cell arrays so that the
% same revision fields are compared
%--------------------------------------------------------------------------
emptyCHK = cellfun('isempty',vCHK);

% parse the requested toolbox version string into a cell array of strings.
%--------------------------------------------------------------------------
vTBX = textscan(tbxVer,'%s %s %s %s','delimiter','.');

% find empty cells. These will be removed from the cell array so that the
% same revision fields are compared
%--------------------------------------------------------------------------
emptyTBX = cellfun('isempty',vTBX);

% combine empty cell indices
%--------------------------------------------------------------------------
emptyCells = find(emptyCHK | emptyTBX);

% and remove them.
%--------------------------------------------------------------------------
vCHK(emptyCells) = [];
vTBX(emptyCells) = [];

% final version fields to compare are converted to numbers
%--------------------------------------------------------------------------
vCHKMat = str2num(strvcat([vCHK{:}]))';
vTBXMat = str2num(strvcat([vTBX{:}]))';


%-Compare versions
%==========================================================================
% This array will be used to decide which version is the later one. The
% differences between versions are computed on each element of the version
% number: major.minor.release.build and then multiplied by the base2array
% and summed. This is equivalent to logically combining the base 2 place
% values and ensures that the fields are evaluated by their significance.
% The sign of the summed values gives the result. Positive values mean the
% toolbox version is later, negative values mean the user supplied version,
% is later, and 0 means they are equal.
%--------------------------------------------------------------------------
base2array             = [8 4 2 1];
base2array(emptyCells) = [];

% take the difference between version elements.
%--------------------------------------------------------------------------
vDiff = vTBXMat - vCHKMat;

% get the sign of the differences
%--------------------------------------------------------------------------
signDiff = sign(vDiff);

% multiply by the base2array
%--------------------------------------------------------------------------
vVal = signDiff.*base2array;

% If the sign of the sum is positive then toolbox version is later; if 0
% they are the same, and if negative the user supplied version is later.
% Think of it this way: MATLAB_TOOLBOX_VERSION - USER_VERSION.
% Remember the comparison is made with respect to the FIELDS SUPPLED IN THE
% INPUT
%--------------------------------------------------------------------------
switch sign(sum(vVal))
    case -1     % MATLB_TOOLBOX_VERSION is EARLIER than USER_VERSION
        status = -1;
    case 0      % MATLB_TOOLBOX_VERSION = EQUALS = USER_VERSION
        status =  0;
    case 1      % MATLB_TOOLBOX_VERSION is LATER than USER_VERSION
        status =  1;
end

varargout = {status};
