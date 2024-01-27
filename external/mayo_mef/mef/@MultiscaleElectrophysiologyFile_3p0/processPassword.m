function pw = processPassword(this, varargin)
% MultiscaleElectrophysiologyFile_3p0.processPassword process password of MEF 3.0 data
% 
% Syntax:
%   pw = processPassword(this)
%   pw = __(__, password)
%   pw = __(__, 'Level1Password', level_1_pw)
%   pw = __(__, 'Level2Password', level_2_pw)
%   pw = __(__, 'AccessLevel', access_level)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile_3p0 object
%   password        - [struct] (opt) password structure
%                     .Level1Password
%                     .Level2Password
%                     .AccessLevel
%   level_1_pw      - [str] (para) password of level 1 (default = '')
%   level_2_pw      - [str] (para) password of level 2 (default = '')
%   access_level    - [str] (para) data decode level to be used
%                     (default = 1)
% 
% Output(s):
%   pw              - [str] password at the selected level
% 
% Note:
%   if 'password' is provided, the parameters will be ignored.
% 
% See also .

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
password = q.password;
level_1_pw = q.Level1Password;
level_2_pw = q.Level2Password;
access_level = q.AccessLevel;

% =========================================================================
% main
% =========================================================================
if isempty(password)
    if isempty(level_1_pw)
        level_1_pw = this.Level1Password;
    end % if
    if isempty(level_2_pw)
        level_2_pw = this.Level2Password;
    end % if
    if isempty(access_level)
        access_level = this.AccessLevel;
    end % if
else
    level_1_pw = password.Level1Password;
    level_2_pw = password.Level2Password;
    access_level = password.AccessLevel;
end % if

switch access_level
    case 1
        pw = level_1_pw;
    case 2
        pw = level_2_pw;
    otherwise
        error(sprintf('%s:%s:invalidAccessLevel', mfilename('class'), mfilename),...
            'invalid access level; access level must be 1 or 2')
end % switch

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_stpw = struct([]);
default_pw = '';
default_al = 1; % access level

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('password', default_stpw, @isstruct);
p.addParameter('Level1Password', default_pw, @isstr);
p.addParameter('Level2Password', default_pw, @isstr);
p.addParameter('AccessLevel', default_al, @isnumeric);

% pase and return results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]