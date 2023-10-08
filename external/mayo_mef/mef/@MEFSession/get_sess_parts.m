function [path_to_sess, sess_name, sess_ext] = get_sess_parts(this, varargin)
% MEFSESSION.GET_SESSEXT get the parts of session path
% 
% Syntax:
%   [path_to_sess, sess_name, sess_ext] = get_sess_parts(this)
%   __ = __(__, sess_path)
% 
% Input(s):
%   this            - [obj] MEFSession object
%   sess_path       - [str] (opt) Session path, including session name and
%                     session extension (default: this.SessionPath)
% 
% Output(s):
%   path_to_sess    - [str] path up to the session, not include session
%                     name and extension
%   sess_name       - [str] session name, no extension
%   sess_ext        - [str] session extension, including '.'
% 
% See also .

% Copyright 2020 Richard J. Cui. Created: Thu 02/20/2020 11:24:47.068 AM
% $Revision: 0.2 $  $Date: Tue 02/21/2023 11:05:42.802 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
sess_path = q.sess_path;
if isempty(sess_path)
    sess_path = this.SessionPath;
end % if

% =========================================================================
% main
% =========================================================================
[path_to_sess, sess_name, sess_ext] = fileparts(sess_path);

% update
this.PathToSession = path_to_sess;
this.SessionName = sess_name;
this.SessionExt = sess_ext;

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_sp = '';

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('sess_path', default_sp, @ischar);

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]