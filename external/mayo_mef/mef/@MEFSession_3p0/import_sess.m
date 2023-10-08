function [X, t] = import_sess(this, varargin)
% MEFSESSION_3P0.IMPORT_SESS import session of MEF 3.0 data
% 
% Syntax:
%   [X, t] = import_sess(this, begin_stop, bs_unit, sel_chan)
%   [X, t] = import_sess(__, pw)
% 
% Input(s):
%   this            - [obj] MEFSession_2p1 object
%   begin_stop      - [num] 1 x 2 array of begin and stop points of
%                     importing the session, absolute time points
%   bs_unit         - [str] unit of begin_stop: 'uUTC','Index', 'Second', 
%                     'Minute', 'Hour', and 'Day'.
%   sel_chan        - [str array] the names of the selected channels
%   pw              - [struct] (para) password structure
%                     .Session      : session password
%                     .Subject      : subject password
%                     .Data         : data password
% 
% Output(s):
%   X               - [num array] M x N array, where M is the number of
%                     channels and N is the number of signals extracted
%   t               - [num] 1 x N array, time indeces of the signals
% 
% Note:
%   Import data from different channels of the session.
% 
% See also importSignal, importSession.

% Copyright 2020 Richard J. Cui. Created: Thu 02/06/2020  3:40:19.634 PM
% $Revision: 0.2 $  $Date: Thu 02/06/2020  7:50:20.972 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
begin_stop = q.begin_stop;
bs_unit = q.bs_unit;
sel_chan = q.sel_chan;
pw = q.pw;

if isempty(pw)
    pw = this.Password;
end % if

sess_path = this.SessionPath;

% =========================================================================
% main
% =========================================================================
num_chan = numel(sel_chan); % number of selected channels
X = [];
for k = 1:num_chan
    fn_k = convertStringsToChars(sel_chan(k) + ".timd"); % filename of channel k
    [x_k, t] = this.importSignal(begin_stop, bs_unit, sess_path, fn_k,...
        'Level1Password', pw.Level1Password,...
        'Level2Password', pw.Level2Password,...
        'AccessLevel', pw.AccessLevel);
    x_k = x_k(:).'; % make sure it is a horizontal vector
    
    X = cat(1, X, x_k);
end % for

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(this, varargin)

% defaults
default_pw = struct([]); % password

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('begin_stop', @(x) isnumeric(x) & numel(x) == 2 & x(1) <= x(2));
p.addRequired('bs_unit', @isstr);
p.addRequired('sel_chan', @isstring) % must be string array
p.addOptional('pw', default_pw, @isstruct);

% parse and return the results
p.parse(this, varargin{:});
q = p.Results;

end % function

% [EOF]