function metadata = read_mef_session_metadata_3p0(this, varargin)
% MEFSESSION_3P0.READ_MEF_SESSION_METADATA_3P0 Retrieves the session metadata from a MEF 3.0 session
%   
% Syntax:
%   metadata = read_mef_session_metadata_3p0(this)
%   metadata = __(__, sess_path)
%   metadata = __(__, sess_path, password)
%   metadata = __(__, sess_path, password, map_indices)
% 
% Input(s):
%   this            - [obj] MEFSession_3p0 object
%   sess_path     	- [str] (opt) path (absolute or relative) to the MEF3 
%                     session folder (default = path in this)
%   password        - [struct] password structure to the MEF3 data (default
%                     = password in this)
%   map_indices     - [logical] flag whether indices should be mapped [true
%                     or false] (default = true)
%
% Output(s): 
%   metadata    	- [struct] structure containing session metadata, 
%                     channels metadata, segments metadata and records
%
% Note:
% 
% See also read_mef_session_metadata.

%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo
%   Clinic, Rochester MN) <https://github.com/MaxvandenBoom/matmef>.
%   Adapted from PyMef (by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and
%   Dan Crepeau)
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. This program is distributed in the hope that
%   it will be useful, but WITHOUT ANY WARRANTY; without even the implied
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
%   the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <https://www.gnu.org/licenses/>.
%

% Copyright 2020 Richard J. Cui. Adapted: Sat 02/01/2020 10:30:50.708 PM
% $Revision: 0.5 $  $Date: Thu 04/09/2020 11:05:29.605 PM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
sess_path = q.sess_path;
password = q.password;
map_indices = q.map_indices;

if isempty(sess_path)
    sess_path = this.SessionPath;
end % if

if isempty(password)
    password = this.Password;
end % if

% =========================================================================
% main
% =========================================================================
pw = this.processPassword(password);
metadata = read_mef_session_metadata(sess_path, pw, map_indices);

end % function

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_sp = ''; % sesseion path
default_pw = struct([]);
default_mi = true;

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('sess_path', default_sp, @isstr);
p.addOptional('password', default_pw, @isstruct);
p.addOptional('map_indices', default_mi, @islogical);

% parse inputs and return results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]