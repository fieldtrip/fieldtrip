%
%   Retrieves the session metadata from a MEF3 file
%   
%   [metadata] = read_mef_session_metadata(fileName, password)
%	
%   sessionPath     = path (absolute or relative) to the MEF3 session folder
%   password        = password to the MEF3 data; Pass empty string/variable if not encrypted
%   mapIndices      = flag whether indices should be mapped [0 or 1; default is 0]
%
%   Returns: 
%       metadata    = structure containing session metadata, channels metadata, segments metadata and records
%
%
%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)
%   Adapted from PyMef (by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and Dan Crepeau)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function metadata = read_mef_session_metadata(sessionPath, password, mapIndices)
