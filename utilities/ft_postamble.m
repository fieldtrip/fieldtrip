function ft_postamble(cmd, varargin)

% FT_POSTAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the end
% of the function.
%
% This ft_postamble m-file is a function, but internally it executes a number of
% private scripts in the callers workspace. This allows the private script to access
% the variables in the callers workspace and behave as if the script were included as
% a header file in C-code.
%
% See also FT_PREAMBLE

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% ideally this would be a script, because the local variables would then be
% shared with the calling function. Instead, this is a function which then
% passes the variables explicitely to another script which is eval'ed.

% the following section ensures that these scripts are included as
% dependencies when using the MATLAB compiler
%
%#function ft_postamble_debug
%#function ft_postamble_provenance
%#function ft_postamble_previous
%#function ft_postamble_history
%#function ft_postamble_savevar
%#function ft_postamble_savefig
%#function ft_postamble_randomseed
%#function ft_postamble_hastoolbox

% this is a trick to pass the input arguments into the ft_postamble_xxx script
assignin('caller', 'postamble_argin', varargin);

full_cmd = ['ft_postamble_' cmd];
cmd_exists = false;

if exist(full_cmd, 'file')
  % Matlab can find commands in a private subdirectory; Octave cannot.
  % If pwd is already the private directory, or if using Matlab then
  % the command can be evaluated directly
  cmd_exists = true;

elseif ~ft_platform_supports('exists-in-private-directory')
  % Octave does not find files by name in a private directory, so the full filename must be specified.
  private_dir = fullfile(fileparts(which(mfilename)),'private');
  full_path   = fullfile(private_dir,[full_cmd '.m']);
  cmd_exists  = exist(full_path, 'file');
  full_cmd_parts = {'ft_tmp_orig_pwd = pwd(); ',...
                    'ft_tmp_orig_pwd_cleaner = onCleanup(@()cd(ft_tmp_orig_pwd)); ',...
                    sprintf('cd(''%s''); ',private_dir),...
                    [full_cmd '; '],...
                    'clear ft_tmp_orig_pwd_cleaner;'};
  full_cmd = sprintf('%s', full_cmd_parts{:});
end

if ~cmd_exists
  % XXX earlier versions would not do anything if ~cmd_exists,
  % but fail silently (without raising an error or warning).
  % Should that behavior be kept?
  ft_error('Could not run %s - does not seem to exist', full_cmd);
end

evalin('caller', full_cmd);
