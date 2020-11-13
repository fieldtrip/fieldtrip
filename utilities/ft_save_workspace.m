function ft_save_workspace(dae1Ot8o_dir)

% FT_SAVE_WORKSPACE saves every variable in the base workspace to a .mat file with
% the same name as the variable in the workspace itself. For example, the variable
% "ans" would be saved to the file "ans.mat". Prior to calling this function, you
% might want to clean up your workspace using CLEAR or KEEP.
%
% Use as
%   ft_save_workspace(dirname)
%
% If the directory does not yet exist, this function will create it for you. If you
% leave it empty, the files will be saved to the present working directory.
%
% For example, the following will save all variables to a time-stamped
% sub-directory that is created inside the present working directory:
%
%   ft_save_workspace(datestr(now))
%
% See also SAVE, LOAD, SAVEFIG, CLEAR, KEEP

% Copyright (C) 2018, Robert Oostenveld
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

% use long and random variable names to avoid clashes with variables in the base workspace
if nargin<1
  dae1Ot8o_dir = '.';
end

if ~isdir(dae1Ot8o_dir)
  ft_info('creating directory ''%s''\n', dae1Ot8o_dir);
  mkdir(dae1Ot8o_dir);
end

% get a list of all the variables in the base workspace
dae1Ot8o_w = evalin('base', 'whos');

% copy each of the variables over to the local workspace
for dae1Ot8o_i=1:numel(dae1Ot8o_w)
  dae1Ot8o_var = dae1Ot8o_w(dae1Ot8o_i).name;
  dae1Ot8o_val = evalin('base', dae1Ot8o_var);
  eval(sprintf('%s = dae1Ot8o_val;', dae1Ot8o_var));
  dae1Ot8o_file = fullfile(dae1Ot8o_dir, [dae1Ot8o_var '.mat']);
  ft_info('saving variable ''%s'' to ''%s''\n', dae1Ot8o_var, dae1Ot8o_file);
  save(dae1Ot8o_file, dae1Ot8o_var);
end
