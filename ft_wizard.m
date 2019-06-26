function [varargout] = ft_wizard(wizard_filename)

% FT_WIZARD is a graphical user interface to evaluate a FieldTrip analysis
% script one step at a time, allowing you to go to the next step if you are
% content with the data so far, or to the previous step if you want to repeat it
% with different configuration settings.
%
% Use as
%   ft_wizard scriptname
% or 
%   ft_wizard('scriptname')
%
% Use the functional form of FT_WIZARD, such as FT_WIZARD('scriptname'), when
% the name of the script is stored in a string, when an output argument is
% requested, or if the name of the script contains spaces. If you do not
% specify an output argument, the results will be stored as variables in
% the main MATLAB workspace. 
% 
% Besides the buttons, you can use the following key combinations
%   Ctrl-O        load a new script from a file
%   Ctrl-S        save the script to a new file
%   Ctrl-E        open the current script in editor
%   Ctrl-P        go to previous step
%   Ctrl-N        go to next step
%   Ctrl-Q        quit, do not save the variables
%   Ctrl-X        exit, save the variables to the workspace
% 
% See also FT_ANALYSISPROTOCOL

% Copyright (C) 2007-2010, Robert Oostenveld
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

% this function is a wrapper around the actual GUI function
% the main purpose of this function is to provide a workspace where the
% intermediate results can be stored

% the wizard_ok variable will be set to 0/1 by the wizard_base function
wizard_ok = 0;

if nargin>0
  wizard_fig = wizard_base(wizard_filename);
else
  wizard_fig = wizard_base;
end

waitfor(wizard_fig);

if wizard_ok
  % get the intermediate variables, do not export the wizard bookkeeping
  % variables (which all start with "wizard_")
  wizard_var = whos;
  wizard_var(strmatch('wizard_', {wizard_var.name})) = [];
  
  if nargout==0
    % store the results in the BASE workspace
    for wizard_i=1:length(wizard_var)
      assignin('base', wizard_var(wizard_i).name, eval(wizard_var(wizard_i).name));
    end
  else
    % store the results in an output structure
    varargout{1} = [];
    for wizard_i=1:length(wizard_var)
      varargout{1}.(wizard_var(wizard_i).name) = eval(wizard_var(wizard_i).name);
    end
  end % if nargout
end % if ok

return % main function 
