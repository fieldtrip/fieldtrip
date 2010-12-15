function varargout = ft_wizard(wizard_filename)

% FT_WIZARD will evaluate a FieldTrip analysis script in steps, allowing you
% to go to the next step if you are content with the data so far, or to the
% previous step if you want to repeat it with different configuration
% settings.
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
% the main Matlab workspace. 
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
% $Log: wizard.m,v $
% Revision 1.3  2008/09/22 20:17:44  roboos
% added call to ft_defaults to the begin of the function
%
% Revision 1.2  2007/05/14 08:30:56  roboos
% renamed wizard_gui to wizard_base, updated help
%
% Revision 1.1  2007/05/10 09:06:21  roboos
% initial version
%

ft_defaults

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
