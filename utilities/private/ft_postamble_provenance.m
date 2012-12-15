% FT_POSTAMBLE_PROVENANCE is a helper script that reports the time and memory that was
% used by the calling function and that stores this information together with user
% name, MATLAB version and other provenance information in the output cfg structure.
% This should be used together with FT_PREAMBLE_PROVENANCE, which records the time and
% memory at the start of the function.
%
% Another aspects of provenance relates to uniquely identifying the input and the
% output data. The code that deals with tracking the information about the input data
% structures is found in ft_preamble_loadvar. The code that deals with tracking the
% information about the output data structures is found in ft_preamble_history.

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% the name of the variables are passed in the preamble field
global ft_default

if isfield(cfg, 'trackcallinfo') && ~istrue(cfg.trackcallinfo)
  % do not track the call information
  return
end

stack = dbstack('-completenames');
% stack(1) is this script
% stack(2) is the calling ft_postamble function
% stack(3) is the main FieldTrip function that we are interested in
stack = stack(3);

% the proctime, procmem and calltime rely on three cryptical variables that were
% created and added to the function workspace by the ft_preamble_callinfo script.
cfg.callinfo.proctime = toc(ftohDiW7th_FuncTimer);
cfg.callinfo.procmem  = memtoc(ftohDiW7th_FuncMem);

if istrue(ft_getopt(cfg, 'showcallinfo', 'yes'))
  % print some feedback on screen, this is meant to educate the user about
  % the requirements of certain computations and to use that knowledge in
  % distributed computing
  
  stack = dbstack('-completenames');
  % stack(1) is this script
  % stack(2) is the calling ft_postamble function
  % stack(3) is the main FieldTrip function that we are interested in
  stack = stack(3);
  
  if ispc()
    % don't print memory usage info under Windows; this does not work (yet)
    fprintf('the call to "%s" took %d seconds\n', stack.name, round(cfg.callinfo.proctime));
  else
    fprintf('the call to "%s" took %d seconds and required the additional allocation of an estimated %d MB\n', stack.name, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));
  end
  
  clear stack
end

% compute the MD5 hash of each of the output arguments
if isequal(ft_default.postamble, {'varargout'})
  cfg.callinfo.outputhash = cellfun(@CalcMD5, cellfun(@mxSerialize, varargout, 'UniformOutput', false), 'UniformOutput', false);
else
  cfg.callinfo.outputhash = cellfun(@CalcMD5, cellfun(@mxSerialize, cellfun(@eval, ft_default.postamble, 'UniformOutput', false), 'UniformOutput', false), 'UniformOutput', false);
end

clear ftohDiW7th_FuncTimer
clear ftohDiW7th_FuncMem
clear stack

