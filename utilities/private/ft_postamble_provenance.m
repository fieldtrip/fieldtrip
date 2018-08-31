% FT_POSTAMBLE_PROVENANCE is a helper script that reports the time and memory used by
% the calling function and that stores this information together with user name,
% MATLAB version and other provenance information in the output cfg structure. This
% script is to be used together with FT_PREAMBLE_PROVENANCE, which records the time
% and memory at the start of the function.
%
% FieldTrip also attempts to uniquely identify the input and the output data. The
% code that deals with tracking the input data structures is found in
% FT_PREAMBLE_LOADVAR. The code that deals with tracking the information about the
% output data structures is found in FT_POSTAMBLE_HISTORY.
%
% Use as
%   ft_preamble provenance
%   .... regular code goes here ...
%   ft_postamble provenance
%
% See also FT_PREAMBLE_PROVENANCE

% Copyright (C) 2011-2016, Robert Oostenveld, DCCN
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

if isfield(cfg, 'trackcallinfo') && ~istrue(cfg.trackcallinfo)
  % do not track the call information
  return
end

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
end % if showcallinfo=yes
clear stack

if isfield(cfg, 'trackdatainfo') && istrue(cfg.trackdatainfo)
  % compute the MD5 hash of each of the output arguments
  % temporarily remove the cfg field for getting the hash (creating a duplicate of the data, but still has the same mem ref, so no extra mem needed)
  if isequal(iW1aenge_postamble, {'varargin'})
    tmpargout = varargout;
  else
    tmpargout = cellfun(@eval, iW1aenge_postamble, 'UniformOutput', false);
  end
  cfg.callinfo.outputhash = cell(1,numel(tmpargout));
  for iargout = 1:numel(tmpargout)
    tmparg = tmpargout{iargout}; % can't get number of bytes with whos unless taken out of it's cell
    if isfield(tmparg,'cfg')
      tmparg = rmfield(tmparg, 'cfg');
    end
    % only calculate md5 when below 2^31 bytes (CalcMD5 can't handle larger input)
    bytenum = whos('tmparg');
    bytenum = bytenum.bytes;
    if bytenum<2^31
      try
        cfg.callinfo.outputhash{iargout} = ft_hash(tmparg);
      catch
        % the mxSerialize mex file is not available on all platforms, do not compute a hash
        % http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2452
      end
    else
      % the data is too large, do not compute a hash
    end
  end
  clear tmpargout iargout tmparg bytenum; % remove the extra references
end

clear ftohDiW7th_FuncTimer
clear ftohDiW7th_FuncMem

