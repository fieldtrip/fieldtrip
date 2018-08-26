% FT_PREAMBLE_PROVENANCE is a helper script that records the time and memory at the
% start of the function. At the end of the function FT_POSTAMBLE_PROVENANCE will
% record and store the time and memory in the output configuration together with
% information about the environment, such as the name of the user and computer, the
% MATLAB and FieldTrip version, etc.
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
% See also FT_POSTAMBLE_PROVENANCE

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

% Record the start time and memory. These are used by ft_postamble_callinfo, which
% stores them in the output cfg.callinfo.  In the mean time, they are stored in the
% function workspace, which is why they should have cryptical names to prevent any
% variable name clashes.

if (isfield(cfg, 'trackcallinfo') && ~istrue(cfg.trackcallinfo))
  % do not track the call information
  return
end

% add the user-specified cfg (before any defaults handling etc.) to the callinfo
% some fields are for internal use only and should not be stored
cfg.callinfo.usercfg = removefields(cfg, ignorefields('provenance'));

if isfield(cfg, 'trackdatainfo') && istrue(cfg.trackdatainfo)
  % compute the MD5 hash of each of the input arguments
  % temporarily remove the cfg field for getting the hash (creating a duplicate of the data, but still has the same mem ref, so no extra mem needed)
  if isequal(iW1aenge_preamble, {'varargin'})
    tmpargin = varargin;
  else
    isvar = cellfun(@(x) exist(x, 'var')==1, iW1aenge_preamble);
    tmpargin = cellfun(@eval, iW1aenge_preamble(isvar), 'UniformOutput', false);
    tmpargin( isvar) = tmpargin;
    tmpargin(~isvar) = {[]};
    clear isvar
  end
  cfg.callinfo.inputhash = cell(1,numel(tmpargin));
  for iargin = 1:numel(tmpargin)
    tmparg = tmpargin{iargin}; % can't get number of bytes with whos unless taken out of it's cell
    if isfield(tmparg,'cfg')
      tmparg = rmfield(tmparg, 'cfg');
    end
    % only calculate md5 when below 2^31 bytes (CalcMD5 can't handle larger input)
    bytenum = whos('tmparg');
    bytenum = bytenum.bytes;
    if bytenum<2^31
      try
        cfg.callinfo.inputhash{iargin} = ft_hash(tmparg);
      catch
        % the mxSerialize mex file is not available on all platforms, do not compute a hash
        % http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2452
      end
    else
      % the data is too large, do not compute a hash
    end
  end
  clear tmpargin iargin tmparg bytenum; % remove the extra references
end

stack = dbstack('-completenames');
% stack(1) is this script
% stack(2) is the calling ft_postamble function
% stack(3) is the main FieldTrip function that we are interested in
stack = stack(3);

% add information about the FieldTrip and MATLAB version used to the configuration
try
  cfg.callinfo.fieldtrip = ft_version();
catch
  cfg.callinfo.fieldtrip = 'unknown';
end
cfg.callinfo.matlab    = version();
cfg.callinfo.computer  = lower(computer); % for example maci64, glnx86, ...

% add information about the execution environment to the configuration
cfg.callinfo.hostname = gethostname();
cfg.callinfo.user     = getusername();
cfg.callinfo.pwd      = pwd;
cfg.callinfo.calltime = clock();

% add information about the function filename and revision to the configuration
cfg.version.name = stack.file;
clear stack

% the revision number is maintained by SVN in the ft_revision variable in the calling function
if ~exist('ft_revision', 'var')
  cfg.version.id   = 'unknown';
else
  cfg.version.id   = ft_revision;
end

ftohDiW7th_FuncTimer = tic();
ftohDiW7th_FuncMem   = memtic();
