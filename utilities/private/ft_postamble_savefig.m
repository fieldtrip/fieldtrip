% FT_POSTAMBLE_SAVEVAR is a helper script that optionally saves the output
% FieldTrip data structures to a *.mat file on disk. This is useful for
% batching and for distributed processing. This makes use of the
% cfg.outputfile variable.
%
% Use as
%   ft_postamble savevar data
%   ft_postamble savevar source mri
%
% See also FT_PREAMBLE, FT_POSTAMBLE, FT_POSTAMBLE_SAVEVAR

% Copyright (C) 2019, Robert Oostenveld, DCCN
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

if exist('Fief7bee_reproducescript', 'var')
  cfg.reproducescript = Fief7bee_reproducescript;
end

% the output data should be saved to a MATLAB file
if (isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile)) || (isfield(cfg, 'reproducescript') && ~isempty(cfg.reproducescript))
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexlock(cfg.outputlock);
  end
  
  % WARNING: the following code is shared between ft_preamble_savefig and ft_preamble_savevar
  if isfield(cfg, 'reproducescript') && ~isempty(cfg.reproducescript)
    iW1aenge_now = datestr(now, 30);
    cfg.outputfile = sprintf('%s_output.fig', iW1aenge_now);
    
    % write the script that reproduces the analysis
    iW1aenge_cfg = removefields(cfg.callinfo.usercfg, ignorefields('reproducescript'));
    iW1aenge_cfg = copyfields(cfg, iW1aenge_cfg, {'outputfile'});
    iW1aenge_cfg = printstruct('cfg', iW1aenge_cfg);
    iW1aenge_fid = fopen(cfg.reproducescript, 'a+');
    fprintf(iW1aenge_fid, "%%%%\n\n");
    fprintf(iW1aenge_fid, "cfg = [];\n");
    fprintf(iW1aenge_fid, "%s\n", iW1aenge_cfg);
    iW1aenge_st = dbstack(2);
    fprintf(iW1aenge_fid, '%s(cfg);\n\n', iW1aenge_st(1).name);
    fclose(iW1aenge_fid);
  end
  
  if isequal(iW1aenge_postamble, {'varargout'}) && ~iscell(cfg.outputfile)
    % this should be a cell-array, oterwise it cannot be matched with varargout
    cfg.outputfile = {cfg.outputfile};
  end
  
  % save the current figure to a MATLAB .fig file
  ft_info('writing output figure to file ''%s''\n', cfg.outputfile);
  savefig(gcf, cfg.outputfile, 'compact')
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexunlock(cfg.outputlock);
  end
  
end
