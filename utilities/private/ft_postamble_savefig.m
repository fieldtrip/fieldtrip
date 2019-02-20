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

% the output data should be saved to a MATLAB file
if (isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile)) || exist('Fief7bee_reproducescript', 'var')
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexlock(cfg.outputlock);
  end
  
  if exist('Fief7bee_reproducescript', 'var')
    % write the output figure to a MATLAB file
    iW1aenge_now = datestr(now, 30);
    cfg.outputfile = fullfile(Fief7bee_reproducescript, sprintf('%s_output', iW1aenge_now));
    % write a snippet of MATLAB code with the configuration and function call
    reproducescript(fullfile(Fief7bee_reproducescript, 'script.m'), cfg, false)
  else
    cfg.outputfile = [];
  end
  
  if ~isempty(cfg.outputfile)
    % save the current figure to a MATLAB .fig and to a .png file
    [Fief7bee_p, Fief7bee_f, Fief7bee_x] = fileparts(cfg.outputfile);
    Fief7bee_outputfile = fullfile(Fief7bee_p, [Fief7bee_f, '.fig']);
    ft_info('writing output figure to file ''%s''\n', Fief7bee_outputfile);
    savefig(gcf, cfg.outputfile, 'compact')
    Fief7bee_outputfile = fullfile(Fief7bee_p, [Fief7bee_f, '.png']);
    ft_info('writing screenshot to file ''%s''\n', Fief7bee_outputfile);
    set(gcf, 'PaperOrientation', 'portrait');
    print(Fief7bee_outputfile, '-dpng');
  end
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexunlock(cfg.outputlock);
  end
  
end
