function [sourceout] = ft_source2sparse(sourcein)

% FT_SOURCE2SPARSE removes the grid locations outside the brain from the source
% reconstruction, thereby saving memory.
%
% This invalidates the fields that describe the grid, and also makes it
% more difficult to make a plot of each of the slices of the source volume.
% The original source structure can be recreated using FT_SOURCE2FULL.
%
% Use as
%   [source] = ft_source2sparse(source)
%
% See also FT_SOURCE2FULL

% Copyright (C) 2004-2021, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% this function does not take a cfg as input
cfg = [];

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance sourcein

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

%%

sourcein = ft_checkdata(sourcein, 'datatype', 'source', 'insidestyle', 'logical');

inside  = find( sourcein.inside(:));
outside = find(~sourcein.inside(:));

fprintf('total number of dipoles        : %d\n', length(inside)+length(outside));
fprintf('number of dipoles inside  brain: %d\n', length(inside));
fprintf('number of dipoles outside brain: %d\n', length(outside));

%% construct the output data

% loop over parameters
fn = fieldnames(sourcein);
fn = setdiff(fn, {'dim', 'transform'});

% keep the descriptive fields that remain valid
sourceout = keepfields(sourcein, {'dim', 'transform', 'time', 'freq'});

for i=1:numel(fn)
  dimord = getdimord(sourcein, fn{i});
  
  clear tmp
  if startsWith(dimord, 'pos_pos') % this should go first
    tmp = sourcein.(fn{i});
    tmp = tmp(inside,inside,:,:,:,:,:);
    sourceout.(fn{i}) = tmp;
  elseif startsWith(dimord, 'pos')
    tmp = sourcein.(fn{i});
    tmp = tmp(inside,:,:,:,:,:,:);
    sourceout.(fn{i}) = tmp;
  elseif startsWith(dimord, '{pos}')
    tmp = sourcein.(fn{i});
    tmp = tmp(inside);
    sourceout.(fn{i}) = tmp;
  elseif startsWith(dimord, '{pos_pos}')
    tmp = sourcein.(fn{i});
    tmp = tmp(inside,inside);
    sourceout.(fn{i}) = tmp;
  else
    % skipping parameters with unknown dimord
  end
  
  if exist('tmp', 'var')
    ft_info('making "%s" with dimord "%s" sparse', fn{i}, dimord);
  else
    ft_info('skipping "%s" with dimord "%s"', fn{i}, dimord);
  end
  
end % for each field

%% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   sourcein
ft_postamble provenance sourceout
ft_postamble history    sourceout
