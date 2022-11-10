function [conjunction] = ft_conjunctionanalysis(cfg, varargin)

% FT_CONJUNCTIONANALYSIS finds the minimum statistic common across two or
% more contrasts, i.e. data following ft_xxxstatistics. Furthermore, it
% finds the overlap of sensors/voxels that show statistically significant
% results (a logical AND on the mask fields).
%
% Alternatively, it finds minimalistic mean power values in the
% input datasets. Here, a type 'relative change' baselinecorrection
% prior to conjunction is advised.
%
% Use as
%   [stat] = ft_conjunctionanalysis(cfg, stat1, stat2, .., statN)
%
% where the input data is the result from either FT_TIMELOCKSTATISTICS,
% FT_FREQSTATISTICS, or FT_SOURCESTATISTICS
%
% No configuration options are yet implemented.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS

% Copyright (C) 2010-2014, Arjen Stolk
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar varargin
ft_preamble provenance varargin

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% input check
ndatasets = length(varargin);
if ndatasets<2
  ft_error('not enough input arguments; there should be at least two');
end
% check if the input data is valid for this function
for i = 1:ndatasets
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'freq', 'source'}, 'feedback', 'yes');
end
fprintf('performing conjunction analysis on %d input datasets \n', ndatasets);
conjunction = [];

% determine datatype
isfreq     = ft_datatype(varargin{1}, 'freq');
istimelock = ft_datatype(varargin{1}, 'timelock');
issource   = ft_datatype(varargin{1}, 'source');

% conjunction loop, in case ndatasets > 2
for i = 1:ndatasets-1

  % align input arguments for conjunction
  if isempty(conjunction)
    data1 = varargin{i};
    data2 = varargin{i+1};
  else
    data1 = conjunction; % use already conjunct output
    data2 = varargin{i+1};
  end

  %% SOURCE DATA
  if issource

    if isfield(data1, 'stat')
      fprintf('minimum statistics on source level data \n');

      % equal size input check
      if ~isequal(size(data1.stat), size(data2.stat))
        ft_error('the input arguments have different sizes');
      end

      % prepare the output data structure
      conjunction = data1;

      if isfield(data1, 'posclusters') % remove cluster details
        fprintf('removing information about positive clusters\n');
        conjunction = rmfield(conjunction, 'posclusters');
        conjunction = rmfield(conjunction, 'posclusterslabelmat');
      end

      if isfield(data1, 'negclusters') % remove cluster details
        fprintf('removing information about negative clusters\n');
        conjunction = rmfield(conjunction, 'negclusters');
        conjunction = rmfield(conjunction, 'negclusterslabelmat');
      end

      fprintf('minimum statistics on stat fields \n');
      conjunction.stat = minimumstatistics(data1.stat, data2.stat);

      if isfield(data1, 'prob') && isfield(data2, 'prob') % conjunction on probabilities
        fprintf('minimum statistics on prob fields \n');
        conjunction.prob = maximumprobabilities(data1.prob, data2.prob);
      end

      if isfield(data1, 'mask') && isfield(data2, 'mask') % conjunction on mask parameters
        fprintf('logical AND on mask fields \n');
        conjunction.mask = logicalAND(data1.mask, data2.mask);
      end

    elseif isfield(data1, 'avg') && isfield(data2, 'avg') % conjunction on mean power values
      fprintf('minimum statistics on mean voxel power \n');

      % equal size input check
      if ~isequal(size(data1.avg.pow), size(data2.avg.pow))
        ft_error('the input arguments have different sizes');
      end

      conjunction = data1;
      conjunction.avg.pow = minimumstatistics(data1.avg.pow, data2.avg.pow);

    elseif isfield(data1, 'trial')
      fprintf('please first compute the averages with ft_sourcedescriptives \n');
    else
      fprintf('this source level data does not fit conjunction analysis \n');
    end
  end % end of source level conjunction

  %% SENSOR DATA
  if isfreq || istimelock

    if isfield(data1, 'stat') % conjunction on t-values
      fprintf('minimum statistics on sensor level data \n');

      % equal size input check
      if ~isequal(size(data1.stat), size(data2.stat))
        ft_error('the input arguments have different sizes');
      end

      % prepare the output data structure
      conjunction = data1;

      if isfield(data1, 'posclusters') % remove cluster details
        fprintf('removing information about positive clusters\n');
        conjunction = rmfield(conjunction, 'posclusters');
        conjunction = rmfield(conjunction, 'posclusterslabelmat');
      end

      if isfield(data1, 'negclusters') % remove cluster details
        fprintf('removing information about negative clusters\n');
        conjunction = rmfield(conjunction, 'negclusters');
        conjunction = rmfield(conjunction, 'negclusterslabelmat');
      end

      fprintf('minimum statistics on stat fields \n');
      conjunction.stat = minimumstatistics(data1.stat, data2.stat);

      if isfield(data1, 'prob') && isfield(data2, 'prob') % conjunction on probabilities
        fprintf('minimum statistics on prob fields \n');
        conjunction.prob = maximumprobabilities(data1.prob, data2.prob);
      end

      if isfield(data1, 'mask') && isfield(data2, 'mask') % conjunction on mask parameters
        fprintf('logical AND on mask fields \n');
        conjunction.mask = logicalAND(data1.mask, data2.mask);
      end

    elseif isfield(data1, 'powspctrm') && isfield(data2, 'powspctrm') % conjunction on mean power values
      fprintf('minimum statistics on mean sensor power \n');

      % equal size input check
      if ~isequal(size(data1.powspctrm), size(data2.powspctrm))
        ft_error('the input arguments have different sizes');
      end

      conjunction = data1;
      conjunction.powspctrm = minimumstatistics(data1.powspctrm, data2.powspctrm);

    elseif isfield(data1, 'avg') && isfield(data2, 'avg') % conjunction on mean signal amplitudes
      fprintf('minimum statistics on mean sensor amplitudes \n');

      % equal size input check
      if ~isequal(size(data1.avg), size(data2.avg))
        ft_error('the input arguments have different sizes');
      end

      conjunction = data1;
      conjunction.avg = minimumstatistics(data1.avg, data2.avg);

    elseif isfield(data1, 'trial')
      fprintf('please first compute the averages with ft_timelockdescriptives/ft_freqdescriptives \n');
    else
      fprintf('this sensor level data does not fit conjunction analysis \n');
    end
  end % end of sensor level conjunction

  clear data1; clear data2;
end % end of conjunction loop

%% UNIDENTIFIED DATA
if istimelock == 0 && isfreq == 0 && issource == 0
  fprintf('this data is not appropriate for conjunction analysis\n');
  conjunction = [];
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   varargin
ft_postamble provenance conjunction
ft_postamble history    conjunction
ft_postamble savevar    conjunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minstat] = minimumstatistics(variable1, variable2)
minAbsT   = min(abs(variable1), abs(variable2));  % minimum of the absolute values
equalSign = (sign(variable1) == sign(variable2)); % 1 is signs are equal, 0 otherwise
origSign  = sign(variable1);                      % sign(varagin2) gives same result
minstat   = minAbsT.*equalSign.*origSign;

function [maxprob] = maximumprobabilities(variable1, variable2)
maxprob = max(variable1, variable2); % maximum of the probabilities

function [logic] = logicalAND(variable1, variable2)
logic = (variable1 & variable2); % compute logical AND
