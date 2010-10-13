function [Cval] = ft_conjunctionanalysis(cfg, varargin)

% FT_CONJUNCTIONANALYSIS finds minimalistic (maximum common) T values,
% clusters and probabilities in the input contrasts (n>=2).
% Furthermore, the output mask structure shows the overlap of sensors/voxels
% that are significantly different (logical AND).
%
% Alternatively, it finds minimalistic mean power values in the input datasets.
% Here, a type 'relative change' baselinecorrection prior to conjunction is advised.
%
% Use as:
%   [stat] = ft_conjunctionanalysis(cfg, data1, data2, .., dataN)
%
%   where data comes from
%   - ft_sourcestatistics post ft_sourceanalysis
%     with or without ft_sourceinterpolate
%   - ft_freqstatistics post ft_freqanalysis
%   - ft_sourceanalysis
%   - ft_freqanalysis
%
%   No configuration options are yet implemented.
%
% Copyright (C) 2010-2011, Arjen Stolk
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

% input check
if nargin < 3
    error('not enough input arguments; there should be at least two');
elseif nargin >= 3
    ndatasets = nargin - 1;
    fprintf('performing conjunction analysis on %s input datasets \n', num2str(ndatasets));
    Cval = [];
end

% output check
voxflag = 0;
sensflag = 0;

% conjunction loop, in case ndatasets > 2
for i = 1:ndatasets-1
    
    % align input arguments for conjunction
    if isempty(Cval)
        data1 = varargin{i};
        data2 = varargin{i+1};
    else
        data1 = Cval; % use already conjunct output
        data2 = varargin{i+1};
    end
    
    %% SOURCE DATA
    if isfield(data1, 'inside') % check
        voxflag = 1;
        
        if isfield(data1, 'stat') % conjunction on t-values
            fprintf('minimum statistics on voxel T values \n');
            
            % equal size input check
            if ~isequal(size(data1.stat), size(data2.stat))
                error('the input arguments have different sizes');
            end
            
            Cval = data1;
            Cval.stat = minimumstatistics(data1.stat, data2.stat);
            
            if isfield(data1, 'posclusterslabelmat') % conjunction on cluster values
                fprintf('minimum statistics on positive clusters \n');
                
                Cval.posclusterslabelmat = minimumstatistics(data1.posclusterslabelmat, data2.posclusterslabelmat);
            end
            
            if isfield(data1, 'negclusterslabelmat') % conjunction on cluster values
                fprintf('minimum statistics on negative clusters \n');
                
                Cval.negclusterslabelmat = minimumstatistics(data1.negclusterslabelmat, data2.negclusterslabelmat);
            end
            
            if isfield(data1, 'prob') % conjunction on probabilities
                fprintf('minimum statistics on probabilities \n');
                
                Cval.prob = maximumprobabilities(data1.prob, data2.prob);
            end
            
            if isfield(data1, 'mask') % conjunction on mask parameters
                fprintf('logical AND on masking parameters \n');
                
                Cval.mask = logicalAND(data1.mask, data2.mask);
            end
            
        elseif isfield(data1, 'avg') % conjunction on mean power values
            fprintf('minimum statistics on mean voxel power \n');
            
            % equal size input check
            if ~isequal(size(data1.avg.pow), size(data2.avg.pow))
                error('the input arguments have different sizes');
            end
            
            Cval = data1;
            Cval.avg.pow = minimumstatistics(data1.avg.pow, data2.avg.pow);
            
        elseif isfield(data1, 'trial')
            fprintf('please first compute the averages with ft_sourcedescriptives \n');
        else
            fprintf('this source level data does not fit conjunction analysis \n');
        end
    end % end of source level conjunction
    
    %% SENSOR DATA
    if isfield(data1, 'dimord') && ...
            (strcmp(data1.dimord, 'chan_freq') || ...
            strcmp(data1.dimord, 'chan_freq_time') || ...
            strcmp(data1.dimord, 'subj_chan_freq') || ...
            strcmp(data1.dimord, 'subj_chan_freq_time') || ...
            strcmp(data1.dimord, 'rpt_chan_freq') || ...
            strcmp(data1.dimord, 'rpt_chan_freq_time')); % check
        sensflag = 1;
        
        if isfield(data1, 'stat') % conjunction on t-values
            fprintf('minimum statistics on sensor T values \n');
            
            % equal size input check
            if ~isequal(size(data1.stat), size(data2.stat))
                error('the input arguments have different sizes');
            end
            
            Cval = data1;
            Cval.stat = minimumstatistics(data1.stat, data2.stat);
            
            if isfield(data1, 'posclusterslabelmat') % conjunction on cluster values
                fprintf('minimum statistics on positive clusters \n');
                
                Cval.posclusterslabelmat = minimumstatistics(data1.posclusterslabelmat, data2.posclusterslabelmat);
            end
            
            if isfield(data1, 'negclusterslabelmat') % conjunction on cluster values
                fprintf('minimum statistics on negative clusters \n');
                
                Cval.negclusterslabelmat = minimumstatistics(data1.negclusterslabelmat, data2.negclusterslabelmat);
            end
            
            if isfield(data1, 'prob') % conjunction on probabilities
                fprintf('minimum statistics on probabilities \n');
                
                Cval.prob = maximumprobabilities(data1.prob, data2.prob);
            end
            
            if isfield(data1, 'mask') % conjunction on mask parameters
                fprintf('logical AND on masking parameters \n');
                
                Cval.mask = logicalAND(data1.mask, data2.mask);
            end
            
        elseif isfield(data1, 'powspctrm') % conjunction on mean power values
            fprintf('minimum statistics on mean sensor power \n');
            
            % equal size input check
            if ~isequal(size(data1.powspctrm), size(data2.powspctrm))
                error('the input arguments have different sizes');
            end
            
            Cval = data1;
            Cval.powspctrm = minimumstatistics(data1.powspctrm, data2.powspctrm);
        else
            fprintf('this sensor level data does not fit conjunction analysis \n');
        end
    end % end of sensor level conjunction
    
    clear data1; clear data2;
end % end of conjunction loop

%% UNIDENTIFIED DATA
if voxflag == 0 && sensflag == 0
    fprintf('this data does not fit conjunction analysis \n');
    Cval = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minstat] = minimumstatistics(variable1, variable2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minAbsT   = min(abs(variable1), abs(variable2));  % minimum of the absolute values
equalSign = (sign(variable1) == sign(variable2)); % 1 is signs are equal, 0 otherwise
origSign  = sign(variable1);                      % sign(varagin2) gives same result
minstat   = minAbsT.*equalSign.*origSign;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxprob] = maximumprobabilities(variable1, variable2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxprob = max(variable1, variable2); % maximum of the probabilities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [logic] = logicalAND(variable1, variable2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logic = (variable1 & variable2); % compute logical AND