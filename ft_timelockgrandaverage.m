function [grandavg] = ft_timelockgrandaverage(cfg, varargin)

% FT_TIMELOCKGRANDAVERAGE computes ERF/ERP average and variance
% over multiple subjects
%
% Use as
%   [grandavg] = ft_timelockgrandaverage(cfg, avg1, avg2, avg3, ...)
%
% where
%   avg1..N are the ERF/ERP averages as obtained from FT_TIMELOCKANALYSIS
%
% and cfg is a configuration structure with
%  cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                       see FT_CHANNELSELECTION for details
%  cfg.latency        = [begin end] in seconds or 'all' (default = 'all')
%  cfg.keepindividual = 'yes' or 'no' (default = 'no')
%  cfg.normalizevar   = 'N' or 'N-1' (default = 'N-1')
%  cfg.parameter      = string or cell-array of strings indicating which
%                        parameter(s) to average. default is set to
%                        'avg', if it is present in the data.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. For this particular function, the input should be
% structured as a cell array.
%
% See also FT_TIMELOCKANALYSIS, FT_TIMELOCKSTATISTICS

% Copyright (C) 2003-2006, Jens Schwarzbach
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble loadvar varargin

% return immediately after distributed execution
if ~isempty(ft_getopt(cfg, 'distribute'))
    return
end

% check if the input data is valid for this function
for i=1:length(varargin)
    varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
end

% set the defaults
cfg.inputfile      = ft_getopt(cfg, 'inputfile',  []);
cfg.outputfile     = ft_getopt(cfg, 'outputfile', []);
cfg.keepindividual = ft_getopt(cfg, 'keepindividual', 'no');
cfg.channel        = ft_getopt(cfg, 'channel',    'all');
cfg.latency        = ft_getopt(cfg, 'latency',    'all');
cfg.normalizevar   = ft_getopt(cfg, 'normalizevar', 'N-1');
cfg.parameter      = ft_getopt(cfg, 'parameter',  'avg');

if ischar(cfg.parameter)
    cfg.parameter = {cfg.parameter};
end

% FIXME: ft_timelockgrandaverage either has .avg or .individual for output,
% and therefore yet supports only 1 input parameter
if numel(cfg.parameter) > 1
    fprintf('supporting 1 parameter, choosing the first parameter: %s\n', cfg.parameter{1});
    temp = cfg.parameter{1}; % remove the other parameters
    cfg = rmfield(cfg, 'parameter');
    cfg.parameter{1} = temp;
end

Nsubj    = length(varargin);
dimord   = varargin{1}.dimord;
hastime  = ~isempty(strfind(varargin{1}.dimord, 'time'));
hasrpt   = ~isempty(strfind(varargin{1}.dimord, 'rpt'));

% check whether the input data is suitable
if hasrpt
    fprintf('ignoring the single-subject repetition dimension\n');
    if strcmp(cfg.parameter, 'trial')
        error('not supporting averaging over the repetition dimension');
    end
end

if ischar(cfg.latency) && strcmp(cfg.latency, 'all')
    tbeg = min(varargin{1}.time);
    tend = max(varargin{1}.time);
else
    tbeg = cfg.latency(1);
    tend = cfg.latency(2);
end

% select the data in all inputs
for k=1:numel(cfg.parameter)
    
    % determine which channels, and latencies are available for all inputs
    for i=1:Nsubj
        cfg.channel = ft_channelselection(cfg.channel, varargin{i}.label);
        if hastime
            tbeg = max(tbeg, varargin{i}.time(1  ));
            tend = min(tend, varargin{i}.time(end));
        end
    end
    cfg.latency = [tbeg tend];
    
    % pick the selections
    for i=1:Nsubj
        if ~isfield(varargin{i}, cfg.parameter{k})
            error('the field %s is not present in data structure %d', cfg.parameter{k}, i);
        end
        chansel = match_str(varargin{i}.label, cfg.channel);
        varargin{i}.label = varargin{i}.label(chansel);
        if hastime
            timesel = nearest(varargin{i}.time, cfg.latency(1)):nearest(varargin{i}.time, cfg.latency(2));
            varargin{i}.time = varargin{i}.time(timesel);
        end
        % select the overlapping samples in the data matrix
        switch dimord
            case 'chan'
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(chansel);
            case 'chan_time'
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(chansel,timesel);
            case {'rpt_chan' 'subj_chan'}
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(chansel);
                varargin{i}.dimord = 'chan';
            case {'rpt_chan_time' 'subj_chan_time'}
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(chansel,timesel);
                varargin{i}.dimord = 'chan_time';
            otherwise
                error('unsupported dimord');
        end
    end % for i = subject
end % for k = parameter

% determine the size of the data to be averaged
dim = cell(1,numel(cfg.parameter));
for k=1:numel(cfg.parameter)
    dim{k} = size(varargin{1}.(cfg.parameter{k}));
end

% give some feedback on the screen
if strcmp(cfg.keepindividual, 'no')
    for k=1:numel(cfg.parameter)
        fprintf('computing average %s over %d subjects\n', cfg.parameter{k}, Nsubj);
    end
else
    for k=1:numel(cfg.parameter)
        fprintf('not computing average, but keeping individual %s for %d subjects\n', cfg.parameter{k}, Nsubj);
    end
end

% allocate memory to hold the data and collect it
for k=1:numel(cfg.parameter)
    avgmat = zeros([Nsubj, dim{k}]);
    for s=1:Nsubj
        avgmat(s, :, :) = varargin{s}.(cfg.parameter{k});
    end
    if strcmp(cfg.keepindividual, 'no')
        % average across subject dimension
        ResultGrandavg = mean(avgmat, 1);
        grandavg.avg = reshape(ResultGrandavg, [dim{k}]);  % Nchan x Nsamples
        switch cfg.normalizevar
            case 'N-1'
                sdflag = 0;
            case 'N'
                sdflag = 1;
        end
        ResultVar = std(avgmat, sdflag, 1).^2;
        grandavg.var = reshape(ResultVar, [dim{k}]);  % Nchanof each subjeof each subject should be an average, use FT_TIMELOCKANALYSIS firsct should be an average, use FT_TIMELOCKANALYSIS firs x Nsamples
    elseif strcmp(cfg.keepindividual, 'yes')
        grandavg.individual = avgmat;         % Nsubj x Nchan x Nsamples
    end
    % grandavg.(cfg.parameter{k}) = tmp; % not yet supported; now
    % implemented as individual vs. avg field
end

% collect the output data
grandavg.label = varargin{1}.label;
if isfield(varargin{1}, 'time')
    % remember the time axis
    grandavg.time = varargin{1}.time;
end
if isfield(varargin{1}, 'labelcmb')
    grandavg.labelcmb = varargin{1}.labelcmb;
end
if isfield(varargin{1}, 'grad')
    warning('discarding gradiometer information because it cannot be averaged');
end
if isfield(varargin{1}, 'elec')
    warning('discarding electrode information because it cannot be averaged');
end
if strcmp(cfg.keepindividual, 'yes')
    grandavg.dimord = ['subj_',varargin{1}.dimord];
elseif strcmp(cfg.keepindividual, 'no')
    grandavg.dimord = varargin{1}.dimord;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
ft_postamble history grandavg
ft_postamble savevar grandavg
