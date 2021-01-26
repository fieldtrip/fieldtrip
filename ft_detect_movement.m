function [cfg, movement] = ft_detect_movement(cfg, data)

% FT_SACCADE_DETECTION performs detection of movements such as saccades and
% microsaccades, but also joystick movements, from time series data over multiple
% trials. Different methods for detecting movements are implemented, which are
% described in detail below:
%
% VELOCITY2D - detects micro/saccades using a two-dimensional (2D) velocity according
% to "Engbert R, Kliegl R (2003) Vision Res 43:1035-1045". The vertical and the
% horizontal eyetracker time series (for one eye) are transformed into velocities and
% microsaccades are indentified as "outlier" eye movements that exceed a given
% threshold for velocity and duration. This method has the additional options
%     cfg.velocity2D.kernel   = vector 1 x nsamples, kernel to compute velocity (default = [1 1 0 -1 -1].*(data.fsample/6);
%     cfg.velocity2D.demean   = 'no' or 'yes', whether to apply centering correction (default = 'yes')
%     cfg.velocity2D.mindur   = minimum microsaccade durantion in samples (default = 3);
%     cfg.velocity2D.velthres = threshold for velocity outlier detection (default = 6);
%
% CLUSTERING - detects movements according to "Otero-Millan et al., (2014) J Vis 14".
%
% Use as
%   [cfg, movement] = ft_detect_movement(cfg, data)
% where the input data should be organised in a structure as obtained from the
% FT_PREPROCESSING function.
%
% The configuration can contain the following options
%   cfg.method  = string representing the method for movement detection
%                 'velocity2D' detects microsaccades using the 2D velocity
%                 'clustering' use unsupervised clustering method to detect microsaccades
%   cfg.channel = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details, (default = 'all')
%   cfg.trials  = 'all' or a selection given as a 1xN vector (default = 'all')
%
% The output argument "movement" is a Nx3 matrix. The first and second columns
% specify the begining and end samples of a movement period (saccade, joystick, ...),
% and the third column contains the peak velocity/acceleration movement. The thrid
% column allows to convert movements into spike data representation, making it
% compatible with the spike toolbox functions.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_DATABROWSER, FT_DATATYPE_SPIKE

% Copyright (C) 2014, Diego Lozano-Soldevilla
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

% FIXME the help mentioned the

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
data = ft_checkdata(data, 'datatype', {'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

if isfield(data, 'fsample')
  fsample = getsubfield(data, 'fsample');
else
  fsample = 1./(mean(diff(data.time{1})));
end

% set the defaults
cfg.method   = ft_getopt(cfg, 'method',   'velocity2D');
cfg.feedback = ft_getopt(cfg, 'feedback', 'yes');

% set the defaults for the various microsaccade detection methods
switch cfg.method
  case 'velocity2D'
    % Engbert R, Kliegl R (2003) Microsaccades uncover the orientation of
    % covert attention. Vision Res 43:1035-1045.
    kernel = [1 1 0 -1 -1].*(fsample/6); % this is equivalent to Engbert et al (2003) Vis Res, eqn. (1)
    if ~isfield(cfg.velocity2D, 'kernel'),   cfg.velocity2D.kernel  = kernel; end
    if ~isfield(cfg.velocity2D, 'demean'),   cfg.velocity2D.demean  = 'yes';  end
    if ~isfield(cfg.velocity2D, 'mindur'),   cfg.velocity2D.mindur  =  3;     end % minimum microsaccade duration in samples
    if ~isfield(cfg.velocity2D, 'velthres'), cfg.velocity2D.velthres = 6;     end
  case 'clustering'
    ft_error('not implemented yet');
    % Otero-Millan J, Castro JLA, Macknik SL, Martinez-Conde S (2014)
    % Unsupervised clustering method to detect microsaccades. J Vis 14.
  otherwise
    ft_error('unsupported option for cfg.method');
end

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'showcallinfo'});
data = ft_selectdata(tmpcfg, data);
[cfg, data] = rollback_provenance(cfg, data);

% determine the size of the data
ntrial = length(data.trial);
nchan  = length(data.label);   % number of channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movement = [];

ft_progress('init', cfg.feedback, 'processing trials');
% do all the computations
for i=1:ntrial
  ft_progress(i/ntrial, 'finding microsaccades trial %d of %d\n', i, ntrial);
  
  dat = data.trial{i};
  nsample = size(dat,2);
  
  switch cfg.method
    case 'velocity2D'
      
      % demean horizontal and vertical time courses
      if strcmp(cfg.velocity2D.demean, 'yes')
        dat = ft_preproc_polyremoval(dat, 0, 1, nsample);
      end
      
      %% eye velocity computation
      % deal with padding
      n = size(cfg.velocity2D.kernel,2);
      pad = ceil(n/2);
      dat = ft_preproc_padding(dat, 'localmean', pad);
      
      % convolution. See Engbert et al (2003) Vis Res, eqn. (1)
      if n<100
        % heuristic: for large kernel the convolution is faster when done along
        % the columns, weighing against the costs of doing the transposition.
        % the threshold of 100 is a bit ad hoc.
        vel = convn(dat,   cfg.velocity2D.kernel,   'same');
      else
        vel = convn(dat.', cfg.velocity2D.kernel.', 'same').';
      end
      % cut the eges
      vel = ft_preproc_padding(vel, 'remove', pad);
      
      %% microsaccade detection
      % compute velocity thresholds as in Engbert et al (2003) Vis Res, eqn. (2)
      medianstd = sqrt( median(vel.^2,2) - (median(vel,2)).^2 );
      
      % Engbert et al (2003) Vis Res, eqn. (3)
      radius = cfg.velocity2D.velthres*medianstd;
      
      % compute test criterion: ellipse equation
      test = sum((vel./radius(:,ones(1,nsample))).^2,1);
      sacsmp = find(test>1); % microsaccade's indexing
      
      %% determine microsaccades per trial
      % first find eye movements of n-consecutive time points
      j = find(diff(sacsmp)==1);
      j1 = [j; j+1];
      com = intersect(j,j+1);
      cut = ~ismember(j1,com);
      sacidx = reshape(j1(cut),2,[]);
      
      for k=1:size(sacidx,2)
        duration = sacidx(1,k):sacidx(2,k);
        if size(duration,2) >= cfg.velocity2D.mindur
          % finding peak velocity by Pitagoras
          begtrl = sacsmp(duration(1,1));
          endtrl = sacsmp(duration(1,end));
          
          [peakvel, smptrl] = max(sqrt(sum(vel(:,begtrl:endtrl).^2,1)));
          veltrl = sacsmp(duration(1,smptrl)); % peak velocity microsaccade sample -> important for spike conversion
          
          trlsmp = data.sampleinfo(i,1):data.sampleinfo(i,2);
          begsample = trlsmp(1, begtrl); % begining microsaccade sample
          endsample = trlsmp(1, endtrl); % end microsaccade sample
          velsample = trlsmp(1, veltrl); % velocity peak microsaccade sample
          movement(end+1,:) = [begsample endsample velsample];
        end
      end
      
    case 'clustering'
      % not implemented yet
  end
end
ft_progress('close');

ft_postamble trackconfig
ft_postamble provenance
ft_postamble debug
ft_postamble previous data
