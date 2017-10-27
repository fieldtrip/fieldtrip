function [cfg movement] = ft_detect_movement(cfg, data)

% FT_SACCADE_DETECTION performs micro/saccade detection on time series data
% over multiple trials
%
% Use as
%   movement = ft_detect_movement(cfg, data)
%
% The input data should be organised in a structure as obtained from the
% FT_PREPROCESSING function. The configuration depends on the type of
% computation that you want to perform.
%
% The configuration should contain:
%  cfg.method   = different methods of detecting different movement types
%                'velocity2D', Micro/saccade detection based on Engbert R,
%                   Kliegl R (2003) Vision Res 43:1035-1045. The method
%                   computes thresholds based on velocity changes from
%                   eyetracker data (horizontal and vertical components).
%                'clustering', Micro/saccade detection based on
%                   Otero-Millan et al., (2014) J Vis 14 (not implemented
%                   yet)
%   cfg.channel = Nx1 cell-array with selection of channels, see
%                 FT_CHANNELSELECTION for details, (default = 'all')
%   cfg.trials  = 'all' or a selection given as a 1xN vector (default = 'all')
%
% METHOD SPECIFIC OPTIONS AND DESCRIPTIONS
%
%  VELOCITY2D
%   VELOCITY2D detects micro/saccades using a two-dimensional (2D) velocity
%   space velocity. The vertical and the horizontal eyetracker time series
%   (one eye) are transformed into velocities and microsaccades are
%   indentified as "outlier" eye movements that exceed a given velocity and
%   duration threshold.
%     cfg.velocity2D.kernel   = vector 1 x nsamples, kernel to compute velocity (default = [1 1 0 -1 -1].*(data.fsample/6);
%     cfg.velocity2D.demean   = 'no' or 'yes', whether to apply centering correction (default = 'yes')
%     cfg.velocity2D.mindur   = minimum microsaccade durantion in samples (default = 3);
%     cfg.velocity2D.velthres = threshold for velocity outlier detection (default = 6);
%
% The output argument "movement" is a Nx3 matrix. The first and second
% columns specify the begining and end samples of a movement period
% (saccade, joystic...), and the third column contains the peak
% velocity/acceleration movement. This last thrid column will allow to
% convert movements into spike data representation, making the spike
% toolbox functions compatible (not implemented yet).
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PLOT_MOVEMENT (not implemented yet)

% Copyright (C) 2014, Diego Lozano-Soldevilla, Robert Oostenveld
%
% $Id$

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

if isfield(data, 'fsample');
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
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
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
  time = data.time{i};

  ndatsample = size(dat,2);

  switch cfg.method
    case 'velocity2D'

      % demean horizontal and vertical time courses
      if strcmp(cfg.velocity2D.demean, 'yes');
        dat = ft_preproc_polyremoval(dat, 0, 1, ndatsample);
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
      test = sum((vel./radius(:,ones(1,ndatsample))).^2,1);
      sacsmp = find(test>1);% microsaccade's indexing

      %% determine microsaccades per trial
      % first find eye movements of n-consecutive time points
      j = find(diff(sacsmp)==1);
      j1 = [j; j+1];
      com = intersect(j,j+1);
      cut = ~ismember(j1,com);
      sacidx = reshape(j1(cut),2,[]);

      for k=1:size(sacidx,2);
        duration = sacidx(1,k):sacidx(2,k);
        if size(duration,2) >= cfg.velocity2D.mindur;
          % finding peak velocity by Pitagoras
          begtrl = sacsmp(duration(1,1));
          endtrl = sacsmp(duration(1,end));

          [peakvel smptrl] = max(sqrt(sum(vel(:,begtrl:endtrl).^2,1)));
          veltrl = sacsmp(duration(1,smptrl));% peak velocity microsaccade sample -> important for spike conversion

          trlsmp = data.sampleinfo(i,1):data.sampleinfo(i,2);
          begsample = trlsmp(1, begtrl); % begining microsaccade sample
          endsample = trlsmp(1, endtrl); % end microsaccade sample
          velsample = trlsmp(1, veltrl); % velocity peak microsaccade sample
          movement(end+1,:) = [begsample endsample velsample];
        end
      end

    case 'clustering';
      %not implemented yet
  end
end
ft_progress('close');

ft_postamble trackconfig
ft_postamble provenance
ft_postamble debug
ft_postamble previous data
