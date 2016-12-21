function [source] = ft_sourcedescriptives(cfg, source)

% FT_SOURCEDESCRIPTIVES computes descriptive parameters of the source
% analysis results.
%
% Use as
%   [source] = ft_sourcedescriptives(cfg, source)
%
% where cfg is a structure with the configuration details and source is the
% result from a beamformer source estimation. The configuration can contain
%   cfg.cohmethod        = 'regular', 'lambda1', 'canonical'
%   cfg.powmethod        = 'regular', 'lambda1', 'trace', 'none'
%   cfg.supmethod        = 'chan_dip', 'chan', 'dip', 'none' (default)
%   cfg.projectmom       = 'yes' or 'no' (default = 'no')
%   cfg.eta              = 'yes' or 'no' (default = 'no')
%   cfg.kurtosis         = 'yes' or 'no' (default = 'no')
%   cfg.keeptrials       = 'yes' or 'no' (default = 'no')
%   cfg.keepcsd          = 'yes' or 'no' (default = 'no')
%   cfg.keepnoisecsd     = 'yes' or 'no' (default = 'no')
%   cfg.keepmom          = 'yes' or 'no' (default = 'yes')
%   cfg.keepnoisemom     = 'yes' or 'no' (default = 'yes')
%   cfg.resolutionmatrix = 'yes' or 'no' (default = 'no')
%   cfg.feedback         = 'no', 'text' (default), 'textbar', 'gui'
%
% The following option only applies to LCMV single-trial timecourses.
%   cfg.fixedori         = 'within_trials' or 'over_trials' (default = 'over_trials')
%
% If repeated trials are present that have undergone some sort of
% resampling (i.e. jackknife, bootstrap, singletrial or rawtrial), the mean,
% variance and standard error of mean will be computed for all source
% parameters. This is done after applying the optional transformation
% on the power and projected noise.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_SOURCEANALYSIS, FT_SOURCESTATISTICS, FT_MATH

% Copyright (C) 2004-2015, Robert Oostenveld & Jan-Mathijs Schoffelen
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

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar source
ft_preamble provenance source
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
% source = ft_checkdata(source, 'datatype', 'source', 'feedback', 'yes');

cfg = ft_checkconfig(cfg, 'forbidden',   {'trials'});    % trial selection is not implented here, you may want to consider ft_selectdata

% DEPRECATED by roboos on 13 June 2013
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2199 for more details
% support for this functionality can be removed at the end of 2013
cfg = ft_checkconfig(cfg, 'deprecated',  {'transform'}); % please use ft_math instead

% set the defaults
cfg.transform        = ft_getopt(cfg, 'transform',        []);
cfg.projectmom       = ft_getopt(cfg, 'projectmom',       'no');% if yes -> svdfft
cfg.numcomp          = ft_getopt(cfg, 'numcomp',          1);
cfg.powmethod        = ft_getopt(cfg, 'powmethod',        []);% see below
cfg.cohmethod        = ft_getopt(cfg, 'cohmethod',        []);% see below
cfg.feedback         = ft_getopt(cfg, 'feedback',         'textbar');
cfg.supmethod        = ft_getopt(cfg, 'supmethod',        'none');
cfg.resolutionmatrix = ft_getopt(cfg, 'resolutionmatrix', 'no');
cfg.eta              = ft_getopt(cfg, 'eta',              'no');
cfg.fa               = ft_getopt(cfg, 'fa',               'no');
cfg.kurtosis         = ft_getopt(cfg, 'kurtosis',         'no');
cfg.keeptrials       = ft_getopt(cfg, 'keeptrials',       'no');
cfg.keepcsd          = ft_getopt(cfg, 'keepcsd',          'no');
cfg.keepmom          = ft_getopt(cfg, 'keepmom',          'yes');
cfg.keepnoisecsd     = ft_getopt(cfg, 'keepnoisecsd',     'no');
cfg.keepnoisemom     = ft_getopt(cfg, 'keepnoisemom',     'yes');
cfg.fwhm             = ft_getopt(cfg, 'fwhm',             'no');
cfg.fwhmremovecenter = ft_getopt(cfg, 'fwhmremovecenter', 0);
cfg.fwhmmethod       = ft_getopt(cfg, 'fwhmmethod',       'barnes');
cfg.fwhmmaxdist      = ft_getopt(cfg, 'fwhmmaxdist',      []);
cfg.fixedori         = ft_getopt(cfg, 'fixedori',         'over_trials');

% only works for minimumnormestimate
cfg.demean         = ft_getopt(cfg, 'demean',         'yes');
cfg.baselinewindow = ft_getopt(cfg, 'baselinewindow', [-inf 0]);
cfg.zscore         = ft_getopt(cfg, 'zscore',         'yes');

zscore = strcmp(cfg.zscore, 'yes');
demean = strcmp(cfg.demean, 'yes');

% get desired method from source structure
source.method = ft_getopt(source,'method',[]);

% this is required for backward compatibility with the old sourceanalysis
if isfield(source, 'method') && strcmp(source.method, 'randomized')
  source.method = 'randomization';
elseif isfield(source, 'method') && strcmp(source.method, 'permuted')
  source.method = 'permutation';
elseif isfield(source, 'method') && strcmp(source.method, 'jacknife')
  source.method = 'jackknife';
end

% determine the type of data, this is only relevant for a few specific types
ispccdata = isfield(source, 'avg')   && isfield(source.avg, 'csdlabel');
islcmvavg = isfield(source, 'avg')   && isfield(source, 'time') && isfield(source.avg,   'mom') && any(size(source.avg.pow)==1);
islcmvtrl = isfield(source, 'trial') && isfield(source, 'time') && isfield(source.trial, 'mom');
ismneavg  = isfield(source, 'avg')   && isfield(source, 'time') && isfield(source.avg,   'mom') && size(source.avg.pow, 2)==numel(source.time);

% check the consistency of the defaults
if strcmp(cfg.projectmom, 'yes')
  if isempty(cfg.powmethod)
    cfg.powmethod = 'regular'; % set the default
  elseif ~strcmp(cfg.powmethod, 'regular')
    error('unsupported powmethod in combination with projectmom');
  end
  if isempty(cfg.cohmethod)
    cfg.cohmethod = 'regular';% set the default
  elseif ~strcmp(cfg.cohmethod, 'regular')
    error('unsupported cohmethod in combination with projectmom');
  end
else
  if isempty(cfg.powmethod)
    cfg.powmethod = 'lambda1'; % set the default
  end
  if isempty(cfg.cohmethod)
    cfg.cohmethod = 'lambda1'; % set the default
  end
end

% this is required for backward compatibility with an old version of sourcedescriptives
if isfield(cfg, 'singletrial'), cfg.keeptrials = cfg.singletrial;  end

% do a validity check on the input data and specified options
if strcmp(cfg.resolutionmatrix, 'yes')
  if ~isfield(source.avg, 'filter')
    error('The computation of the resolution matrix requires keepfilter=''yes'' in sourceanalysis.');
  elseif ~isfield(source, 'leadfield')
    error('The computation of the resolution matrix requires keepleadfield=''yes'' in sourceanalysis.');
  end
end

if strcmp(cfg.fwhm, 'yes')
  if ~isfield(source.avg, 'filter')
    error('The computation of the fwhm requires keepfilter=''yes'' in sourceanalysis.');
  end
end

if strcmp(cfg.eta, 'yes') && strcmp(cfg.cohmethod, 'svdfft'),
  error('eta cannot be computed in combination with the application of svdfft');
end

if strcmp(cfg.keeptrials, 'yes') && ~strcmp(cfg.supmethod, 'none'),
  error('you cannot keep trials when you want to partialize something');
end

% set some flags for convenience
isnoise    = isfield(source, 'avg') && isfield(source.avg, 'noisecsd');
keeptrials = strcmp(cfg.keeptrials, 'yes');
projectmom = strcmp(cfg.projectmom, 'yes');

% determine the subfunction used for computing power
switch cfg.powmethod
  case 'regular'
    powmethodfun = @powmethod_regular;
  case 'lambda1'
    powmethodfun = @powmethod_lambda1;
  case 'trace'
    powmethodfun = @powmethod_trace;
  case 'none'
    powmethodfun = [];
  otherwise
    error('unsupported powmethod');
end

% represent the selection of sources in the brain as a row-vector with indices
insideindx = find(source.inside(:)');

if ispccdata
  % the source reconstruction was computed using the pcc beamformer

  Ndipole = size(source.pos,1);

  if ischar(source.avg.csdlabel{1}), source.avg.csdlabel = {source.avg.csdlabel}; end
  if numel(source.avg.csdlabel)==1,
    source.avg.csdlabel = repmat(source.avg.csdlabel, [Ndipole 1]);
  end

  dipsel     = find(strcmp(source.avg.csdlabel{1}, 'scandip'));
  refchansel = find(strcmp(source.avg.csdlabel{1}, 'refchan'));
  refdipsel  = find(strcmp(source.avg.csdlabel{1}, 'refdip'));
  supchansel = find(strcmp(source.avg.csdlabel{1}, 'supchan'));
  supdipsel  = find(strcmp(source.avg.csdlabel{1}, 'supdip'));

  % cannot handle reference channels and reference dipoles simultaneously
  if numel(refchansel)>0 && numel(refdipsel)>0
    error('cannot simultaneously handle reference channels and reference dipole');
  end

  % these are only used to count the number of reference/suppression dipoles and channels
  refsel = [refdipsel refchansel];
  supsel = [supdipsel supchansel];


  % first do the projection of the moment, if requested
  if projectmom
    source.avg.ori = cell(1, Ndipole);
    ft_progress('init', cfg.feedback, 'projecting dipole moment');
    for i=insideindx
      ft_progress(i/length(insideindx), 'projecting dipole moment %d/%d\n', i, length(insideindx));

      if numel(source.avg.csdlabel)>1,
        dipsel     = find(strcmp(source.avg.csdlabel{i}, 'scandip'));
        refchansel = find(strcmp(source.avg.csdlabel{i}, 'refchan'));
        refdipsel  = find(strcmp(source.avg.csdlabel{i}, 'refdip'));
        supchansel = find(strcmp(source.avg.csdlabel{i}, 'supchan'));
        supdipsel  = find(strcmp(source.avg.csdlabel{i}, 'supdip'));

        % these are only used to count the number of reference/suppression dipoles and channels
        refsel = [refdipsel refchansel];
        supsel = [supdipsel supchansel];
      end

      mom     = source.avg.mom{i}(dipsel,     :);
      ref     = source.avg.mom{i}(refdipsel,  :);
      sup     = source.avg.mom{i}(supdipsel,  :);
      refchan = source.avg.mom{i}(refchansel, :);
      supchan = source.avg.mom{i}(supchansel, :);
      % compute the projection of the scanning dipole along the direction of the dominant amplitude
      if length(dipsel)>1, [mom, rmom]  = svdfft(mom, cfg.numcomp, source.cumtapcnt); else rmom = []; end
      source.avg.ori{i} = rmom;
      % compute the projection of the reference dipole along the direction of the dominant amplitude
      if length(refdipsel)>1, [ref, rref] = svdfft(ref, 1, source.cumtapcnt); else rref = []; end
      % compute the projection of the supression dipole along the direction of the dominant amplitude
      if length(supdipsel)>1, [sup, rsup] = svdfft(sup, 1, source.cumtapcnt); else rsup = []; end

      % compute voxel-level fourier-matrix
      source.avg.mom{i} = cat(1, mom, ref, sup, refchan, supchan);

      % create rotation-matrix
      rotmat = zeros(0, length(source.avg.csdlabel{i}));
      if ~isempty(rmom),
        rotmat = [rotmat; rmom zeros(numel(refsel)+numel(supsel),1)];
      end
      if ~isempty(rref),
        rotmat = [rotmat; zeros(1, numel(dipsel)), rref, zeros(1,numel(refchansel)+numel(supsel))];
      end
      if ~isempty(rsup),
        rotmat = [rotmat; zeros(1, numel(dipsel)+numel(refdipsel)), rsup, zeros(1,numel(refchansel)+numel(supchansel))];
      end
      for j=1:length(supchansel)
        rotmat(end+1,:) = 0;
        rotmat(end,numel(dipsel)+numel(refdipsel)+numel(supdipsel)+j) = 1;
      end
      for j=1:length(refchansel)
        rotmat(end+1,:) = 0;
        rotmat(end,numel(dipsel)+numel(refdipsel)+numel(supdipsel)+numel(supchansel)+j) = 1;
      end

      % compute voxel-level csd-matrix
      if isfield(source.avg, 'csd'), source.avg.csd{i}           = rotmat * source.avg.csd{i} * rotmat'; end
      % compute voxel-level noisecsd-matrix
      if isfield(source.avg, 'noisecsd'), source.avg.noisecsd{i} = rotmat * source.avg.noisecsd{i} * rotmat'; end
      % compute rotated filter
      if isfield(source.avg, 'filter'),   source.avg.filter{i}   = rotmat * source.avg.filter{i}; end
      if isfield(source.avg, 'csdlabel'),
        % remember what the interpretation is of all CSD output components
        scandiplabel = repmat({'scandip'}, 1, cfg.numcomp);          % only one dipole orientation remains
        refdiplabel  = repmat({'refdip'},  1, length(refdipsel)>0);  % for svdfft at max. only one dipole orientation remains
        supdiplabel  = repmat({'supdip'},  1, length(supdipsel)>0);  % for svdfft at max. only one dipole orientation remains
        refchanlabel = repmat({'refchan'}, 1, length(refchansel));
        supchanlabel = repmat({'supchan'}, 1, length(supchansel));
        % concatenate all the labels
        source.avg.csdlabel{i} = cat(2, scandiplabel, refdiplabel, supdiplabel, refchanlabel, supchanlabel);
      end

      % compute rotated leadfield
      % FIXME in the presence of a refdip and/or supdip, this does not work; leadfield is Nx3
      if isfield(source,  'leadfield'),
        %FIXME this is a proposed dirty fix
        n1 = size(source.leadfield{i},2);
        %n2 = size(rotmat,2) - n1;
        n2 = size(rotmat,2) - n1 +1; %added 1 JM
        source.leadfield{i}    = source.leadfield{i} * rotmat(1:n2, 1:n1)';
      end
    end % for i=insideindx
    ft_progress('close');

    % update the indices
    dipsel     = find(strcmp(source.avg.csdlabel, 'scandip'));
    refchansel = find(strcmp(source.avg.csdlabel, 'refchan'));
    refdipsel  = find(strcmp(source.avg.csdlabel, 'refdip'));
    supchansel = find(strcmp(source.avg.csdlabel, 'supchan'));
    supdipsel  = find(strcmp(source.avg.csdlabel, 'supdip'));
    refsel     = [refdipsel refchansel];
    supsel     = [supdipsel supchansel];
  end % if projectmom

  if keeptrials
    cumtapcnt = source.cumtapcnt(:);
    sumtapcnt = cumsum([0;cumtapcnt]);
    Ntrial = length(cumtapcnt);

    ft_progress('init', cfg.feedback, 'computing singletrial voxel-level cross-spectral densities');
    for triallop = 1:Ntrial
      source.trial(triallop).csd = cell(Ndipole, 1);  % allocate memory for this trial
      source.trial(triallop).mom = cell(Ndipole, 1);  % allocate memory for this trial

      ft_progress(triallop/Ntrial, 'computing singletrial voxel-level cross-spectral densities %d%d\n', triallop, Ntrial);
      for i=insideindx
        dat = source.avg.mom{i};
        tmpmom = dat(:, sumtapcnt(triallop)+1:sumtapcnt(triallop+1));
        tmpcsd = (tmpmom * tmpmom') ./cumtapcnt(triallop);
        source.trial(triallop).mom{i} = tmpmom;
        source.trial(triallop).csd{i} = tmpcsd;
      end % for i=insideindx
    end % for triallop
    ft_progress('close');
    % remove the average, continue with separate trials, but keep track of
    % the csdlabel
    csdlabel = source.avg.csdlabel;
    source   = rmfield(source, 'avg');
  else
    fprintf('using average voxel-level cross-spectral densities\n');
    csdlabel = source.avg.csdlabel;
  end % if keeptrials

  % process the csdlabel for each of the dipoles
  hasrefdip  = true;
  hasrefchan = true;
  hassupdip  = true;
  hassupchan = true;

  dipselcell     = cell(Ndipole,1);
  refdipselcell  = cell(Ndipole,1);
  refchanselcell = cell(Ndipole,1);
  supdipselcell  = cell(Ndipole,1);
  supchanselcell = cell(Ndipole,1);

  for i = insideindx
    dipsel     = find(strcmp(csdlabel{i}, 'scandip'));
    refchansel = find(strcmp(csdlabel{i}, 'refchan'));
    refdipsel  = find(strcmp(csdlabel{i}, 'refdip'));
    supchansel = find(strcmp(csdlabel{i}, 'supchan'));
    supdipsel  = find(strcmp(csdlabel{i}, 'supdip'));

    hasrefdip  = ~isempty(refdipsel)  && hasrefdip; %NOTE: it has to be true for all dipoles!
    hasrefchan = ~isempty(refchansel) && hasrefchan;
    hassupdip  = ~isempty(supdipsel)  && hassupdip;
    hassupchan = ~isempty(supchansel) && hassupchan;

    dipselcell{i}     = dipsel;
    refdipselcell{i}  = refdipsel;
    refchanselcell{i} = refchansel;
    supdipselcell{i}  = supdipsel;
    supchanselcell{i} = supchansel;
  end

  if keeptrials
    % do the processing of the CSD matrices for each trial
    if ~strcmp(cfg.supmethod, 'none')
      error('suppression is only supported for average CSD');
    end
    %dipselcell = mat2cell(repmat(dipsel(:)', [Ndipole 1]), ones(Ndipole,1), length(dipsel));
    %if hasrefdip,  refdipselcell  = mat2cell(repmat(refdipsel(:)',  [Ndipole 1]), ones(Ndipole,1), length(refdipsel));  end
    %if hasrefchan, refchanselcell = mat2cell(repmat(refchansel(:)', [Ndipole 1]), ones(Ndipole,1), length(refchansel)); end
    %if hassupdip,  supdipselcell  = mat2cell(repmat(supdipsel(:)',  [Ndipole 1]), ones(Ndipole,1), length(supdipsel));  end
    %if hassupchan, supchanselcell = mat2cell(repmat(supchansel(:)', [Ndipole 1]), ones(Ndipole,1), length(supchansel)); end

    ft_progress('init', cfg.feedback, 'computing singletrial voxel-level power');
    for triallop = 1:Ntrial
      %initialize the variables
      source.trial(triallop).pow = zeros(Ndipole, 1);
      if hasrefdip,  source.trial(triallop).refdippow     = zeros(Ndipole, 1); end
      if hasrefchan, source.trial(triallop).refchanpow    = zeros(Ndipole, 1); end
      if hassupdip,  source.trial(triallop).supdippow     = zeros(Ndipole, 1); end
      if hassupchan, source.trial(triallop).supchanpow    = zeros(Ndipole, 1); end

      ft_progress(triallop/Ntrial, 'computing singletrial voxel-level power %d%d\n', triallop, Ntrial);
      source.trial(triallop).pow(source.inside) = cellfun(powmethodfun, source.trial(triallop).csd(source.inside), dipselcell(source.inside));
      if hasrefdip,  source.trial(triallop).refdippow(source.inside)  = cellfun(powmethodfun,source.trial(triallop).csd(source.inside), refdipselcell(source.inside));  end
      if hassupdip,  source.trial(triallop).supdippow(source.inside)  = cellfun(powmethodfun,source.trial(triallop).csd(source.inside), supdipselcell(source.inside));  end
      if hasrefchan, source.trial(triallop).refchanpow(source.inside) = cellfun(powmethodfun,source.trial(triallop).csd(source.inside), refchanselcell(source.inside)); end
      if hassupchan, source.trial(triallop).supchanpow(source.inside) = cellfun(powmethodfun,source.trial(triallop).csd(source.inside), supchanselcell(source.inside)); end
      %FIXME kan volgens mij niet
      if isnoise && isfield(source.trial(triallop), 'noisecsd'),
        % compute the power of the noise projected on each source component
        source.trial(triallop).noise = cellfun(powmethodfun,source.trial(triallop).csd, dipselcell);
        if hasrefdip,  source.trial(triallop).refdipnoise  = cellfun(powmethodfun,source.trial(triallop).noisecsd, refdipselcell);  end
        if hassupdip,  source.trial(triallop).supdipnoise  = cellfun(powmethodfun,source.trial(triallop).noisecsd, supdipselcell);  end
        if hasrefchan, source.trial(triallop).refchannoise = cellfun(powmethodfun,source.trial(triallop).noisecsd, refchanselcell); end
        if hassupchan, source.trial(triallop).supchannoise = cellfun(powmethodfun,source.trial(triallop).noisecsd, supchanselcell); end
      end % if isnoise
    end % for triallop
    ft_progress('close');

    if strcmp(cfg.keepcsd, 'no')
      source.trial = rmfield(source.trial, 'csd');
    end

  else
    % do the processing of the average CSD matrix
    for i=insideindx
      switch cfg.supmethod
        case 'chan_dip'
          supindx = [supdipsel supchansel];
          if i==insideindx(1), refsel  = refsel - length(supdipsel); end % adjust index only once
        case 'chan'
          supindx = supchansel;
        case 'dip'
          supindx = supdipsel;
          if i==insideindx(1), refsel  = refsel - length(supdipsel); end
        case 'none'
          % do nothing
          supindx = [];
      end
      tmpcsd  = source.avg.csd{i};
      scnindx = setdiff(1:size(tmpcsd,1), supindx);
      tmpcsd  = tmpcsd(scnindx, scnindx) - tmpcsd(scnindx, supindx)*pinv(tmpcsd(supindx, supindx))*tmpcsd(supindx, scnindx);
      source.avg.csd{i}   = tmpcsd;
    end % for i=insideindx
    %     source.avg.csdlabel = source.avg.csdlabel(scnindx);

    if isnoise && ~strcmp(cfg.supmethod, 'none')
      source.avg = rmfield(source.avg, 'noisecsd');
    end

    % initialize the variables
    source.avg.pow           = nan(Ndipole, 1);
    if hasrefdip,  source.avg.refdippow     = nan(Ndipole, 1); end
    if hasrefchan, source.avg.refchanpow    = nan(Ndipole, 1); end
    if hassupdip,  source.avg.supdippow     = nan(Ndipole, 1); end
    if hassupchan, source.avg.supchanpow    = nan(Ndipole, 1); end
    if isnoise
      source.avg.noise         = nan(Ndipole, 1);
      if hasrefdip,  source.avg.refdipnoise     = nan(Ndipole, 1); end
      if hasrefchan, source.avg.refchannoise    = nan(Ndipole, 1); end
      if hassupdip,  source.avg.supdipnoise     = nan(Ndipole, 1); end
      if hassupchan, source.avg.supchannoise    = nan(Ndipole, 1); end
    end % if isnoise
    if hasrefdip||hasrefchan, source.avg.coh    = nan(Ndipole, 1); end
    if strcmp(cfg.eta, 'yes'),
      source.avg.eta           = nan(Ndipole, 1);
      source.avg.ori             = cell(1, Ndipole);
    end
    if strcmp(cfg.eta, 'yes') && ~isempty(refsel),
      source.avg.etacsd = nan(Ndipole, 1);
      source.avg.ucsd   = cell(1, Ndipole);
    end
    if strcmp(cfg.fa, 'yes'),
      source.avg.fa = nan(Ndipole, 1);
    end

    for i=insideindx
      dipsel = dipselcell{i};
			refsel = [refchanselcell{i} refdipselcell{i}];

      % compute the power of each source component
      if strcmp(cfg.projectmom, 'yes') && cfg.numcomp>1,
        source.avg.pow(i) = powmethodfun(source.avg.csd{i}(dipselcell{i},dipselcell{i}), 1);
      else
        source.avg.pow(i) = powmethodfun(source.avg.csd{i}(dipselcell{i},dipselcell{i}));
      end

      if hasrefdip,  source.avg.refdippow(i)  = powmethodfun(source.avg.csd{i}(refdipsel,refdipsel));   end
      if hassupdip,  source.avg.supdippow(i)  = powmethodfun(source.avg.csd{i}(supdipsel,supdipsel));   end
      if hasrefchan, source.avg.refchanpow(i) = powmethodfun(source.avg.csd{i}(refchansel,refchansel)); end
      if hassupchan, source.avg.supchanpow(i) = powmethodfun(source.avg.csd{i}(supchansel,supchansel)); end
      if isnoise
        % compute the power of the noise projected on each source component
        if strcmp(cfg.projectmom, 'yes') && cfg.numcomp>1,
          source.avg.noise(i) = powmethodfun(source.avg.noisecsd{i}(dipselcell{i},dipselcell{i}), 1);
        else
          source.avg.noise(i) = powmethodfun(source.avg.noisecsd{i}(dipselcell{i},dipselcell{i}));
        end
        if hasrefdip,  source.avg.refdipnoise(i)  = powmethodfun(source.avg.noisecsd{i}(refdipsel,refdipsel));   end
        if hassupdip,  source.avg.supdipnoise(i)  = powmethodfun(source.avg.noisecsd{i}(supdipsel,supdipsel));   end
        if hasrefchan, source.avg.refchannoise(i) = powmethodfun(source.avg.noisecsd{i}(refchansel,refchansel)); end
        if hassupchan, source.avg.supchannoise(i) = powmethodfun(source.avg.noisecsd{i}(supchansel,supchansel)); end
      end % if isnoise

      if ~isempty(refsel)
        % compute coherence
        csd = source.avg.csd{i};
        switch cfg.cohmethod
          case 'regular'
            % assume that all dipoles have been projected along the direction of maximum power
            Pd                = abs(csd(dipsel, dipsel));
            Pr                = abs(csd(refsel, refsel));
            Cdr               = csd(dipsel, refsel);
            source.avg.coh(i) = (Cdr.^2) ./ (Pd*Pr);
          case 'lambda1'
            %compute coherence on Joachim Gross' way
            Pd                = lambda1(csd(dipsel, dipsel));
            Pr                = lambda1(csd(refsel, refsel));
            Cdr               = lambda1(csd(dipsel, refsel));
            source.avg.coh(i) = abs(Cdr).^2 ./ (Pd*Pr);
          case 'canonical'
            [ccoh, c2, v1, v2] = cancorr(csd, dipsel, refsel);
            [cmax, indmax]     = max(ccoh);
            source.avg.coh(i)  = ccoh(indmax);
          otherwise
            error('unsupported cohmethod');
        end % cohmethod
      end

      % compute eta
      if strcmp(cfg.eta, 'yes')
        [source.avg.eta(i), source.avg.ori{i}] = csd2eta(source.avg.csd{i}(dipselcell{i},dipselcell{i}));
        if ~isempty(refsel),
          %FIXME this only makes sense when only a reference signal OR a dipole is selected
          [source.avg.etacsd(i), source.avg.ucsd{i}] = csd2eta(source.avg.csd{i}(dipsel,refsel));
        end
      end

      %compute fa
      if strcmp(cfg.fa, 'yes')
        source.avg.fa(i) = csd2fa(source.avg.csd{i}(dipsel,dipsel));
      end
    end % for diplop

    if strcmp(cfg.keepcsd, 'no')
      source.avg = rmfield(source.avg, 'csd');
    end
    if strcmp(cfg.keepnoisecsd, 'no') && isnoise
      source.avg = rmfield(source.avg, 'noisecsd');
    end

  end

elseif ismneavg
  %the source reconstruction was computed using the minimumnormestimate and contains an average timecourse
  if demean
    begsmp = nearest(source.time, cfg.baselinewindow(1));
    endsmp = nearest(source.time, cfg.baselinewindow(2));
    ft_progress('init', cfg.feedback, 'baseline correcting dipole moments');
    for diplop=1:length(insideindx)
      ft_progress(diplop/length(insideindx), 'baseline correcting dipole moments %d/%d\n', diplop, length(insideindx));
      mom = source.avg.mom{insideindx(diplop)};
      mom = ft_preproc_baselinecorrect(mom, begsmp, endsmp);
      source.avg.mom{insideindx(diplop)} = mom;
    end
    ft_progress('close');
  end

  if projectmom
    if isfield(source, 'tri')
      nrm = normals(source.pos, source.tri, 'vertex');
      source.avg.phi = zeros(size(source.pos,1),1);
    end
    ft_progress('init', cfg.feedback, 'projecting dipole moment');
    for diplop=1:length(insideindx)
      ft_progress(diplop/length(insideindx), 'projecting dipole moment %d/%d\n', diplop, length(insideindx));
      mom = source.avg.mom{insideindx(diplop)};
      [mom, rmom] = svdfft(mom, 1);
      source.avg.mom{insideindx(diplop)} = mom;
      source.avg.ori{insideindx(diplop)} = rmom;
    end
    if isfield(source, 'tri')
      for diplop=insideindx
        source.avg.phi(diplop) = source.avg.ori{diplop}*nrm(diplop,:)';
      end
    end
    if isfield(source.avg, 'noisecov')
      source.avg.noise = nan+zeros(size(source.pos,1),1);
      for diplop=insideindx
        rmom = source.avg.ori{diplop};
        source.avg.noise(diplop) = rmom*source.avg.noisecov{diplop}*rmom';
      end
    end
    ft_progress('close');
  end

  if zscore
    begsmp = nearest(source.time, cfg.baselinewindow(1));
    endsmp = nearest(source.time, cfg.baselinewindow(2));
    % zscore using baselinewindow for power
    ft_progress('init', cfg.feedback, 'computing power');
    %source.avg.absmom = source.avg.pow;
    for diplop=1:length(insideindx)
      ft_progress(diplop/length(insideindx), 'computing power %d/%d\n', diplop, length(insideindx));
      mom = source.avg.mom{insideindx(diplop)};
      mmom = mean(mom(:,begsmp:endsmp),2);
      smom = std(mom(:,begsmp:endsmp),[],2);
      pow  = sum(((mom-mmom(:,ones(size(mom,2),1)))./smom(:,ones(size(mom,2),1))).^2,1);
      source.avg.pow(insideindx(diplop),:) = pow;
      %source.avg.absmom(source.inside(diplop),:) = sum((mom-mmom)./smom,1);
    end
    ft_progress('close');

  else
    % just square for power
    ft_progress('init', cfg.feedback, 'computing power');
    %source.avg.absmom = source.avg.pow;
    for diplop=1:length(insideindx)
      ft_progress(diplop/length(insideindx), 'computing power %d/%d\n', diplop, length(insideindx));
      mom = source.avg.mom{insideindx(diplop)};
      pow = sum(mom.^2,1);
      source.avg.pow(insideindx(diplop),:) = pow;
      %source.avg.absmom(insideindx(diplop),:) = sum(mom,1);
    end
    ft_progress('close');

  end


  if strcmp(cfg.kurtosis, 'yes')
    fprintf('computing kurtosis based on dipole timecourse\n');
    source.avg.k2 = nan(size(source.pos,1),1);
    for diplop=1:length(insideindx)
      mom = source.avg.mom{insideindx(diplop)};
      if length(mom)~=prod(size(mom))
        error('kurtosis can only be computed for projected dipole moment');
      end
      source.avg.k2(insideindx(diplop)) = kurtosis(mom);
    end
  end


elseif islcmvavg
  % the source reconstruction was computed using the lcmv beamformer and contains an average timecourse

  if projectmom
    ft_progress('init', cfg.feedback, 'projecting dipole moment');
    for diplop=1:length(insideindx)
      ft_progress(diplop/length(insideindx), 'projecting dipole moment %d/%d\n', diplop, length(insideindx));
      mom = source.avg.mom{insideindx(diplop)};
      [mom, rmom] = svdfft(mom, 1);
      source.avg.mom{insideindx(diplop)} = mom;
      source.avg.ori{insideindx(diplop)} = rmom;
    end
    ft_progress('close');
  end

  if ~strcmp(cfg.powmethod, 'none')
    fprintf('recomputing power based on dipole timecourse\n')
    source.avg.pow = nan(size(source.pos,1),1);
    for diplop=1:length(insideindx)
      mom = source.avg.mom{insideindx(diplop)};
      cov = mom * mom';
      source.avg.pow(insideindx(diplop)) = powmethodfun(cov);
    end
  end

  if strcmp(cfg.kurtosis, 'yes')
    fprintf('computing kurtosis based on dipole timecourse\n');
    source.avg.k2 = nan(size(source.pos,1),1);
    for diplop=1:length(insideindx)
      mom = source.avg.mom{insideindx(diplop)};
      if length(mom)~=prod(size(mom))
        error('kurtosis can only be computed for projected dipole moment');
      end
      source.avg.k2(insideindx(diplop)) = kurtosis(mom);
    end
  end

elseif islcmvtrl
  % the source reconstruction was computed using the lcmv beamformer and contains a single-trial timecourse
  ntrial = length(source.trial);

  if projectmom && strcmp(cfg.fixedori, 'within_trials')
    % the dipole orientation is re-determined for each trial
    ft_progress('init', cfg.feedback, 'projecting dipole moment');
    for trllop=1:ntrial
      ft_progress(trllop/ntrial, 'projecting dipole moment %d/%d\n', trllop, ntrial);
      for diplop=1:length(insideindx)
        mom = source.trial(trllop).mom{insideindx(diplop)};
        [mom, rmom] = svdfft(mom, 1);
        source.trial(trllop).mom{insideindx(diplop)} = mom;
        source.trial(trllop).ori{insideindx(diplop)} = rmom;  % remember the orientation
      end
    end
    ft_progress('close');
  elseif projectmom && strcmp(cfg.fixedori, 'over_trials')
    ft_progress('init', cfg.feedback, 'projecting dipole moment');
    % compute average covariance over all trials
    for trllop=1:ntrial
      for diplop=1:length(insideindx)
        mom = source.trial(trllop).mom{insideindx(diplop)};
        if trllop==1
          cov{diplop} = mom*mom'./size(mom,2);
        else
          cov{diplop} = mom*mom'./size(mom,2) + cov{diplop};
        end
      end
    end
    % compute source orientation over all trials
    for diplop=1:length(insideindx)
      [dum, ori{diplop}] = svdfft(cov{diplop}, 1);
    end
    % project the data in each trial
    for trllop=1:ntrial
      ft_progress(trllop/ntrial, 'projecting dipole moment %d/%d\n', trllop, ntrial);
      for diplop=1:length(insideindx)
        mom = source.trial(trllop).mom{insideindx(diplop)};
        mom = ori{diplop}*mom;
        source.trial(trllop).mom{insideindx(diplop)} = mom;
        source.trial(trllop).ori{insideindx(diplop)} = ori{diplop};
      end
    end
    ft_progress('close');
  end

  if ~strcmp(cfg.powmethod, 'none')
    fprintf('recomputing power based on dipole timecourse\n')
    for trllop=1:ntrial
      for diplop=1:length(insideindx)
        mom = source.trial(trllop).mom{insideindx(diplop)};
        cov = mom * mom';
        source.trial(trllop).pow(insideindx(diplop)) = powmethodfun(cov);
      end
    end
  end

  if strcmp(cfg.kurtosis, 'yes')
    fprintf('computing kurtosis based on dipole timecourse\n');
    for trllop=1:ntrial
      source.trial(trllop).k2 = nan(size(source.pos,1),1);
      for diplop=1:length(insideindx)
        mom = source.trial(trllop).mom{insideindx(diplop)};
        if length(mom)~=numel(mom)
          error('kurtosis can only be computed for projected dipole moment');
        end
        source.trial(trllop).k2(insideindx(diplop)) = kurtosis(mom);
      end
    end
  end

end % dealing with pcc or lcmv input

if isfield(source, 'avg') && isfield(source.avg, 'pow') && isfield(source.avg, 'noise') && ~ismneavg
  % compute the neural activity index for the average
  source.avg.nai = source.avg.pow(:) ./ source.avg.noise(:);
end

if isfield(source, 'trial') && isfield(source.trial, 'pow') && isfield(source.trial, 'noise')
  % compute the neural activity index for the trials
  ntrials = length(source.trial);
  for trlop=1:ntrials
    source.trial(trlop).nai = source.trial(trlop).pow ./ source.trial(trlop).noise;
  end
end

if strcmp(source.method, 'randomization') || strcmp(source.method, 'permutation')
  % compute the neural activity index for the two randomized conditions
  source.avgA.nai = source.avgA.pow ./ source.avgA.noise;
  source.avgB.nai = source.avgB.pow ./ source.avgB.noise;
  for trlop=1:length(source.trialA)
    source.trialA(trlop).nai = source.trialA(trlop).pow ./ source.trialA(trlop).noise;
  end
  for trlop=1:length(source.trialB)
    source.trialB(trlop).nai = source.trialB(trlop).pow ./ source.trialB(trlop).noise;
  end
end

if ~isempty(cfg.transform)
  fprintf('applying %s transformation on the power and projected noise\n', cfg.transform);
  % apply the specified transformation on the power
  if isfield(source, 'avg'   ) && isfield(source.avg   , 'pow'), source.avg .pow = feval(cfg.transform, source.avg .pow); end
  if isfield(source, 'avgA'  ) && isfield(source.avgA  , 'pow'), source.avgA.pow = feval(cfg.transform, source.avgA.pow); end
  if isfield(source, 'avgB'  ) && isfield(source.avgB  , 'pow'), source.avgB.pow = feval(cfg.transform, source.avgB.pow); end
  if isfield(source, 'trial' ) && isfield(source.trial , 'pow'), for i=1:length(source.trial ), source.trial (i).pow = feval(cfg.transform, source.trial (i).pow); end; end
  if isfield(source, 'trialA') && isfield(source.trialA, 'pow'), for i=1:length(source.trialA), source.trialA(i).pow = feval(cfg.transform, source.trialA(i).pow); end; end
  if isfield(source, 'trialB') && isfield(source.trialB, 'pow'), for i=1:length(source.trialB), source.trialB(i).pow = feval(cfg.transform, source.trialB(i).pow); end; end
  % apply the specified transformation on the projected noise
  if isfield(source, 'avg'   ) && isfield(source.avg   , 'noise'), source.avg .noise = feval(cfg.transform, source.avg .noise); end
  if isfield(source, 'avgA'  ) && isfield(source.avgA  , 'noise'), source.avgA.noise = feval(cfg.transform, source.avgA.noise); end
  if isfield(source, 'avgB'  ) && isfield(source.avgB  , 'noise'), source.avgB.noise = feval(cfg.transform, source.avgB.noise); end
  if isfield(source, 'trial' ) && isfield(source.trial , 'noise'), for i=1:length(source.trial ), source.trial (i).noise = feval(cfg.transform, source.trial (i).noise); end; end
  if isfield(source, 'trialA') && isfield(source.trialA, 'noise'), for i=1:length(source.trialA), source.trialA(i).noise = feval(cfg.transform, source.trialA(i).noise); end; end
  if isfield(source, 'trialB') && isfield(source.trialB, 'noise'), for i=1:length(source.trialB), source.trialB(i).noise = feval(cfg.transform, source.trialB(i).noise); end; end
end

if strcmp(source.method, 'pseudovalue')
  % compute the pseudovalues for the beamformer output
  avg = source.trial(1);        % the first is the complete average
  Ntrials = length(source.trial)-1; % the remaining are the leave-one-out averages
  pseudoval = [];
  if isfield(source.trial, 'pow')
    allavg = getfield(avg, 'pow');
    for i=1:Ntrials
      thisavg = getfield(source.trial(i+1), 'pow');
      thisval = Ntrials*allavg - (Ntrials-1)*thisavg;
      pseudoval(i).pow = thisval;
    end
  end
  if isfield(source.trial, 'coh')
    allavg = getfield(avg, 'coh');
    for i=1:Ntrials
      thisavg = getfield(source.trial(i+1), 'coh');
      thisval = Ntrials*allavg - (Ntrials-1)*thisavg;
      pseudoval(i).coh = thisval;
    end
  end
  if isfield(source.trial, 'nai')
    allavg = getfield(avg, 'nai');
    for i=1:Ntrials
      thisavg = getfield(source.trial(i+1), 'nai');
      thisval = Ntrials*allavg - (Ntrials-1)*thisavg;
      pseudoval(i).nai = thisval;
    end
  end
  if isfield(source.trial, 'noise')
    allavg = getfield(avg, 'noise');
    for i=1:Ntrials
      thisavg = getfield(source.trial(i+1), 'noise');
      thisval = Ntrials*allavg - (Ntrials-1)*thisavg;
      pseudoval(i).noise = thisval;
    end
  end
  % store the pseudovalues instead of the original values
  source.trial = pseudoval;
end

if strcmp(source.method, 'jackknife') || strcmp(source.method, 'bootstrap') || strcmp(source.method, 'pseudovalue') || strcmp(source.method, 'singletrial') || strcmp(source.method, 'rawtrial')
  % compute descriptive statistics (mean, var, sem) for multiple trial data
  % compute these for as many source parameters as possible

  % for convenience copy the trials out of the source structure
  dip = source.trial;

  % determine the (original) number of trials in the data
  if strcmp(source.method, 'bootstrap') %VERANDERD ER ZAT GEEN .RESAMPLE IN SOURCE
    Ntrials = size(source.trial,2);% WAS size(source.resample, 2);
  else
    Ntrials = length(source.trial);
  end
  fprintf('original data contained %d trials\n', Ntrials);

  % allocate memory for all elements in the dipole structure
  sumdip = [];
  if isfield(dip(1), 'var'),   sumdip.var    = zeros(size(dip(1).var  )); sumdip.var(~source.inside) = nan; end
  if isfield(dip(1), 'pow'),   sumdip.pow    = zeros(size(dip(1).pow  )); sumdip.pow(~source.inside) = nan; end
  if isfield(dip(1), 'coh'),   sumdip.coh    = zeros(size(dip(1).coh  )); sumdip.coh(~source.inside) = nan; end
  if isfield(dip(1), 'rv'),    sumdip.rv     = zeros(size(dip(1).rv   )); sumdip.rv(~source.inside) = nan; end
  if isfield(dip(1), 'noise'), sumdip.noise  = zeros(size(dip(1).noise)); sumdip.noise(~source.inside) = nan; end
  if isfield(dip(1), 'nai'),   sumdip.nai    = zeros(size(dip(1).nai  )); sumdip.nai(~source.inside) = nan; end
  sqrdip = [];
  if isfield(dip(1), 'var'),   sqrdip.var    = zeros(size(dip(1).var  )); sqrdip.var(~source.inside) = nan; end
  if isfield(dip(1), 'pow'),   sqrdip.pow    = zeros(size(dip(1).pow  )); sqrdip.pow(~source.inside) = nan; end
  if isfield(dip(1), 'coh'),   sqrdip.coh    = zeros(size(dip(1).coh  )); sqrdip.coh(~source.inside) = nan; end
  if isfield(dip(1), 'rv'),    sqrdip.rv     = zeros(size(dip(1).rv   )); sqrdip.rv(~source.inside) = nan; end
  if isfield(dip(1), 'noise'), sqrdip.noise  = zeros(size(dip(1).noise)); sqrdip.noise(~source.inside) = nan; end
  if isfield(dip(1), 'nai'),   sqrdip.nai    = zeros(size(dip(1).nai  )); sqrdip.nai(~source.inside) = nan; end
  if isfield(dip(1), 'mom')
    sumdip.mom = cell(size(dip(1).mom));
    sqrdip.mom = cell(size(dip(1).mom));
    for i=1:length(dip(1).mom)
      sumdip.mom{i} = zeros(size(dip(1).mom{i}));
      sqrdip.mom{i} = zeros(size(dip(1).mom{i}));
    end
  end
  if isfield(dip(1), 'csd')
    sumdip.csd = cell(size(dip(1).csd));
    sqrdip.csd = cell(size(dip(1).csd));
    for i=1:length(dip(1).csd)
      sumdip.csd{i} = zeros(size(dip(1).csd{i}));
      sqrdip.csd{i} = zeros(size(dip(1).csd{i}));
    end
  end

  for trial=1:length(dip)
    % compute the sum of all values
    if isfield(dip(trial), 'var'),    sumdip.var   = sumdip.var    + dip(trial).var;    end
    if isfield(dip(trial), 'pow'),    sumdip.pow   = sumdip.pow    + dip(trial).pow;    end
    if isfield(dip(trial), 'coh'),    sumdip.coh   = sumdip.coh    + dip(trial).coh;    end
    if isfield(dip(trial), 'rv'),     sumdip.rv    = sumdip.rv     + dip(trial).rv;     end
    if isfield(dip(trial), 'noise'),  sumdip.noise = sumdip.noise  + dip(trial).noise;  end
    if isfield(dip(trial), 'nai'),    sumdip.nai   = sumdip.nai    + dip(trial).nai;    end
    % compute the sum of squared values
    if isfield(dip(trial), 'var'),    sqrdip.var    = sqrdip.var   + (dip(trial).var  ).^2; end
    if isfield(dip(trial), 'pow'),    sqrdip.pow    = sqrdip.pow   + (dip(trial).pow  ).^2; end
    if isfield(dip(trial), 'coh'),    sqrdip.coh    = sqrdip.coh   + (dip(trial).coh  ).^2; end
    if isfield(dip(trial), 'rv'),     sqrdip.rv     = sqrdip.rv    + (dip(trial).rv   ).^2; end
    if isfield(dip(trial), 'noise'),  sqrdip.noise  = sqrdip.noise + (dip(trial).noise).^2; end
    if isfield(dip(trial), 'nai'),    sqrdip.nai    = sqrdip.nai   + (dip(trial).nai  ).^2; end
    % do the same for the cell array with mom
    if isfield(dip(trial), 'mom')
      for i=1:length(dip(1).mom)
        sumdip.mom{i} = sumdip.mom{i} +  dip(trial).mom{i};
        sqrdip.mom{i} = sqrdip.mom{i} + (dip(trial).mom{i}).^2;
      end
    end
    % do the same for the cell array with csd
    if isfield(dip(trial), 'csd')
      for i=1:length(dip(1).csd)
        sumdip.csd{i} = sumdip.csd{i} +  dip(trial).csd{i};
        sqrdip.csd{i} = sqrdip.csd{i} + (dip(trial).csd{i}).^2;
      end
    end
  end

  % compute the mean over all repetitions
  if isfield(sumdip, 'var'),    dipmean.var    = sumdip.var   / length(dip); end
  if isfield(sumdip, 'pow'),    dipmean.pow    = sumdip.pow   / length(dip); end
  if isfield(sumdip, 'coh'),    dipmean.coh    = sumdip.coh   / length(dip); end
  if isfield(sumdip, 'rv'),     dipmean.rv     = sumdip.rv    / length(dip); end
  if isfield(sumdip, 'noise'),  dipmean.noise  = sumdip.noise / length(dip); end
  if isfield(sumdip, 'nai'),    dipmean.nai    = sumdip.nai   / length(dip); end
  % for the cell array with mom, this is done further below
  % for the cell array with csd, this is done further below

  % the estimates for variance and SEM are biased if we are working with the jackknife/bootstrap
  % determine the proper variance scaling that corrects for this bias
  % note that Ntrials is not always the same as the length of dip, especially in case of the bootstrap
  if strcmp(source.method, 'singletrial')
    bias = 1;
  elseif strcmp(source.method, 'rawtrial')
    bias = 1;
  elseif strcmp(source.method, 'jackknife')
    % Effron gives SEM estimate for the jackknife method in equation 11.5 (paragraph 11.2)
    % to get the variance instead of SEM, we also have to multiply with the number of trials
    bias = (Ntrials - 1)^2;
  elseif strcmp(source.method, 'bootstrap')
    % Effron gives SEM estimate for the bootstrap method in algorithm 6.1 (equation 6.6)
    % to get the variance instead of SEM, we also have to multiply with the number of trials
    bias = Ntrials;
  elseif strcmp(source.method, 'pseudovalue')
    % note that I have not put any thought in this aspect yet
    warning('don''t know how to compute bias for pseudovalue resampling');
    bias = 1;
  end

  % compute the variance over all repetitions
  if isfield(sumdip, 'var'),    dipvar.var    = bias*(sqrdip.var    - (sumdip.var   .^2)/length(dip))/(length(dip)-1); end
  if isfield(sumdip, 'pow'),    dipvar.pow    = bias*(sqrdip.pow    - (sumdip.pow   .^2)/length(dip))/(length(dip)-1); end
  if isfield(sumdip, 'coh'),    dipvar.coh    = bias*(sqrdip.coh    - (sumdip.coh   .^2)/length(dip))/(length(dip)-1); end
  if isfield(sumdip, 'rv' ),    dipvar.rv     = bias*(sqrdip.rv     - (sumdip.rv    .^2)/length(dip))/(length(dip)-1); end
  if isfield(sumdip, 'noise' ), dipvar.noise  = bias*(sqrdip.noise  - (sumdip.noise .^2)/length(dip))/(length(dip)-1); end
  if isfield(sumdip, 'nai' ),   dipvar.nai    = bias*(sqrdip.nai    - (sumdip.nai   .^2)/length(dip))/(length(dip)-1); end

  % compute the SEM over all repetitions
  if isfield(sumdip, 'var'),    dipsem.var    = (dipvar.var   /Ntrials).^0.5; end
  if isfield(sumdip, 'pow'),    dipsem.pow    = (dipvar.pow   /Ntrials).^0.5; end
  if isfield(sumdip, 'coh'),    dipsem.coh    = (dipvar.coh   /Ntrials).^0.5; end
  if isfield(sumdip, 'rv' ),    dipsem.rv     = (dipvar.rv    /Ntrials).^0.5; end
  if isfield(sumdip, 'noise' ), dipsem.noise  = (dipvar.noise /Ntrials).^0.5; end
  if isfield(sumdip, 'nai' ),   dipsem.nai    = (dipvar.nai   /Ntrials).^0.5; end

  % compute the mean and SEM over all repetitions for the cell array with mom
  if isfield(dip(trial), 'mom')
    for i=1:length(dip(1).mom)
      dipmean.mom{i} = sumdip.mom{i}/length(dip);
      dipvar.mom{i} = bias*(sqrdip.mom{i} - (sumdip.mom{i}.^2)/length(dip))/(length(dip)-1);
      dipsem.mom{i} = (dipvar.mom{i}/Ntrials).^0.5;
    end
  end

  % compute the mean and SEM over all repetitions for the cell array with csd
  if isfield(dip(trial), 'csd')
    for i=1:length(dip(1).csd)
      dipmean.csd{i} = sumdip.csd{i}/length(dip);
      dipvar.csd{i} = bias*(sqrdip.csd{i} - (sumdip.csd{i}.^2)/length(dip))/(length(dip)-1);
      dipsem.csd{i} = (dipvar.csd{i}/Ntrials).^0.5;
    end
  end

  if strcmp(source.method, 'pseudovalue')
    % keep the trials, since they have been converted to pseudovalues
    % and hence the trials contain the interesting data
  elseif keeptrials
    % keep the trials upon request
  else
    % remove the original trials
    source = rmfield(source, 'trial');
    % assign the descriptive statistics to the output source structure
    source.avg = dipmean;
    source.var = dipvar;
    source.sem = dipsem;
  end
end

if strcmp(cfg.resolutionmatrix, 'yes')
  % this is only implemented for pcc and no refdips/chans at the moment
  Nchan        = size(source.leadfield{insideindx(1)}, 1);
  Ninside      = length(insideindx);
  allfilter    = zeros(Ninside,Nchan);
  allleadfield = zeros(Nchan,Ninside);
  dipsel       = match_str(source.avg.csdlabel, 'scandip');
  ft_progress('init', cfg.feedback, 'computing resolution matrix');
  for diplop=1:length(insideindx)
    ft_progress(diplop/length(insideindx), 'computing resolution matrix %d/%d\n', diplop, length(insideindx));
    % concatenate all filters
    allfilter(diplop,:)    = source.avg.filter{insideindx(diplop)}(dipsel,:);
    % concatenate all leadfields
    allleadfield(:,diplop) = source.leadfield{insideindx(diplop)};
  end
  ft_progress('close');
  % multiply the filters and leadfields to obtain the resolution matrix
  % see equation 1 and 2 in De Peralta-Menendez RG, Gonzalez-Andino SL: A critical analysis of linear inverse solutions to the neuroelectromagnetic inverse problem. IEEE Transactions on Biomedical Engineering 45: 440-448, 1998.
  source.resolution = nan(Ndipole, Ndipole);
  source.resolution(insideindx, insideindx) = allfilter*allleadfield;
end

% compute fwhm
if strcmp(cfg.fwhm, 'yes')
  switch cfg.fwhmmethod
    case 'barnes'
      if ~isfield(source, 'dim')
        error('computation of fwhm is not possible with method ''barnes'' is not possible when the dipoles are not defined on a regular 3D grid');
      end
      fprintf('computing fwhm of spatial filters using method ''barnes''\n');
      source = estimate_fwhm1(source, cfg.fwhmremovecenter);
    case 'gaussfit'
      fprintf('computing fwhm of spatial filters using method ''gaussfit''\n');
      source = estimate_fwhm2(source, cfg.fwhmmaxdist);
    otherwise
      error('unknown method for fwhm estimation');
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   source
ft_postamble provenance source
ft_postamble history    source
ft_postamble savevar    source


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute eta from a csd-matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta, u] = csd2eta(csd)
[u,s,v] = svd(real(csd));
eta     = s(2,2)./s(1,1);
u       = u'; %orientation is defined in the rows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute fa from a csd-matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fa] = csd2fa(csd)
s  = svd(real(csd));
ns = rank(real(csd));
s  = s(1:ns);
ms = mean(s);
fa = sqrt( (ns./(ns-1)) .* (sum((s-ms).^2))./(sum(s.^2)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = powmethod_lambda1(x, ind)

if nargin==1,
  ind = 1:size(x,1);
end
s = svd(x(ind,ind));
p = s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = powmethod_trace(x, ind)

if nargin==1,
  ind = 1:size(x,1);
end
p = trace(x(ind,ind));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = powmethod_regular(x, ind)

if nargin==1,
  ind = 1:size(x,1);
end
p = abs(x(ind,ind));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to obtain the largest singular value or trace of the
% source CSD matrices resulting from DICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = lambda1(x)
s = svd(x);
s = s(1);
