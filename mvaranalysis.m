function [mvardata] = mvaranalysis(cfg, data)

% MVARANALYSIS performs multivariate autoregressive modeling on any
% time series trial data.
%
% Use as
%   [mvardata] = mvaranalysis(cfg, data)
%
%%% The bivariate quantities should be interpreted according to the
%%% in which the labels occur in labelcmb. Here the second column
%%% causes the first column

% Copyright (C) 2009, Jan-Mathijs Schoffelen
% $Log: mvaranalysis.m,v $
% Revision 1.1  2009/10/02 13:55:23  jansch
% first implementation in fieldtrip
%

if ~isfield(cfg, 'channel'),    cfg.channel    = 'all';          end
if ~isfield(cfg, 'channelcmb'), cfg.channelcmb = {'all' 'all'};  end
if ~isfield(cfg, 'toi'),        cfg.toi        = 'all';          end
if ~isfield(cfg, 't_ftimwin'),  cfg.t_ftimwin  = [];             end
if ~isfield(cfg, 'order'),      cfg.order      = 10;             end
if ~isfield(cfg, 'keeptrials'), cfg.keeptrials = 'no';           end
if ~isfield(cfg, 'jackknife'),  cfg.jackknife  = 'no';           end
if ~isfield(cfg, 'keeptapers'), cfg.keeptapers = 'yes';          end
if ~isfield(cfg, 'taper'),      cfg.taper      = 'rectwin';      end
if ~isfield(cfg, 'mvarmethod'), cfg.mvarmethod = 2;              end %is Mode==2 in biosig-toolbox
if ~isfield(cfg, 'zscore'),     cfg.zscore     = 'no';           end
if ~isfield(cfg, 'blc'),        cfg.blc        = 'yes';          end
if ~isfield(cfg, 'feedback'),   cfg.feedback   = 'textbar';      end

if strcmp(cfg.toi,    'all'),
  %FIXME do something?
end  

cfg.channel    = channelselection(cfg.channel,      data.label);
%cfg.channelcmb = channelcombination(cfg.channelcmb, data.label);

keeprpt  = strcmp(cfg.keeptrials, 'yes');
keeptap  = strcmp(cfg.keeptapers, 'yes');
dojack   = strcmp(cfg.jackknife,  'yes');
dozscore = strcmp(cfg.zscore,     'yes');

if ~keeptap, error('not keeping tapers is not possible yet'); end
if dojack && keeprpt, error('you cannot simultaneously keep trials and do jackknifing'); end

tfwin    = round(data.fsample.*cfg.t_ftimwin);
ntrl     = length(data.trial);
ntoi     = length(cfg.toi);
chanindx = match_str(data.label, cfg.channel);
nchan    = length(chanindx);
label    = data.label(chanindx);

ncmb     = nchan*nchan;
cmbindx1 = repmat(chanindx(:),  [1 nchan]);
cmbindx2 = repmat(chanindx(:)', [nchan 1]);
labelcmb = [data.label(cmbindx1(:)) data.label(cmbindx2(:))];

if strcmp(cfg.taper, 'dpss')
  % create a sequence of DPSS (Slepian) tapers
  % ensure that the input arguments are double precision
  tap = double_dpss(tfwin,tfwin*(cfg.tapsmofrq./data.fsample))';
  tap = tap(1,:); %only use first 'zero-order' taper
elseif strcmp(cfg.taper, 'sine')
  tap = sine_taper(tfwin, tfwin*(cfg.tapsmofrq./data.fsample))';
  tap = tap(1,:);
else
  % create a single taper according to the window specification as a
  % replacement for the DPSS (Slepian) sequence
  tap = window(cfg.taper, tfwin)';
  tap = tap./norm(tap);
  % freqanalysis_mvar throws away the last taper of the Slepian sequence, add a dummy taper
end
ntap = size(tap,1);

%---allocate memory
coeffs   = zeros(nchan, nchan, cfg.order, ntoi, ntap);
noisecov = zeros(nchan, nchan,            ntoi, ntap);

%---preprocess data if necessary
%---blc
if strcmp(cfg.blc, 'yes'),
  tmpcfg           = [];
  tmpcfg.blc       = 'yes';
  tmpcfg.blcwindow = cfg.toi([1 end]) + cfg.t_ftimwin.*[-0.5 0.5];
  data             = preprocessing(tmpcfg, data);
else
  %do nothing
end

%---zscore
if dozscore,
  zwindow = cfg.toi([1 end]) + cfg.t_ftimwin.*[-0.5 0.5];
  sumval  = 0;
  sumsqr  = 0;
  numsmp  = 0;
  trlindx = [];
  for k = 1:ntrl
    begsmp = nearest(data.time{k}, zwindow(1));
    endsmp = nearest(data.time{k}, zwindow(2));
    if endsmp>=begsmp,
      sumval  = sumval + sum(data.trial{k}(:, begsmp:endsmp),    2);
      sumsqr  = sumsqr + sum(data.trial{k}(:, begsmp:endsmp).^2, 2);
      numsmp  = numsmp + endsmp - begsmp + 1;
      trlindx = [trlindx; k];
    end
  end
  datavg = sumval./numsmp;
  datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);
  
  data.trial = data.trial(trlindx);
  data.time  = data.time(trlindx);
  ntrl       = length(trlindx);
  for k = 1:ntrl
    rvec          = ones(1,size(data.trial{k},2));
    data.trial{k} = (data.trial{k} - datavg*rvec)./(datstd*rvec);
  end
else
  %do nothing
end

%---generate time axis
maxtim = -inf;
mintim = inf;
for k = 1:ntrl
  maxtim = max(maxtim, data.time{k}(end));
  mintim = min(mintim, data.time{k}(1));
end
timeaxis = mintim:1./data.fsample:maxtim;

%---loop over the tois
progress('init', cfg.feedback, 'computing MAR-model');
for j = 1:ntoi
  progress(j/ntoi, 'processing timewindow %d from %d\n', j, ntoi);
  sample        = nearest(timeaxis, cfg.toi(j));
  cfg.toi(j)    = timeaxis(sample);
  
  tmpcfg        = [];
  tmpcfg.toilim = [timeaxis(sample-floor(tfwin/2)) timeaxis(sample+ceil(tfwin/2)-1)];
  tmpcfg.feedback = 'no';
  tmpcfg.minlength= 'maxperlen';
  tmpdata       = redefinetrial(tmpcfg, data);
  
  cfg.toi(j)    = mean(tmpdata.time{1}([1 end]))+0.5./data.fsample;
  
  for m = 1:ntap
    %---construct data-matrix
    dat  = zeros(nchan, 0);
    for k = 1:length(tmpdata.trial)
      tmpdat  = tmpdata.trial{k}(chanindx, :).*tap(m.*ones(nchan,1),:);
      tmptime = tmpdata.time{k};
      dat     = [dat tmpdat nan*ones(nchan, cfg.order)];
    end

    tmpdat = dat';
    
    %---estimate autoregressive model
    [ar,rc,pe] = mvar(tmpdat, cfg.order, cfg.mvarmethod);
    coeffs(:,:,:,j,m) = reshape(ar, [nchan nchan cfg.order]);

    %---compute noise covariance and rescale if necessary
    tmpnoisecov     = pe(:,nchan*cfg.order+1:nchan*(cfg.order+1));
    if dozscore,
      noisecov(:,:,j,m) = diag(datstd)*tmpnoisecov*diag(datstd);
    else
      noisecov(:,:,j,m) = tmpnoisecov;
    end
    dof(:,j)     = length(tmpdata.trial);
  end

end  
progress('close');

if dojack,
  error('this has to be implemented still');
  %---allocate memory
  jh         = complex(zeros(ntrl, ncmb,  nfoi, ntoi)+nan, zeros(ntrl, ncmb,  nfoi, ntoi)+nan);
  ja         = complex(zeros(ntrl, ncmb,  nfoi, ntoi)+nan, zeros(ntrl, ncmb,  nfoi, ntoi)+nan);
  jcrsspctrm = complex(zeros(ntrl, ncmb,  nfoi, ntoi)+nan, zeros(ntrl, ncmb,  nfoi, ntoi)+nan);
  jpowspctrm = zeros(ntrl, nchan, nfoi, ntoi)+nan;
  jnoisecov  = zeros(ntrl, ncmb, ntoi)+nan;
    
  %---loop over the tois
  progress('init', cfg.feedback, 'computing jackknife estimate of MAR-model based TFR');
  for j = 1:ntoi
    tmpcfg        = [];
    tmpcfg.toilim = [cfg.toi(j)-cfg.t_ftimwin/2 cfg.toi(j)+cfg.t_ftimwin/2-1./data.fsample];
    tmpcfg.feedback = 'no';    
    tmpdata       = redefinetrial(tmpcfg, data);
    
    %---construct data-matrix
    dat = zeros(nchan, 0);
    for k = 1:length(tmpdata.trial)
      tmpdat  = tmpdata.trial{k};
      tmptime = tmpdata.time{k};
      dat     = [dat tmpdat(chanindx, :) nan*ones(nchan, cfg.order)];
    end
  
    %---locate trial boundaries
    tmp    = isnan(dat(1,:));
    begsmp = [find(diff([1 tmp])==-1) size(dat,2)+1];
    endsmp = find(diff([tmp 1])==1);
    njtrl  = length(endsmp);
   
    %---perform jackknife
    for k = 1:njtrl
      jdat = dat;
      jdat(:,begsmp(k):endsmp(k)) = nan;
      
      jdat = jdat';
      
      [jar,jrc,jpe]     = mvar(jdat, cfg.order, cfg.mvarmethod);
      tmpnoisecov       = jpe(:,nchan*cfg.order+1:nchan*(cfg.order+1));
      if dozscore,
        jnoisecov(k,:,j)  = reshape(diag(datstd)*tmpnoisecov*diag(datstd), [ncmb 1]);
      else
        jnoisecov(k,:,j)  = reshape(tmpnoisecov, [ncmb 1]);
      end
      [tmpjh, tmpja]    = ar2h(jar, cfg.foi, data.fsample);
      jh(k,:,:,j)       = reshape(tmpjh, [ncmb nfoi]);
      ja(k,:,:,j)       = reshape(tmpja, [ncmb nfoi]);

      %---compute cross-spectra
      tmpnoisecov = reshape(jnoisecov(k,:,j), [nchan nchan]);
      tmpjcrsspctrm = complex(zeros(ncmb, nfoi), zeros(ncmb, nfoi));
      tmpjpowspctrm = complex(zeros(nchan, nfoi), zeros(nchan, nfoi));
      for m = 1:nfoi
        tmph               = reshape(jh(k,:,m,j), [nchan nchan]);
        tmpcrsspctrm       = tmph*tmpnoisecov*tmph';
        tmpjcrsspctrm(:,m) = tmpcrsspctrm(:);
        tmpjpowspctrm(:,m) = abs(diag(tmpcrsspctrm));
      end
      jcrsspctrm(k, :, :, j) = tmpjcrsspctrm;
      jpowspctrm(k, :, :, j) = tmpjpowspctrm;
    end  
  end
  progress('close');
end

%---create output-structure
mvardata          = [];
mvardata.label    = label;
mvardata.time     = cfg.toi;
mvardata.fsampleorig = data.fsample;

if dojack,
  mvardata.dimord    = 'rpt_chan_chan_lag_time';
  mvardata.coeffs    = jcoeffs;
  mvardata.noisecov  = jnoisecov;
  mvardata.method    = 'jackknife';
else
  mvardata.dimord   = 'chan_chan_lag_time';
  mvardata.coeffs   = coeffs;
  mvardata.noisecov = noisecov;
  mvardata.dof      = dof;
  mvardata.method   = 'alltrial';
end

try,
  cfg.previous = data.cfg;
end
mvardata.cfg     = cfg; 

%---SUBFUNCTION to ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
function [tap] = double_dpss(a, b, varargin);
tap = dpss(double(a), double(b), varargin{:});
