function [simulated] = ft_connectivitysimulation(cfg)

% FT_CONNECTIVITYSIMULATION simulates channel-level time-series data with a
% specified connectivity structure. This function returns an output data
% structure that resembles the output of FT_PREPROCESSING.
%
% Use as
%   [data] = ft_connectivitysimulation(cfg)
%
% where the configuration structure should contain:
%   cfg.method      = string, can be 'linear_mix', 'mvnrnd', 'ar', 'ar_reverse' (see below)
%   cfg.nsignal     = scalar, number of signals
%   cfg.ntrials     = scalar, number of trials
%   cfg.triallength = in seconds
%   cfg.fsample     = in Hz
%
% Depending on the specific method that is selected, the configuration
% may also contain:
%
% Method 'linear_mix' implements a linear mixing with optional time shifts
% where the number of unobserved signals can be different from the number
% of observed signals
%
% Required cfg options:
%   cfg.mix    = matrix, [nsignal x number of unobserved signals]
%                specifying the mixing from the unobserved signals to
%                the observed signals, or
%              = matrix, [nsignal x number of unobserved signals x number of
%                samples] specifying the mixing from the
%                unobserved signals to the observed signals which
%                changes as a function of time within the trial
%              = cell-arry, [1 x ntrials] with each cell a matrix as
%                specified above, when a trial-specific mixing is
%                required
%   cfg.delay  = matrix, [nsignal x number of unobserved signals]
%                specifying the time shift (in samples) between the
%                unobserved signals and the observed signals
%
% Optional cfg options:
%   cfg.bpfilter  = 'yes' (or 'no')
%   cfg.bpfreq    = [bplow bphigh] (default: [15 25])
%   cfg.demean    = 'yes' (or 'no')
%   cfg.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.absnoise  = scalar (default: 1), specifying the standard deviation of
%                   white noise superimposed on top of the simulated signals
%   cfg.randomseed = 'yes' or a number or vector with the seed value (default = 'yes')
%
% Method 'mvnrnd' implements a linear mixing with optional timeshifts in
% where the number of unobserved signals is equal to the number of observed
% signals. This method used the MATLAB function mvnrnd. The implementation
% is a bit ad-hoc and experimental, so users are discouraged to apply it.
% The time shift occurs only after the linear mixing, so the effect of the
% parameters on the simulation is not really clear. This method will be
% disabled in the future.
%
% Required cfg options:
%   cfg.covmat    = covariance matrix between the signals
%   cfg.delay     = delay vector between the signals in samples
%
% Optional cfg options:
%   cfg.bpfilter  = 'yes' (or 'no')
%   cfg.bpfreq    = [bplow bphigh] (default: [15 25])
%   cfg.demean    = 'yes' (or 'no')
%   cfg.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.absnoise  = scalar (default: 1), specifying the standard
%                   deviation of white noise superimposed on top
%                   of the simulated signals
%
% Method 'ar' implements a multivariate autoregressive model to generate
% the data.
%
% Required cfg options:
%   cfg.params   = matrix, [nsignal x nsignal x number of lags] specifying the
%                  autoregressive coefficient parameters. A non-zero
%                  element at cfg.params(i,j,k) means a
%                  directional influence from signal j onto
%                  signal i (at lag k).
%   cfg.noisecov = matrix, [nsignal x nsignal] specifying the covariance
%                  matrix of the innovation process
%
% Method 'ar_reverse' implements a multivariate autoregressive
% autoregressive model to generate the data, where the model coefficients
% are reverse-computed, based on the interaction pattern specified.
%
% Required cfg options:
%   cfg.coupling = nxn matrix, specifying coupling strength, rows causing
%                   column
%   cfg.delay    = nxn matrix, specifying the delay, in seconds, from one
%                   signal's spectral component to the other signal, rows
%                   causing column
%   cfg.ampl     = nxn matrix, specifying the amplitude 
%   cfg.bpfreq   = nxnx2 matrix, specifying the lower and upper frequencies
%                   of the bands that are transmitted, rows causing column
%
% The generated signals will have a spectrum that is 1/f + additional
% band-limited components, as specified in the cfg.
%
% See also FT_FREQSIMULATION, FT_DIPOLESIMULATION, FT_SPIKESIMULATION,
% FT_CONNECTIVITYANALYSIS

% Copyright (C) 2009-2015, Donders Institute for Brain, Cognition and Behaviour
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
ft_preamble provenance
ft_preamble randomseed
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check input configuration for the generally applicable options
cfg = ft_checkconfig(cfg, 'required', {'nsignal' 'ntrials' 'triallength' 'fsample' 'method'});
cfg = ft_checkconfig(cfg, 'rename',   {'blc', 'demean'});

% method specific defaults
switch cfg.method
  case {'ar'}
    cfg.absnoise = ft_getopt(cfg, 'absnoise', zeros(cfg.nsignal,1));
    cfg          = ft_checkconfig(cfg, 'required', {'params' 'noisecov'});
  case {'linear_mix'}
    cfg.bpfilter = ft_getopt(cfg, 'bpfilter', 'yes');
    cfg.bpfreq   = ft_getopt(cfg, 'bpfreq',   [15 25]);
    cfg.demean   = ft_getopt(cfg, 'demean',   'yes');
    cfg.absnoise = ft_getopt(cfg, 'absnoise', 1);
    cfg          = ft_checkconfig(cfg, 'required', {'mix' 'delay'});
  case {'mvnrnd'}
    cfg.bpfilter = ft_getopt(cfg, 'bpfilter', 'yes');
    cfg.bpfreq   = ft_getopt(cfg, 'bpfreq',   [15 25]);
    cfg.demean   = ft_getopt(cfg, 'demean',   'yes');
    cfg.absnoise = ft_getopt(cfg, 'absnoise', 1);
    cfg          = ft_checkconfig(cfg, 'required', {'covmat' 'delay'});
  case {'ar_reverse'}
    % reverse engineered high order ar-model
    cfg = ft_checkconfig(cfg, 'required', {'coupling' 'delay' 'ampl' 'bpfreq'});
  otherwise
end

trial = cell(1, cfg.ntrials);
time  = cell(1, cfg.ntrials);
nsmp  = round(cfg.triallength*cfg.fsample);
tim   = (0:nsmp-1)./cfg.fsample;

% create the labels
label = cell(cfg.nsignal,1);
for k = 1:cfg.nsignal
  label{k,1} = ['signal',num2str(k, '%03d')];
end

switch cfg.method
  case {'ar'}

    nlag    = size(cfg.params,3);
    nsignal = cfg.nsignal;
    params  = zeros(nlag*nsignal, nsignal);
    for k = 1:nlag
      %params(((k-1)*nsignal+1):k*nsignal,:) = cfg.params(:,:,k);
      params(((k-1)*nsignal+1):k*nsignal,:) = cfg.params(:,:,k)';
      % Use the transposition to make the implementation consistent with what
      % comes out of ft_mvaranalysis. The transposition is introduced on May
      % 13, 2011. This swaps the directional influence for existing scripts.
    end
    for k = 1:cfg.ntrials
      tmp   = zeros(nsignal, nsmp+ceil(nlag*1.05));
      noise  = mvnrnd(zeros(nsignal,1), cfg.noisecov, ceil(nsmp+nlag*1.05))';
      state0 = zeros(nsignal*nlag, 1);
      for m = 1:nlag
        indx = ((m-1)*nsignal+1):m*nsignal;
        state0(indx) = params(indx,:)'*noise(:,m);
      end
      tmp(:,1:nlag) = flip(reshape(state0, [nsignal nlag]),2);

      for m = (nlag+1):(nsmp+ceil(nlag*1.05))
        state0    = reshape(flip(tmp(:,(m-nlag):(m-1)),2), [nlag*nsignal 1]);
        tmp(:, m) = params'*state0 + noise(:,m);
      end

      trial{k} = tmp(:,(ceil(nlag*1.05)+1):end);
      if any(cfg.absnoise>0)
        trial{k} = trial{k} + diag(cfg.absnoise)*randn(size(trial{k}));
      end
      time{k}  = tim;
    end

    % create the output data
    simulated         = [];
    simulated.trial   = trial;
    simulated.time    = time;
    simulated.fsample = cfg.fsample;
    simulated.label   = label;

  case {'linear_mix'}

    fltpad = 50; %hard coded to avoid filtering artifacts
    delay  = cfg.delay;
    delay  = delay - min(delay(:)); %make explicitly >= 0
    maxdelay = max(delay(:));

    if iscell(cfg.mix)
      %each trial has different mix
      mix = cfg.mix;
    else
      %make cell-array out of mix
      tmpmix = cfg.mix;
      mix    = cell(1,cfg.ntrials);
      for tr = 1:cfg.ntrials
        mix{1,tr} = tmpmix;
      end
    end

    nmixsignal = size(mix{1}, 2); %number of "mixing signals"
    nsignal    = size(mix{1}, 1);

    if numel(size(mix{1}))==2
      %mix is static, no function of time
      for tr = 1:cfg.ntrials
        mix{tr} = mix{tr}(:,:,ones(1,nsmp+maxdelay));
      end
    elseif numel(size(mix{1}))==3 && size(mix{1},3)==nsmp
      %mix changes with time
      for tr = 1:cfg.ntrials
        mix{tr} = cat(3,mix{tr},mix{tr}(:,:,nsmp*ones(1,maxdelay)));
      end
      %FIXME think about this
      %due to the delay the mix cannot be defined instantaneously with respect to all signals
    end

    for tr = 1:cfg.ntrials
      mixsignal = randn(nmixsignal,  nsmp + 2*fltpad + maxdelay);
      mixsignal = preproc(mixsignal, label, offset2time(-fltpad, cfg.fsample, size(mixsignal,2)), cfg, fltpad, fltpad);
      tmp       = zeros(cfg.nsignal, nsmp);
      for i=1:cfg.nsignal
        for j=1:nmixsignal
          begsmp   = 1    + delay(i,j);
          endsmp   = nsmp + delay(i,j);
          tmpmix   = reshape(mix{tr}(i,j,:),[1 nsmp+maxdelay]) .* mixsignal(j,:);
          tmp(i,:) = tmp(i,:) + tmpmix(begsmp:endsmp);
        end
      end
      trial{tr} = tmp;

      % add some noise
      trial{tr} = ft_preproc_baselinecorrect(trial{tr} + cfg.absnoise*randn(size(trial{tr})));

      % define time axis for this trial
      time{tr}  = tim;
    end

  case {'mvnrnd'}
    fltpad = 100; %hard coded

    shift = max(cfg.delay(:,1)) - cfg.delay(:,1);
    for k = 1:cfg.ntrials
      % create the multivariate time series plus some padding
      tmp = mvnrnd(zeros(1,cfg.nsignal), cfg.covmat, nsmp+2*fltpad+max(shift))';

      % add the delays
      newtmp = zeros(cfg.nsignal, nsmp+2*fltpad);
      for kk = 1:cfg.nsignal
        begsmp =      + shift(kk) + 1;
        endsmp = nsmp + 2*fltpad + shift(kk);
        newtmp(kk,:) = ft_preproc_baselinecorrect(tmp(kk,begsmp:endsmp));
      end

      % apply preproc
      newtmp = preproc(newtmp, label, offset2time(-fltpad, cfg.fsample, size(newtmp,2)), cfg, fltpad, fltpad);

      trial{k} = newtmp;

      % add some noise
      trial{k} = ft_preproc_baselinecorrect(trial{k} + cfg.absnoise*randn(size(trial{k})));

      % define time axis for this trial
      time{k}  = tim;
    end
    
    % create the output data
    simulated         = [];
    simulated.trial   = trial;
    simulated.time    = time;
    simulated.fsample = cfg.fsample;
    simulated.label   = label;

  case 'ar_reverse'
    % generate a spectral transfer matrix, and a cross-spectral matrix
    % according to the specifications
    
    % predefine some variables
    fstep = 1/5;
    fs    = cfg.fsample;
    Nyq   = fs./2;
    foi   = (0:fstep:Nyq);
    omega = foi./fs;
    n     = numel(foi);
    
    % local renaming
    nsignal = cfg.nsignal;
    fband   = cfg.bpfreq;
    coupling = cfg.coupling;
    ampl     = cfg.ampl;
    delay    = cfg.delay;
    
    % create a 1/f spectrum
    slope    = 0.5;
    oneoverf = sqrt(max(omega(2)./10,omega).^-slope); % takes sqrt for amplitude
    oneoverf = oneoverf./oneoverf(1);
    %oneoverf(1) = 0;
    %z = firws_filter(5.*fs, fs, Nyq./1.01);
    %z = z(1:numel(foi));%.*exp(-1i.*pi.*foi.*rand(1)./100);
    %oneoverf = z.*oneoverf;
    
    % convert into indices
    findx = fband;
    for k = 1:numel(fband)
      if isfinite(fband(k))
        findx(k) = nearest(foi, fband(k));
      end
    end
    
    % allocate some memory
    mask = false(nsignal, nsignal, n);
    krn = zeros(size(mask));
    phi = zeros(size(krn));
    dat = zeros(size(krn));
    coupling_ampl = zeros(size(krn));
    
    for k = 1:nsignal
      for m = 1:nsignal
        if all(isfinite(squeeze(findx(k,m,:))))
          mask(k,m,findx(k,m,1):findx(k,m,2)) = true;
        end
        krn(k,m,mask(k,m,:))  = hanning(sum(mask(k,m,:)))';
        
        phi(k,m,:) = 2.*pi.*delay(k,m).*foi;
        %phi(k,m,:) = phi(k,m,:).*mask(k,m,:);
        %phi(k,m,mask(k,m,:)) = phi(k,m,mask(k,m,:))-mean(phi(k,m,mask(k,m,:)));
        if all(isfinite(squeeze(findx(k,m,:))))
          phi(k,m,1:findx(k,m,1)) = phi(k,m,findx(k,m,1));
          phi(k,m,findx(k,m,2):end) = phi(k,m,findx(k,m,2));
          phi(k,m,:) = phi(k,m,:)-mean(phi(k,m,:));
        end
        
        coupling_ampl(k,m,:) = coupling(k,m).*krn(k,m,:);
      end
    end
    
    % this matrix contains the intrinsic amplitude spectra on the diagonal
    for k = 1:nsignal
      if all(isfinite(squeeze(fband(k,k,:))))      
        z = firws_filter((1/fstep).*fs, fs, [fband(k,k,1) fband(k,k,2)]);
        z = z(1:numel(foi));%.*exp(-1i.*pi.*foi.*rand(1)./100); 
        z = z.*ampl(k,k);
        
        plateau = nearest(foi,fband(k,k,1)):nearest(foi,fband(k,k,2));
        oneoverf(plateau) = mean(abs(oneoverf(plateau)));
        dat(k,k,:) = -(abs(oneoverf)+abs(z)).*exp(1i.*(angle(z)+angle(oneoverf)));
      else
        dat(k,k,:) = oneoverf;
      end
    end
    
    % now we can create a spectral transfer matrix
    tf = zeros(nsignal,nsignal,n)+1i.*zeros(nsignal,nsignal,n);
    for k = 1:nsignal
      for m = 1:nsignal
        if k~=m && all(isfinite(squeeze(fband(k,m,:))))
          z = firws_filter((1/fstep).*fs, fs, [fband(k,m,1) fband(k,m,2)]);
          z = z(1:numel(foi));
          tf(m,k,:) = coupling(k,m).*exp(-1i.*phi(k,m,:)).*shiftdim(z,-1); % deliberate index swap!
        
        elseif k==m
          tf(k,m,:) = dat(k,m,:);
        end
      end
    end
    
    % create the cross spectral matrix
    c = zeros(size(tf));
    for k = 1:n
      c(:,:,k) = tf(:,:,k)*tf(:,:,k)'; % assume noise to be I, i.e. the tf to swallow the amplitudes
    end
    
    % scale the Nyquist and DC bins
    c(:,:,1)   = real(c(:,:,1)./2);
    c(:,:,end) = real(c(:,:,end)./2);
    
    % create a freq-structure
    freq           = [];
    freq.crsspctrm = c;
    freq.label     = label;
    freq.freq      = foi;
    freq.dimord    = 'chan_chan_freq';
   
    % estimate the transfer-matrix non-parametrically
    tmpcfg        = [];
    tmpcfg.method = 'transfer';
    tmpcfg.granger.stabilityfix = true;
    t             = ft_connectivityanalysis(tmpcfg, freq);
    
    % estimate the ar-model coefficients
     a = transfer2coeffs(t.transfer,t.freq);
    
    % recursively call this function to generate the data, this is
    % somewhate tricky with respect to keeping the provenance info. Here,
    % it is solved by removing from the cfg the original user-specified
    % fields
    cfgorig      = cfg;
    cfg          = removefields(cfgorig, {'coupling' 'ampl' 'delay' 'bpfreq'});
    cfg.method   = 'ar';
    cfg.params   = a;
    cfg.noisecov = diag(diag(t.noisecov.*cfg.fsample./2));
    simulated    = ft_connectivitysimulation(cfg);
    cfg.previous = keepfields(cfgorig, {'coupling' 'ampl' 'delay' 'bpfreq'});
    
  otherwise
    ft_error('unknown method');
end


% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble randomseed
ft_postamble provenance
ft_postamble history simulated
ft_postamble savevar simulated


%%%%%%
% helper function
function A = transfer2coeffs(H, freq, labelcmb, maxlag)

% TRANSFER2COEFFS converts a spectral transfer matrix into the time domain
% equivalent multivariate autoregressive coefficients up to a specified
% lag, starting from lag 1.

if nargin<3
  labelcmb = [];
end
if nargin<4
  maxlag = [];
end

% do a check on the input data
siz = size(H);
if numel(siz)==3 && siz(1)==siz(2)
  % assume chan_chan_freq
  isfull = true;
elseif numel(siz)==2
  % assume chancmb_freq
  isfull = false;
  %assert(~isempty(labelcmb), 'input data appears to be chancmb_freq, but labelcmb is missing');
else
  ft_error('dimensionality of input data is not supported');
end

dfreq = round(diff(freq)*1e5)./1e5; % allow for some numeric issues
if ~all(dfreq==dfreq(1))
  ft_error('the frequency axis is not evenly spaced');
end

if freq(1)~=0
  ft_warning('when converting the transfer function to coefficients, the frequency axis should ideally start at 0, zero padding the spectral density'); 
  dfreq = mean(dfreq);
  npad  = freq(1)./dfreq;
  
  % update the freq axis and keep track of the frequency bins that are
  % expected in the output
  selfreq  = (1:numel(freq)) + npad;
  freq     = [(0:(npad-1))./dfreq freq];
  if isfull
    H = cat(3, zeros(siz(1),siz(2),npad), H);
  else
    H = cat(2, zeros(siz(1),npad), H);
  end
else
  selfreq  = 1:numel(freq);
end

% ensure H to be double precision
H = double(H);

% deal with the two different types of input
if isfull
  % check whether the last frequency bin is strictly real-valued.
  % if that's the case, then it is assumed to be the Nyquist frequency
  % and the two-sided spectral density will have an even number of
  % frequency bins. if not, in order to preserve hermitian symmetry,
  % the number of frequency bins needs to be odd.
  Hend = H(:,:,end);
  N    = numel(freq);
  m    = size(H,1);
  if all(imag(Hend(:))<abs(trace(Hend)./size(Hend,1)*1e-9))
    N2 = 2*(N-1);
  else
    N2 = 2*(N-1)+1;
  end

  % preallocate memory for efficiency
  Harr   = zeros(m,m,N2) + 1i.*zeros(m,m,N2);

  % the input cross-spectral density is assumed to be weighted with a
  % factor of 2 in all non-DC and Nyquist bins, therefore weight the
  % DC-bin with a factor of sqrt(2) to get a correct two-sided representation
  Harr(:,:,1) = H(:,:,1).*2;
  for k = 2:N
    Harr(:,:,       k) = H(:,:,k);
    Harr(:,:,(N2+2)-k) = conj(H(:,:,k));
  end
  
  % the input cross-spectral density is assumed to be weighted with a
  % factor of 2 in all non-DC and Nyquist bins, therefore weight the
  % Nyquist bin with a factor of sqrt(2) to get a correct two-sided representation
  if mod(size(Harr,3),2)==0
    Harr(:,:,N) = Harr(:,:,N).*sqrt(2);
  end
  
  % invert the transfer matrix to get the fourier representation of the
  % coefficients, and add an identity matrix 
  I = eye(siz(1));
  for k = 1:size(Harr,3)
    Harr(:,:,k) = I-inv(Harr(:,:,k));
  end
  
  % take the inverse fft to get the coefficients
  A = ifft(reshape(permute(Harr, [3 1 2]), N2, []), 'symmetric');
  A = A(2:end,:);
  A = ipermute(reshape(A, [N2-1 siz(1) siz(1)]), [3 1 2]);
  
  if ~isempty(maxlag)
    A = A(:,:,1:maxlag);
  end
else
  % check whether the last frequency bin is strictly real-valued.
  % if that's the case, then it is assumed to be the Nyquist frequency
  % and the two-sided spectral density will have an even number of
  % frequency bins. if not, in order to preserve hermitian symmetry,
  % the number of frequency bins needs to be odd.
  Hend = H(:,end);
  N    = numel(freq);
  m    = size(H,1);
  if all(imag(Hend(:))<max(abs(Hend))*1e-9)
    % the above heuristic may be a bit silly, FIXME
    N2 = 2*(N-1);
  else
    N2 = 2*(N-1)+1;
  end

  % preallocate memory for efficiency
  Harr   = zeros(m,N2) + 1i.*zeros(m,N2);

  % the input cross-spectral density is assumed to be weighted with a
  % factor of 2 in all non-DC and Nyquist bins, therefore weight the
  % DC-bin with a factor of sqrt(2) to get a correct two-sided representation
  Harr(:,1) = H(:,1).*sqrt(2);
  for k = 2:N
    Harr(:,       k) = H(:,k);
    Harr(:,(N2+2)-k) = conj(H(:,k));
  end
  
  % the input cross-spectral density is assumed to be weighted with a
  % factor of 2 in all non-DC and Nyquist bins, therefore weight the
  % Nyquist bin with a factor of sqrt(2) to get a correct two-sided representation
  if mod(size(Harr,3),2)==0
    Harr(:,N) = Harr(:,N).*sqrt(2);
  end
  
  % invert the transfer matrix to get the fourier representation of the
  % coefficients, and add an identity matrix 
  %
  % this assumes Harr to be in the rows quadruplets of pairwise
  % decompositions, i.e. reshapable, without checking the labelcmb
  ncmb = size(Harr,1)./4;
  I = eye(2);
  for k = 1:N2
    Htmp = reshape(Harr(:,k), [2 2 ncmb]);
    Htmp = repmat(I, [1 1 ncmb]) - inv2x2(Htmp);
    Harr(:,k) = Htmp(:);
  end
    
  % take the inverse fft to get the coefficients
  A = ifft(permute(Harr, [2 1]), 'symmetric');
  A = A(2:end,:);
  A = ipermute(A, [2 1]);
  
  if ~isempty(maxlag)
    A = A(:,1:maxlag);
  end

end

function z = firws_filter(N, Fs, Fbp)

switch numel(Fbp)
  case 1
    [~, B, ~] = ft_preproc_lowpassfilter(randn(1,N), Fs, Fbp, [], 'firws', 'onepass-minphase');
    z  = fft(B, N);

  case 2
    [~, B, ~] = ft_preproc_bandpassfilter(randn(1,N), Fs, Fbp, [], 'firws', 'onepass-minphase');
    z  = fft(B, N);
  
end
