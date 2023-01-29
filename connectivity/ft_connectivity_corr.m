function [c, v, outcnt] = ft_connectivity_corr(inputdata, varargin)

% FT_CONNECTIVITY_CORR computes correlation, coherence or a related quantity from a
% data-matrix containing a covariance or cross-spectral density. This implements the
% methods as described in the following papers:
%
% Coherence: Rosenberg et al, The Fourier approach to the identification of
% functional coupling between neuronal spike trains. Prog Biophys Molec
% Biol 1989; 53; 1-31
%
% Partial coherence: Rosenberg et al, Identification of patterns of
% neuronal connectivity - partial spectra, partial coherence, and neuronal
% interactions. J. Neurosci. Methods, 1998; 83; 57-72
%
% Phase locking value: Lachaux et al, Measuring phase sychrony in brain
% signals. Human Brain Mapping, 1999; 8; 194-208.
%
% Imaginary part of coherency: Nolte et al, Identifying true brain
% interaction from EEG data using the imaginary part of coherence. Clinical
% Neurophysiology, 2004; 115; 2292-2307
%
% Use as
%   [c, v, n] = ft_connectivity_corr(inputdata, ...)
%
% The input data should be a covariance or cross-spectral density array organized as
%   Repetitions x Channel x Channel (x Frequency) (x Time)
% or
%   Repetitions x Channelcombination (x Frequency) (x Time)
%
% If the input already contains an average, the first dimension must be singleton.
% Furthermore, the input data can be complex-valued cross spectral densities, or
% real-valued covariance estimates. If the former is the case, the output will be
% coherence (or a derived metric), if the latter is the case, the output will be the
% correlation coefficient.
%
% The output represents
%   c = the correlation/coherence
%   v = variance estimate, this can only be computed if the data contains leave-one-out samples
%   n = the number of repetitions in the input data
%
% Additional optional input arguments come as key-value pairs:
%   'dimord'    = string, specifying how the input matrix should be interpreted
%   'hasjack'   = boolean flag that specifies whether the repetitions represent leave-one-out samples
%   'complex'   = 'abs', 'angle', 'real', 'imag', 'complex', 'logabs' for post-processing of coherency
%   'powindx'   = required if the input data contain linearly indexed channel pairs. This
%                 should be an Nx2 matrix indexing on each row for the respective channel
%                 pair the indices of the corresponding auto-spectra.
%   'pownorm'   = boolean flag that specifies whether normalisation with the product
%                 of the power should be performed (thus should be true when
%                 correlation/coherence is requested, and false when covariance
%                 or cross-spectral density is requested).
%   'feedback'  = 'none', 'text', 'textbar', 'dial', 'etf', 'gui' type of feedback showing progress of computation, see FT_PROGRESS
%
% Partialisation can be performed when the input data is (chan x chan). The following
% option needs to be specified:
%   'pchanindx' = index-vector to the channels that need to be partialised
%
% See also CONNECTIVITY, FT_CONNECTIVITYANALYSIS

% Copyright (C) 2009-2010 Donders Institute, Jan-Mathijs Schoffelen
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

% FiXME: If output is angle, then jack-knifing should be done
% differently since it is a circular variable

hasjack     = ft_getopt(varargin, 'hasjack',  false);
cmplx       = ft_getopt(varargin, 'complex',  'abs');
feedback    = ft_getopt(varargin, 'feedback', 'none');
dimord      = ft_getopt(varargin, 'dimord');
powindx     = ft_getopt(varargin, 'powindx');
pownorm     = ft_getopt(varargin, 'pownorm', 0);
pchanindx   = ft_getopt(varargin, 'pchanindx');

if isempty(dimord)
  ft_error('input parameters should contain a dimord');
end

siz = [size(inputdata) 1];

% do partialisation if necessary
if ~isempty(pchanindx) && isempty(powindx)
  % partial spectra are computed as in Rosenberg JR et al (1998) J.Neuroscience Methods, equation 38
  
  npchanindx = numel(pchanindx);
  chan   = setdiff(1:size(inputdata,2), pchanindx);
  nchan  = numel(chan);
  newsiz = siz;
  newsiz(2:3) = numel(chan); % size of partialised csd
  
  A = zeros(newsiz);
  
  for j = 1:siz(1) % loop over rpt
    AA = reshape(inputdata(j, chan,  chan, : ),         [nchan  nchan      prod(siz(4:end))]); % fold freq_time into one dimension
    AB = reshape(inputdata(j, chan,  pchanindx,: ),     [nchan  npchanindx prod(siz(4:end))]);
    BA = reshape(inputdata(j, pchanindx, chan, : ),     [npchanindx nchan  prod(siz(4:end))]);
    BB = reshape(inputdata(j, pchanindx, pchanindx, :), [npchanindx npchanindx prod(siz(4:end))]);
    for k = 1:prod(siz(4:end)) % loop over freq or freq_time
      A(j,:,:,k) = AA(:,:,k) - AB(:,:,k)/(BB(:,:,k))*BA(:,:,k);
    end
  end
  inputdata = A;
  siz = [size(inputdata) 1];
  
elseif ~isempty(pchanindx)
  % linearly indexed crossspectra require some more complicated handling
  if numel(pchanindx)>1
    ft_error('more than one channel for partialisation with linearly indexed crossspectra is currently not implemented');
  end
  
  p_input = inputdata;
  for k = 1:size(powindx,1)
    % we need to look for the combi's (and take conjugates if needed), to
    % achieve F_ab\p = F_ab - F_ap*inv(F_p)*F_pb;
    
    this = powindx(k,:);
    sela = powindx(:,1)==this(1)&powindx(:,2)==pchanindx;
    if any(sela)
      F_ap = inputdata(:,sela,:,:);
    else
      sela = powindx(:,2)==this(1)&powindx(:,1)==pchanindx;
      F_ap = conj(inputdata(:,sela,:,:));
    end
    
    selb = powindx(:,2)==this(2)&powindx(:,1)==pchanindx;
    if any(selb)
      F_pb = inputdata(:,selb,:,:);
    else
      selb = powindx(:,1)==this(2)&powindx(:,2)==pchanindx;
      F_pb = conj(inputdata(:,selb,:,:));
    end
    selp = powindx(:,1)==pchanindx&powindx(:,2)==pchanindx;
    F_pp = inputdata(:,selp,:,:);
    
    p_input(:,k,:,:) = inputdata(:,k,:,:) - F_ap.*(1./F_pp).*F_pb;
    
  end
  inputdata = p_input; clear p_input;
else
  % do nothing
end

% compute the metric
if (length(strfind(dimord, 'chan'))~=2 || contains(dimord, 'pos')) && ~isempty(powindx)
  % crossterms are not described with chan_chan_therest, but are linearly indexed
  
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  outcnt = zeros(siz(2:end));
  ft_progress('init', feedback, 'computing metric...');
  for j = 1:siz(1)
    ft_progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    if pownorm
      p1    = reshape(inputdata(j,powindx(:,1),:,:,:), siz(2:end));
      p2    = reshape(inputdata(j,powindx(:,2),:,:,:), siz(2:end));
      denom = sqrt(p1.*p2); clear p1 p2
    else
      denom = 1;
    end
    tmp    = complexeval(reshape(inputdata(j,:,:,:,:), siz(2:end))./denom, cmplx);
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
    outcnt = outcnt + double(~isnan(tmp));
  end
  ft_progress('close');
  
elseif length(strfind(dimord, 'chan'))==2 || length(strfind(dimord, 'pos'))==2
  % crossterms are described by chan_chan_therest
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  outcnt = zeros(siz(2:end));
  ft_progress('init', feedback, 'computing metric...');
  for j = 1:siz(1)
    ft_progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    if pownorm
      p1  = zeros([siz(2) 1 siz(4:end)]);
      p2  = zeros([1 siz(3) siz(4:end)]);
      for k = 1:siz(2)
        p1(k,1,:,:,:,:) = inputdata(j,k,k,:,:,:,:);
        p2(1,k,:,:,:,:) = inputdata(j,k,k,:,:,:,:);
      end
      p1    = p1(:,ones(1,siz(3)),:,:,:,:);
      p2    = p2(ones(1,siz(2)),:,:,:,:,:);
      denom = sqrt(p1.*p2); clear p1 p2;
    else
      denom = 1;
    end
    tmp    = complexeval(reshape(inputdata(j,:,:,:,:,:,:), siz(2:end))./denom, cmplx); % added this for nan support marvin
    %tmp(isnan(tmp)) = 0; % added for nan support
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
    outcnt = outcnt + double(~isnan(tmp));
  end
  ft_progress('close');
  
else
  ft_error('unsupported dimord "%s"', dimord);
end

n  = siz(1);
if all(outcnt(:)==n)
  outcnt = n;
end

%n1 = shiftdim(sum(~isnan(input),1),1);
%c  = outsum./n1; % added this for nan support marvin
c = outsum./outcnt;

% correct the variance estimate for the under-estimation introduced by the jackknifing
if n>1
  if hasjack
    %bias = (n1-1).^2; % added this for nan support marvin
    bias = (outcnt-1).^2;
  else
    bias = 1;
  end
  %v = bias.*(outssq - (outsum.^2)./n1)./(n1 - 1); % added this for nan support marvin
  v = bias.*(outssq - (outsum.^2)./outcnt)./(outcnt-1);
else
  v = [];
end

function [c] = complexeval(c, str)

switch str
  case 'complex'
    % do nothing
  case 'abs'
    c = abs(c);
  case 'angle'
    c = angle(c); % negative angle means first row leads second row
  case 'imag'
    c = imag(c);
  case 'absimag'
    c = abs(imag(c));
  case 'real'
    c = real(c);
  case '-logabs'
    c = -log(1 - abs(c).^2);
  otherwise
    ft_error('complex = ''%s'' not supported', str);
end
