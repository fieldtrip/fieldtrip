function [c, v, outcnt] = ft_connectivity_corr(input, varargin)

% FT_CONNECTIVITY_CORR computes correlation, coherence or a related
% quantity from a data-matrix containing a covariance or cross-spectral
% density.
%
% Use as
%   [c, v, n] = ft_connectivity_corr(input, varargin)
%
% The input data input should be an array organized as:
%   Repetitions x Channel x Channel (x Frequency) (x Time)
% or
%   Repetitions x Channelcombination (x Frequency) (x Time)
%
% If the input already contains an average, the first dimension should be
% singleton. Furthermore, the input data can be complex-valued cross
% spectral densities, or real-valued covariance estimates. If the former is
% the case, the output will be coherence (or a derived metric), if the
% latter is the case, the output will be the correlation coefficient.
%
% Additional input arguments come as key-value pairs:
%
%   hasjack   = 0 or 1 specifying whether the Repetitions represent
%               leave-one-out samples
%   complex   = 'abs', 'angle', 'real', 'imag', 'complex', 'logabs' for
%               post-processing of coherency
%   feedback  = 'none', 'text', 'textbar' type of feedback showing progress of
%               computation
%   dimord    = specifying how the input matrix should be interpreted
%   powindx   = required if the input data contain linearly indexed
%               channel pairs. should be an Nx2 matrix indexing on each
%               row for the respective channel pair the indices of the
%               corresponding auto-spectra
%   pownorm   = flag that specifies whether normalisation with the
%               product of the power should be performed (thus should
%               be true when correlation/coherence is requested, and
%               false when covariance or cross-spectral density is
%               requested).
%
% Partialisation can be performed when the input data is (chan x chan). The
% following options need to be specified:
%
%   pchanindx   = index-vector to the channels that need to be
%                 partialised
%   allchanindx = index-vector to all channels that are used
%                 (including the "to-be-partialised" ones).
%
% The output c contains the correlation/coherence, v is a variance estimate
% which only can be computed if the data contains leave-one-out samples,
% and n is the number of repetitions in the input data.
%
% It implements the methods as described in the following papers:
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
% See also FT_CONNECTIVITYANALYSIS

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

hasjack     = ft_getopt(varargin, 'hasjack', 0);
cmplx       = ft_getopt(varargin, 'complex', 'abs');
feedback    = ft_getopt(varargin, 'feedback', 'none');
dimord      = ft_getopt(varargin, 'dimord');
powindx     = ft_getopt(varargin, 'powindx');
pownorm     = ft_getopt(varargin, 'pownorm', 0);
pchanindx   = ft_getopt(varargin, 'pchanindx');
allchanindx = ft_getopt(varargin, 'allchanindx');

if isempty(dimord)
  error('input parameters should contain a dimord');
end

siz = [size(input) 1];

% do partialisation if necessary
if ~isempty(pchanindx),
  % partial spectra are computed as in Rosenberg JR et al (1998) J.Neuroscience Methods, equation 38
  
  chan   = allchanindx;
  nchan  = numel(chan);
  pchan  = pchanindx;
  npchan = numel(pchan);
  newsiz = siz;
  newsiz(2:3) = numel(chan); % size of partialised csd
  
  A  = zeros(newsiz);
  
  % FIXME this only works for data without time dimension
  if numel(siz)==5 && siz(5)>1, error('this only works for data without time'); end
  for j = 1:siz(1) %rpt loop
    AA = reshape(input(j, chan,  chan, : ), [nchan  nchan  siz(4:end)]);
    AB = reshape(input(j, chan,  pchan,: ), [nchan  npchan siz(4:end)]);
    BA = reshape(input(j, pchan, chan, : ), [npchan nchan  siz(4:end)]);
    BB = reshape(input(j, pchan, pchan, :), [npchan npchan siz(4:end)]);
    for k = 1:siz(4) %freq loop
      %A(j,:,:,k) = AA(:,:,k) - AB(:,:,k)*pinv(BB(:,:,k))*BA(:,:,k);
      A(j,:,:,k) = AA(:,:,k) - AB(:,:,k)/(BB(:,:,k))*BA(:,:,k);
    end
  end
  input = A;
  siz = size(input);
else
  % do nothing
end

% compute the metric
if (length(strfind(dimord, 'chan'))~=2 || ~isempty(strfind(dimord, 'pos'))) && ~isempty(powindx),
  % crossterms are not described with chan_chan_therest, but are linearly indexed
  
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  outcnt = zeros(siz(2:end));
  ft_progress('init', feedback, 'computing metric...');
  for j = 1:siz(1)
    ft_progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    if pownorm
      p1    = reshape(input(j,powindx(:,1),:,:,:), siz(2:end));
      p2    = reshape(input(j,powindx(:,2),:,:,:), siz(2:end));
      denom = sqrt(p1.*p2); clear p1 p2
    else
      denom = 1;
    end
    tmp    = complexeval(reshape(input(j,:,:,:,:), siz(2:end))./denom, cmplx);
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
    outcnt = outcnt + double(~isnan(tmp));
  end
  ft_progress('close');
  
elseif length(strfind(dimord, 'chan'))==2 || length(strfind(dimord, 'pos'))==2,
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
        p1(k,1,:,:,:,:) = input(j,k,k,:,:,:,:);
        p2(1,k,:,:,:,:) = input(j,k,k,:,:,:,:);
      end
      p1    = p1(:,ones(1,siz(3)),:,:,:,:);
      p2    = p2(ones(1,siz(2)),:,:,:,:,:);
      denom = sqrt(p1.*p2); clear p1 p2;
    else
      denom = 1;
    end
    tmp    = complexeval(reshape(input(j,:,:,:,:,:,:), siz(2:end))./denom, cmplx); % added this for nan support marvin
    %tmp(isnan(tmp)) = 0; % added for nan support
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
    outcnt = outcnt + double(~isnan(tmp));
  end
  ft_progress('close');
  
end

n  = siz(1);
if all(outcnt(:)==n)
  outcnt = n;
end

%n1 = shiftdim(sum(~isnan(input),1),1);
%c  = outsum./n1; % added this for nan support marvin
c = outsum./outcnt;

% correct the variance estimate for the under-estimation introduced by the jackknifing
if n>1,
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
    %do nothing
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
    error('complex = ''%s'' not supported', str);
end
