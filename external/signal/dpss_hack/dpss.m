function [dps_seq, lambda] = dpss(seq_length,time_halfbandwidth)

% DPSS approximate the discrete prolate spheroidal (DPSS), or Slepian sequences
% using a precomputed version on disk

persistent w
if isempty(w)
  % load the precomputed tapers only once
  load precompute_dpss
end

n = 1000;
s = 1:0.5:50;

% find the nearest match
i = nearest(s, time_halfbandwidth, true, true);

% interpolate onto the desired number of samples
dps_seq = interp1(linspace(0,1,n), w{i}, linspace(0,1,seq_length));
lambda  = ones(1, size(dps_seq, 2));

ft_warning('concentration eigenvalues cannot be computed, and are set to 1 for all tapers, this affects the correctness of the results when using adaptively weighted multitapering');
