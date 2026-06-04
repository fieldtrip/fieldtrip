function test_pull2593

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_freqstatistics ft_statistics_montecarlo tfcestat findcluster
% DATA no

% This test covers the exact TFCE (eTFCE) implementation in tfcestat, added in
% pull request 2593, and its 'discrete' (original) counterpart. It runs the
% high-level ft_freqstatistics with cfg.correctm='tfce' on simulated data with
% a fixed random seed, for cfg.tfce_method = 'exact' and 'discrete', and checks:
%
%   1) the exact method equals the discrete method in the limit of many height
%      steps (exact == nsteps->Inf), and is closer to it than the default 100
%      steps -> the exact method removes the discretisation bias;
%   2) the observed TFCE statistic reproduces hard reference values (these are
%      deterministic, independent of the randomization);
%   3) the test detects the simulated effect (significant cluster);
%   4) an unsupported cfg.tfce_method errors.
%
% The reference values were determined by running this configuration once.

% ----------------------------------------------------------------------------
% simulated data (fixed seed so the outcome is reproducible)
% ----------------------------------------------------------------------------
rng(42,'twister');
nchan = 12; nfreq = 1; ntime = 15; nsubj = 10;

% 3 x 4 channel grid with 4-connected (Manhattan) neighbours
[gx,gy] = ndgrid(1:3,1:4);
pos     = [gx(:) gy(:)];
label   = arrayfun(@(k)sprintf('chan%02d',k),(1:nchan)','UniformOutput',false);
neighbours = struct('label',{},'neighblabel',{});
for i = 1:nchan
  d  = abs(pos(:,1)-pos(i,1)) + abs(pos(:,2)-pos(i,2));
  nb = find(d==1);
  neighbours(i).label       = label{i};
  neighbours(i).neighblabel = label(nb);
end

tmpl.label  = label;
tmpl.freq   = 10;
tmpl.time   = linspace(0,0.3,ntime);
tmpl.dimord = 'chan_freq_time';

sigchan = [1 2 4 5];   % a 2x2 block in one corner
sigtime = 6:10;
A = cell(1,nsubj); B = cell(1,nsubj);
for s = 1:nsubj
  a = randn(nchan,nfreq,ntime);
  b = randn(nchan,nfreq,ntime);
  b(sigchan,1,sigtime) = b(sigchan,1,sigtime) + 1.2;   % a real effect in condition B
  A{s} = tmpl; A{s}.powspctrm = a;
  B{s} = tmpl; B{s}.powspctrm = b;
end

% ----------------------------------------------------------------------------
% common statistics configuration
% ----------------------------------------------------------------------------
design = zeros(2,2*nsubj);
design(1,:) = [1:nsubj 1:nsubj];
design(2,:) = [ones(1,nsubj) 2*ones(1,nsubj)];

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'tfce';
cfg.numrandomization = 100;
cfg.randomseed       = 666;   % fixed -> reproducible randomization distribution
cfg.tail             = 0;
cfg.alpha            = 0.05;
cfg.design           = design;
cfg.uvar             = 1;
cfg.ivar             = 2;
cfg.neighbours       = neighbours;
cfg.tfce_H           = 2;
cfg.tfce_E           = 0.5;

% ----------------------------------------------------------------------------
% run the exact method and the discrete method with 100 and 2000 height steps
% ----------------------------------------------------------------------------
cfg.tfce_method = 'exact';                            stat_exact = ft_freqstatistics(cfg, A{:}, B{:});
cfg.tfce_method = 'discrete'; cfg.tfce_nsteps = 100;  stat_d100  = ft_freqstatistics(cfg, A{:}, B{:});
cfg.tfce_method = 'discrete'; cfg.tfce_nsteps = 2000; stat_d2000 = ft_freqstatistics(cfg, A{:}, B{:});

e = stat_exact.stattfce(:);
f = stat_d2000.stattfce(:);
d = stat_d100.stattfce(:);
scale = max(abs(e));

rel_d2000 = max(abs(e-f))/scale;   % distance exact <-> 2000-step discrete
rel_d100  = max(abs(e-d))/scale;   % distance exact <-> 100-step discrete

% ----------------------------------------------------------------------------
% 1) exact == discrete in the many-steps limit, and beats the 100-step default
% ----------------------------------------------------------------------------
assert(rel_d2000 < 2e-3, ...
  'exact TFCE does not match the 2000-step discrete approximation (relmax=%.3e)', rel_d2000);
assert(rel_d2000 < rel_d100, ...
  'increasing the number of discrete steps did not move the result towards the exact one (%.3e vs %.3e)', rel_d2000, rel_d100);

% ----------------------------------------------------------------------------
% 2) hard reference values for the observed exact TFCE (deterministic)
% ----------------------------------------------------------------------------
REF_SUM   = -1000.095014;
REF_MAX   =  15.92646833;
REF_MIN   = -209.8410122;
REF_SUMSQ =  106587.4964;
reltol    = 1e-5;

assert(abs(sum(e)   - REF_SUM)   <= reltol*abs(REF_SUM),   'sum of exact TFCE changed: %.10g (expected %.10g)', sum(e), REF_SUM);
assert(abs(max(e)   - REF_MAX)   <= reltol*abs(REF_MAX),   'max of exact TFCE changed: %.10g (expected %.10g)', max(e), REF_MAX);
assert(abs(min(e)   - REF_MIN)   <= reltol*abs(REF_MIN),   'min of exact TFCE changed: %.10g (expected %.10g)', min(e), REF_MIN);
assert(abs(sum(e.^2)- REF_SUMSQ) <= reltol*abs(REF_SUMSQ), 'sum-of-squares of exact TFCE changed: %.10g (expected %.10g)', sum(e.^2), REF_SUMSQ);

% ----------------------------------------------------------------------------
% 3) the simulated effect is detected (randomization-dependent, kept robust)
% ----------------------------------------------------------------------------
assert(any(stat_exact.mask(:)),        'exact TFCE did not find any significant sample for a clear simulated effect');
assert(min(stat_exact.prob(:)) <= 0.05,'exact TFCE smallest p-value is not significant for a clear simulated effect');

% ----------------------------------------------------------------------------
% 4) an unsupported method must error
% ----------------------------------------------------------------------------
cfg.tfce_method = 'nonexistent';
ok = false;
try
  ft_freqstatistics(cfg, A{:}, B{:});
catch
  ok = true;
end
assert(ok, 'an unsupported cfg.tfce_method should raise an error but did not');

fprintf('test_pull2593 passed: exact-vs-2000step relmax=%.3e, 100step relmax=%.3e\n', rel_d2000, rel_d100);
