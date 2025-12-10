function test_issue2085

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_resampledata ft_specest_irasa
% DATA no

% test function to evaluate the reported scale difference when resampling
% with the drop in compat replacement of the resample function (and lower
% level dependencies).

global ft_default
[ftver, ftpath] = ft_version;


% simulate data
t = (1:1000)/1000;
fn = cumsum(randn(1,length(t))); 
in_dat = fn + cos(2*pi*10*t) + cos(2*pi*60*t);

in_dat = in_dat - mean(in_dat);

% resample, using a variety of p/q ratios
p = [1:50 1:50]; q = [10*ones(1,50) 11*ones(1,50)];

restoredefaultpath
addpath(ftpath);
ft_default.toolbox.signal = 'compat';
ft_defaults;
fprintf(sprintf('resampling with: %s\n', which('resample')));
for k = 1:numel(p)
  out_ft_rspu{k,1} = resample(in_dat, p(k), q(k));
  out_ft_rspd{k,1} = resample(in_dat, q(k), p(k));
end

restoredefaultpath
addpath(ftpath);
ft_default.toolbox.signal = 'matlab';
ft_defaults;
fprintf(sprintf('resampling with: %s\n', which('resample')));
for k = 1:numel(p)
  out_mat_rspu{k,1} = resample(in_dat, p(k), q(k));
  out_mat_rspd{k,1} = resample(in_dat, q(k), p(k));
end
  
for k = 1:numel(p)
  bd(k,1) = out_ft_rspd{k}/out_mat_rspd{k};
  bu(k,1) = out_ft_rspu{k}/out_mat_rspu{k};
end
figure;plot(bu);hold on;plot(bd);

assert(all(bu>0.999) && all(bd>0.999));
% there's a tiny scale difference between the legacy matlab and the compat
% version, which could be titrated to be closer to 1, but that would not be
% based on a mathematically informed choice, so I think we may accept this <0.1% amplitude diff.


% figure(); hold on;
% plot(in_dat, '-b');
% plot(out_ft_rspu, '-r');
% plot(out_mat_rspu, '-g');
% legend({'original data', 'ft resampled', 'mat resampled'});
% title('upsample')
% hold off;
% 
% figure(); hold on;
% plot(in_dat, '-b');
% plot(out_ft_rspd, '-r');
% plot(out_mat_rspd, '-g');
% legend({'original data', 'ft resampled', 'mat resampled'});
% title('downsample');
% hold off;