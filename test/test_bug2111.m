function test_bug2111

% TEST test_bug2111 
% TEST ft_prepare_layout

% ORIGINAL BUGREPORT:
% Hi,
% 
% I'm getting a weird results projecting a CTF grad to 2D. It looks like there is
% a mismatch between locations and labels (as the reference sensors are also
% there and they shouldn't be) and also some kind of misalignment. I attach the
% cfg.
% 
% Thanks,
% 
% Vladimir

load test_bug2111.mat

cfg.feedback = 'yes';
cfg.channel = 'MEG';
% this verifies the bug visually
ft_prepare_layout(cfg);

% conclusion: all fine, fiducials have probably been swapped



