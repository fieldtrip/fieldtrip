function nmt_animate(cfg)
% designed to animate activations over time in nutmegtrip viewer

global st

cfg=[];

Fs = 1/(st.nmt.time(2)-st.nmt.time(1));

cfg.stepsize = ft_getopt(cfg, 'stepsize', round(10e-3 * Fs)); % default is 10ms step
cfg.winsize = ft_getopt(cfg,'winsize',1);

cfg.winsize = cfg.winsize - 1; % simplifies index tracking

for ii = 1:cfg.stepsize:(length(st.nmt.time) - cfg.winsize)
    st.nmt.cfg.time_idx = [ii ii+cfg.winsize];
    nmt_timeselect;
end