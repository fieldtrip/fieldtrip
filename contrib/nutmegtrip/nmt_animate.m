function nmt_animate(cfg)
% designed to animate activations over time in nutmegtrip viewer

global st

if(~exist('cfg','var'))
    cfg=[];
end

switch(get(st.nmt.gui.animate,'String'))
    case 'Stop'
        set(st.nmt.gui.animate,'String','Animate');
        return
    case 'Animate'
        set(st.nmt.gui.animate,'String','Stop');
        
        
        Fs = 1/(st.nmt.time(2)-st.nmt.time(1));
        
        cfg.stepsize = ft_getopt(cfg, 'stepsize', round(1e-3 * Fs)); % default is 5ms step
        cfg.winsize = ft_getopt(cfg,'winsize',1);
        
        cfg.winsize = cfg.winsize - 1; % simplifies index tracking
        
        for ii = st.nmt.cfg.time_idx:cfg.stepsize:(length(st.nmt.time) - cfg.winsize)
            if(strcmp(get(st.nmt.gui.animate,'String'),'Animate'))
                break
            else
                st.nmt.cfg.time_idx = [ii ii+cfg.winsize];
                nmt_timeselect;
            end
        end
        
        set(st.nmt.gui.animate,'String','Animate');
end
