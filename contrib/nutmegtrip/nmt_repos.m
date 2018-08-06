function nmt_repos(pos,time)
% nmt_repos([x y z], t)
% nmt_repos([x y z], [t0 t1])
% nmt_repos(pos_idx, t_idx)
% nmt_repos(pos_idx, [t0_idx t1_idx])
%
% moves to selected voxel/time, then resyncs functional data

global st

if(length(pos)==1) % voxel/time indices, rather than coordinates
    pos = st.nmt.pos(pos,:);
    if(exist('time','var'))
        st.nmt.cfg.time_idx = time;
    end
else
    if(exist('time','var')) % search for index corresponding to requested time
        st.nmt.cfg.time_idx = (dsearchn(st.nmt.time',time'))';
    end
end

if(length(st.nmt.cfg.time_idx) == 1)
    st.nmt.cfg.time_idx(2) = st.nmt.cfg.time_idx(1);
end

spm_orthviews('Reposition',pos);
nmt_image('shopos');
