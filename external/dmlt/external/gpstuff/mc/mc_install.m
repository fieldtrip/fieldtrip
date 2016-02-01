function mc_install

    
    if ispc   % A windows version of Matlab
        mex -O -output bbmean winCSource\bbmean.c
        mex -O -output resampres winCsource\resampres.c winCsource\binsgeq.c
        mex -O -output resampsim winCsource\resampsim.c winCsource\binsgeq.c
        mex -O -output resampstr winCsource\resampstr.c winCsource\binsgeq.c
        mex -O -output resampdet winCsource\resampdet.c winCsource\binsgeq.c        
    else
        mex -O -output bbmean linuxCsource/bbmean.c
        mex -O -output resampres linuxCsource/resampres.c linuxCsource/binsgeq.c
        mex -O -output resampsim linuxCsource/resampsim.c linuxCsource/binsgeq.c
        mex -O -output resampstr linuxCsource/resampstr.c linuxCsource/binsgeq.c
        mex -O -output resampdet linuxCsource/resampdet.c linuxCsource/binsgeq.c
    end
