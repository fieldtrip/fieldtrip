function diag_install

    if ispc   % A windows version of Matlab
        mex -O -output bbprctile winCsource\bbprctile.c
    else
        mex -O -output bbprctile linuxCsource/bbprctile.c
    end