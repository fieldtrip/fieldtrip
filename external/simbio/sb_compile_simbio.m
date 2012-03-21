function sb_compile_simbio(force)

% FT_COMPILE_VISTA is used for compiling most of the Vista MEX files that are used in FieldTrip
%
% Please note that this script does NOT set up your MEX environment for you, so in case
% you haven't installed a fortran compiler yet, see the list at
% http://www.mathworks.de/support/compilers/R2011b/win32.html
% or similar for your matlab version and operating system for supported
% compilers.

% The logic in this script is to first build a list of files that actually need compilation for the
% particular platform that Matlab is running on, and then to go through that list.
% Functions are added to the list by giving their destination directory and (relative to that) the 
% name of the source file (without the .c). Optionally, you can specify a list of platform this
% file needs to be compiled on only, and a list of platforms where you don't compile it on.
% Finally, you can give extra arguments to the MEX command, e.g., for including other c-sources or
% giving compiler flags.

% Copyright (C) 2010, Stefan Klanke, 2011 Cristiano Micheli, Johannes
% Vorwerk
%
% $Log$

% ft_compile_vista compiles the files belonging to the fileio Vista format library
% 
% vista/libvista.a
% vista/read_vista_mesh.cpp
% vista/write_vista_mesh.cpp
% vista/write_vista_vol.cpp
% vista/vistaprimitive.h
% vista/vistaprimitive.cpp
% vista/write_vista_vol.m

if nargin<1
   force=false;
end

% Possible COMPUTER types
% GLNX86
% GLNXA64
% PCWIN
% PCWIN64
% MAC
% MACI
% MACI64

% At the moment this works only for Linux 32-64 bit
% Vista folder has to be in the same folder as sb_compile_vista

% Step 1
% compile the libvista.a shared library for UNIX/LINUX
L = [];
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dadd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'drot.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsymv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtrsm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ismax.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sgemv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sspr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'stpmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dasum.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'drotg.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsyr2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtrsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ismin.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sger.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssum.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'stpsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'daxpy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsbmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsyr2k.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dvcal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'isum.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'snrm2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sswap.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'strmm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dcopy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dscal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsyr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'iadd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'blas'],'lsame.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'spmpar.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssymm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'strmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ddot.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dset.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsyrk.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'iasum.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sadd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'srot.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssymv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'strsm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dgbmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dspmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtbmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'icopy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sasum.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'srotg.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssyr2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'strsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dgemm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dspr2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtbsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'idamax.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'saxpy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssbmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssyr2k.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'svcal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dgemv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dspr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtpmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'idmax.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'scopy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sscal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssyr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'blas'],'xerbla.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dger.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsum.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtpsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'idmin.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sdot.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sset.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'ssyrk.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dnrm2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dswap.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtrmm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'isamax.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sgbmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sspmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'stbmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dpmpar.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dsymm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'dtrmv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'iset.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sgemm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'sspr2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'blas'],'stbsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'elknda.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofini.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofu02.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofu12.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofu23.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'getadj.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'getmet.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'itpdat.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofagn.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofu01.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofu03.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofu13.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'fofulb.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'getjac.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'itpanz.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'fofu'],'loknko.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'axunsu.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'axunsy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'cgstab.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'esv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'freund.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'illcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'kograd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'matpro.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'partch.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'pccgsy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'prcgst.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'precnd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'ptfqmr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'ransol.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'scalen.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'scalun.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'solchs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'solvel.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'solveu.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'solvlu.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'unscal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'libsol'],'unscun.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'femutl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem11.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem1.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem5.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem9.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'findobfl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem12.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem6.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'newelement.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'linpac.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem13.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem3.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem7.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'numutl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem10.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem4.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'neurofem8.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'neurofem'],'readascii.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'libtest.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfbase.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qidbl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qsinpu.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcctof.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfcopy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qinte.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qsmnam.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcdtoa.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfhist.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qiprot.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qsstam.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcfill.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfilao.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qireal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qssysv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcftoc.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfmove.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qitori.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qsvlen.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcitoa.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfname.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qiyeno.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qsznam.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qclen.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfpath.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qlbini.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qtpath.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcllen.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfprot.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qlbver.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qvfind.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcrtoa.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfreeu.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qlisin.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qvread.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qctolo.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qftemp.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qsdefi.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qvset.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qctoup.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qfunit.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qsdefo.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qvshow.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qcvgl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qichar.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'neutral'],'qshist.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qdtext.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfilda.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfredi.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qsdati.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qsprnt.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfclal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qflegn.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfsize.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qsecht.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qsuser.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfclos.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qflegp.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfwrdi.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qshell.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfcrdi.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfopdi.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
% L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qsbzzz.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qspara.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfcrse.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfopse.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qsclr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qspdbl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfdel.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qfrecl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qscput.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'ime' filesep 'linux'],'qsprea.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dbdsqr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeqrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlacpy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlapy3.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlauu2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbtf2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dspcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsytrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgerfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dladiv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaqgb.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlauum.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbtrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dspev.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtbcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbequ.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgerq2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlae2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaqge.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlazro.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbtrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dspevx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtbrfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbrfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgerqf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaebz.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaqsb.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dopgtr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpocon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dspgst.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtbtrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgesvd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaein.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaqsp.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dopmtr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpoequ.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dspgv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtgevc.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbsvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgesv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaev2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaqsy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorg2l.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dporfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsprfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtgsja.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbtf2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgesvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaexc.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaqtr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorg2r.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dposv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dspsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtpcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbtrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgetf2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlag2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlar2v.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorgbr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dposvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dspsvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtprfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgbtrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgetrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlags2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlarfb.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorghr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpotf2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsptrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtptri.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgebak.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgetri.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlagtf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlarf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorgl2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpotrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsptrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtptrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgebal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgetrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlagtm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlarfg.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorglq.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpotri.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsptri.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgebd2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dggbak.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlagts.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlarft.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorgql.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpotrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsptrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrevc.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgebrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dggbal.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlahqr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlarfx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorgqr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dppcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dstebz.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrexc.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgecon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dggglm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlahrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlargv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorgr2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dppequ.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dstein.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrrfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeequ.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgghrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaic1.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlarnv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorgrq.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpprfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsteqr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrsen.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgees.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgglse.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaln2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlartg.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorgtr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dppsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsterf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrsna.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeesx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dggqrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlamch.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlartv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorm2l.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dppsvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dstev.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrsyl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeev.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dggrqf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlangb.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaruv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorm2r.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpptrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dstevx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrti2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeevx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dggsvd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlange.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlas2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormbr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpptri.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsycon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrtri.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgegs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dggsvp.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlangt.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlascl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormhr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpptrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsyev.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtrtrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgegv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgtcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlanhs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaset.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dorml2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dptcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsyevx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dtzrqf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgehd2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgtrfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlansb.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlasr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormlq.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpteqr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsygs2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'icmax1.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgehrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgtsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlansp.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlassq.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormql.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dptrfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsygst.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'ilaenv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgelq2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgtsvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlanst.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlasv2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormqr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dptsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsygv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'lsame.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgelqf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgttrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlansy.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlaswp.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormr2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dptsvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsyrfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'lsamen.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgels.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgttrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlantb.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlasy2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormrq.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpttrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsysv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'xerbla.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgelss.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dhgeqz.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlantp.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlasyf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dormtr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpttrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsysvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgelsx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dhsein.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlantr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlatbs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbcon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'drscl.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsytd2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeql2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dhseqr.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlanv2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlatps.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbequ.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsbev.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsytf2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeqlf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlabad.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlapll.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlatrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbrfs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsbevx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsytrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeqpf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlabrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlapmt.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlatrs.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbsv.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsbtrd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsytrf.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dgeqr2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlacon.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlapy2.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dlatzm.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dpbsvx.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsecnd.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');
L = add_mex_source(L,['simbio-src' filesep 'lapack'],'dsytri.f',{'GLNX86', 'GLNXA64','PCWIN','PCWIN64','MACI','MACI64'},[],'-c');

oldDir = pwd;
[baseDir, myName] = fileparts(mfilename('fullpath'));
try
  compile_mex_list_sb(L, baseDir, force); 
catch
  me = lasterror;
  cd(oldDir);
  rethrow(me);
end
cd(baseDir);
% FIXME: the following should be done inside a function called add_shared_library
if isunix || ismac
 system('ar rcs simbio-src/libblas_sb.a simbio-src/blas/*.o');
 system('ar rcs simbio-src/liblapack_sb.a simbio-src/lapack/*.o');
 system('ar rcs simbio-src/libime_sb.a simbio-src/ime/neutral/*.o simbio-src/ime/linux/*.o');
 system('ar rcs simbio-src/libneurofem_sb.a simbio-src/neurofem/*.o simbio-src/fofu/*.o simbio-src/libsol/*.o');
 % Step 2
 % compile the mesh/vol functions
 L = [];
 L = add_mex_source(L,'.','calc_stiff_matrix_val.F',{'GLNX86','GLNXA64','MACI','MACI64'},[],'-Isimbio-src/neurofem -Lsimbio-src -lneurofem_sb -llapack_sb -lblas_sb -lime_sb -fortran');
elseif ispc
 % Step 2
 % compile the mesh/vol functions
 L = [];
 L = add_mex_source(L,'.','calc_stiff_matrix_val.F',{'PCWIN','PCWIN64'},[],'-Isimbio-src\neurofem simbio-src\blas\*.obj simbio-src\fofu\*.obj simbio-src\ime\neutral\*.obj simbio-src\ime\linux\*.obj simbio-src\lapack\*.obj simbio-src\libsol\*.obj simbio-src\neurofem\*.obj');
end


try
  compile_mex_list_sb(L, baseDir, force); 
catch
  me = lasterror;
  cd(oldDir);
  rethrow(me);
end
cd(oldDir);

function compile_mex_list_sb(L, baseDir, force)
% function compile_mex_list(L, baseDir)
%
% Compile a list of MEX files as determined by the input argument L.
% The second argument 'baseDir' is the common base directory for the
% files listed in L. The third argument is a flag that determines
% whether to force (re-)compilation even if the MEX file is up-to-date.
%
% See also ft_compile_mex, add_mex_source.

% (C) 2010 S. Klanke, 2011, Cristiano Micheli

for i=1:length(L)
   [~, name] = fileparts(L(i).relName);

   sfname = [baseDir filesep L(i).dir filesep L(i).relName ];
   SF = dir(sfname);
   if numel(SF)<1
      fprintf(1,'Error: source file %s cannot be found.\n', sfname);
      continue;
   end
   
   if ~force
      if(~isempty(strfind(L(i).extras,'-c')))
          if isunix || ismac
              mfname = [baseDir filesep L(i).dir filesep name '.o'];
              MF = dir(mfname);
          elseif ispc
              mfname = [baseDir filesep L(i).dir filesep name '.obj'];
              MF = dir(mfname);
          end
      else
          mfname = [baseDir filesep L(i).dir filesep name '.' mexext];
          MF = dir(mfname);
      end
      if numel(MF)==1 && datenum(SF.date) <= datenum(MF.date)
            fprintf(1,'Skipping up-to-date file %s%s%s ...\n', L(i).dir, filesep, name);
         continue;
      end 
   end
   fprintf(1,'Compiling file %s%s%s ...\n', L(i).dir, filesep, name);
   cd([baseDir filesep L(i).dir]);
   cmd = sprintf('mex %s %s', L(i). relName, L(i).extras);
   eval(cmd);
end

function L = add_mex_source(L, directory, relName, matchPlatform, excludePlatform, extras)
% function L = add_mex_source(L, directory, relName, matchPlatform, excludePlatform, extras)
%
% Input + output argument L is a structure array of directory names, source file names,
% and extra arguments required for the compilation of MEX files. This function will
% create a new element of this structure and append it to L.
%
% Further inputs:
%   directory
%      target directory of the mex-file
%   relName
%      source file relative to 'directory'
%   matchPlatform
%      list of platforms this MEX file should only be compiled for.
%      use an empty matrix [] to compile for all platforms
%   excludePlatform
%      list of platforms this MEX file should NOT be compiled for.
%   extras
%      extra arguments to the MEX command, e.g. additional source files

% (C) 2010 S. Klanke

% Check if this file only needs compilation on certain platforms (including this one)
if nargin>3 && ~isempty(matchPlatform) 
   ok = false;
   for k=1:numel(matchPlatform)
       if strcmp(matchPlatform{k}, computer)
		  ok = true;
		  break;
	   end
	end
	if ~ok
	   return
	end
end

% Check if this file cannot be compiled on certain platforms (including this one)
if nargin>4 && ~isempty(excludePlatform) 
   ok = true;
   for k=1:numel(excludePlatform)
      if strcmp(excludePlatform{k}, computer)
         return;
      end
   end
end

L(end+1).dir   = directory;
L(end).relName = relName;
if nargin>5
  L(end).extras = extras;
end
