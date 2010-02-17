% Script to compile C code for EEG-BEM

% Christophe Phillips
% $Id: make.m 2718 2009-02-09 19:14:20Z vladimir $

mex -O bem_Cii_cog.c
mex -O bem_Cii_cst.c
mex -O bem_Cii_lin.c
mex -O bem_Cij_cog.c
mex -O bem_Cij_cst.c
mex -O bem_Cij_lin.c
mex -O bem_Gi_cog.c
mex -O bem_Gi_vert.c
