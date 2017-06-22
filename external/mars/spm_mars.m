function deltaFEF = spm_mars(J,q,l,ext_eng,lkp0,rLabel)
% The core updating step (E-step) for the posterior q_ik - a compiled routine
% FORMAT deltaFEF = spm_mars(J,q,l,ext_eng,lkp0,rLabel)
%
% See spm_mars_core.m for details on inputs and output.
%
% This function implements the E-step updating for the segmentation. The
% corresponding equation is the Eq. 2 in the following paper:
%
% Huang Y, Parra LC (2015) Fully Automated Whole-Head Segmentation with Improved
% Smoothness and Continuity, with Theory Reviewed. PLoS ONE 10(5): e0125477.
% doi:10.1371/journal.pone.0125477
%
%
% Yu (Andy) Huang
% $Id: spm_mars.m 2015-07-27 andy$
% Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% yhuang16@citymail.cuny.edu
%
%
% This is only the help file for the compiled routine
error('spm_mars.c not compiled - see Makefile');