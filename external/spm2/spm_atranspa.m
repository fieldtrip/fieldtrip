function varargout = spm_atranspa(varargin)
% Multiplies the transpose of a matrix by itself - a compiled routine
% FORMAT C = spm_atranspa(A)
% A - real matrix
% C - real symmetric matrix resulting from A'A
%_______________________________________________________________________
%
% This routine was written to save both memory and CPU time.
% The memory saving is achieved by not having to generate A'.
% CPU saving is by only generating half of C, and filling the
% rest in later.
%_______________________________________________________________________
% @(#)spm_atranspa.m	2.1 John Ashburner MRCCU/FIL 99/04/19

%-This is merely the help file for the compiled routine
error('spm_atranspa.c not compiled - see spm_MAKE.sh')
