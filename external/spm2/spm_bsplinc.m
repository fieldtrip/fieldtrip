function varargout = spm_bsplinc(varargin)
% Return volume of B-spline coefficients
% FORMAT c = spm_bsplinc(V,d)
% 	V - a structure of image volume information (from spm_vol)
% 	d(1:3) - degree of B-spline (from 0 to 7) along different dimensions
%       d(4:6) - 1/0 to indicate wrapping along the dimensions
% 	c - returned volume of B-spline coefficients
%
% This function deconvolves B-splines from volume V, returning
% coefficients, c.  These coefficients are then passed to spm_bsplins
% in order to sample the data using B-spline interpolation.
%
%_______________________________________________________________________
%
% References:
%	M. Unser, A. Aldroubi and M. Eden.
%	"B-Spline Signal Processing: Part I-Theory,"
%	IEEE Transactions on Signal Processing 41(2):821-832 (1993).
%
%	M. Unser, A. Aldroubi and M. Eden.
%	"B-Spline Signal Processing: Part II-Efficient Design and
%	Applications,"
%	IEEE Transactions on Signal Processing 41(2):834-848 (1993).
%
%	M. Unser.
%	"Splines: A Perfect Fit for Signal and Image Processing,"
%	IEEE Signal Processing Magazine, 16(6):22-38 (1999)
%
%	P. Thevenaz and T. Blu and M. Unser.
%	"Interpolation Revisited"
%	IEEE Transactions on Medical Imaging 19(7):739-758 (2000).
%_______________________________________________________________________
% @(#)spm_bsplinc.m	2.3 John Ashburner 02/07/30

%-This is merely the help file for the compiled routine
error('spm_bsplinc.c not compiled.');

