function [err, bf] = spm_bf (TR)
%
%   [err,] = spm_bf ()
%
%Purpose:
%
%
%
%Input Parameters:
%
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%
%
%
%Key Terms:
%
%More Info :
%
%
%
%
%     Author : Gang Chen
%     Date : Mon Oct 27 17:13:26 EST 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda MD 20892


%Define the function name for easy referencing
FuncName = 'spm_bf.m';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;


%	p(1) - delay of response (relative to onset)	   6        -gamd  in waver -GAM
%	p(2) - delay of undershoot (relative to onset)    16
%	p(3) - dispersion of response			   1
%	p(4) - dispersion of undershoot			   1   X
%	p(5) - ratio of response to undershoot		   6   X
%	p(6) - onset (seconds)				   0     -gamd
%	p(7) - length of kernel (seconds)		  32

      p      = [6 16 1 1 6 0 32];   %default parameters
	   [bf1 p]         = spm_hrf(TR,p);      %get gamma hrf from SPM function spm_hrf
	   dp     = 1;
	   p(6)   = p(6) + dp;    %for calculating time derivative
	   bf2    = (bf1 - spm_hrf(TR,p))/dp;
%		D      = (bf(:,1) - spm_hrf(TR,p))/dp;   %time derivative
		bf     = [bf1 bf2];     %combine the two

err = 0;
return;

