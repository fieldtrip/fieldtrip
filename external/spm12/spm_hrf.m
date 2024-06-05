function [hrf,p] = spm_hrf(RT,P,T)
% Haemodynamic response function
% FORMAT [hrf,p] = spm_hrf(RT,p,T)
% RT   - scan repeat time
% p    - parameters of the response function (two Gamma functions)
%
%                                                           defaults
%                                                          {seconds}
%        p(1) - delay of response (relative to onset)          6
%        p(2) - delay of undershoot (relative to onset)       16
%        p(3) - dispersion of response                         1
%        p(4) - dispersion of undershoot                       1
%        p(5) - ratio of response to undershoot                6
%        p(6) - onset {seconds}                                0
%        p(7) - length of kernel {seconds}                    32
%
% T    - microtime resolution [Default: 16]
%
% hrf  - haemodynamic response function
% p    - parameters of the response function
%__________________________________________________________________________
%
% The parameters p(1:4) correspond to the shape and scale parameters of two
% probability density functions of the Gamma distribution (see spm_Gpdf.m),
% one corresponding to the main response and the other one to the
% undershoot.
% Note that the mean of the Gamma distribution is shape*scale and its mode
% is (shape-1)*scale.  This means that with the default values of the
% parameters the peak of the heamodynamic response function will be around
% 5 seconds.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging


%-Parameters of the response function
%--------------------------------------------------------------------------
try
    p = spm_get_defaults('stats.fmri.hrf');
catch
    p = [6 16 1 1 6 0 32];
end
if nargin > 1
    p(1:length(P)) = P;
end

%-Microtime resolution
%--------------------------------------------------------------------------
if nargin > 2
    fMRI_T = T;
else
    fMRI_T = spm_get_defaults('stats.fmri.t');
end

%-Modelled haemodynamic response function - {mixture of Gammas}
%--------------------------------------------------------------------------
dt  = RT/fMRI_T;
u   = [0:ceil(p(7)/dt)] - p(6)/dt;
hrf = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
hrf = hrf([0:floor(p(7)/RT)]*fMRI_T + 1);
hrf = hrf'/sum(hrf);
