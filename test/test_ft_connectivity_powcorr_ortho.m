function test_ft_connectivity_powcorr_ortho

% MEM 1500mb
% WALLTIME 00:10:00

% TEST: test_ft_connectivity_powcorr_ortho
% TEST: ft_connectivity_powcorr_ortho

mom1 = randn(1,100)+1i*randn(1,100);
mom2 = randn(1,100)+1i*randn(1,100);

c = ft_connectivity_powcorr_ortho([mom1;mom2], 'refindx', 1);
c = c(2,:);

[c1, c2] = hipp_testfunction(mom1, mom2);
%assert(all(abs(c-[c1 c2])<10*eps));
assert(all(abs(c-(c1+c2)./2)<10*eps));

% subfunction that does the computation according to the paper
% Nat Neuro 2012 Hipp et al.
function [c1, c2] = hipp_testfunction(mom1, mom2)

% normalise the amplitudes
mom1n = mom1./abs(mom1);
mom2n = mom2./abs(mom2);

% rotate 
mom12 = mom1.*conj(mom2n);
mom21 = mom2.*conj(mom1n);

% take the projection along the imaginary axis
mom12i = abs(imag(mom12));
mom21i = abs(imag(mom21));

% compute the correlation on the log transformed power values
c1 = corr(log10(mom12i.^2'), log10(abs(mom2).^2')); 
c2 = corr(log10(mom21i.^2'), log10(abs(mom1).^2'));
