function [Phi,DPhi] = spm_dartel_integrate(U,t,K)
% Integrate a Dartel flow field
% FORMAT [Phi,DPhi] = spm_dartel_exp(U,t,K)
%     U    - name of flow field (nx x ny x nz x nt x 3)
%     t    - [t0 t1] Start and end time (values between 0 and 1)
%     K    - log2 of the Euler time steps to integrate the
%            flow field.
%
%     Phi  - deformation field (nx x ny x nz x 3)
%     DPhi - Jacobian determinant field (nx x ny x nz)
%
% The function integrates
%     Phi(x,t) = \int_{t_0}^{t_1} U(Phi(x,t),t) dt
% where U is a piecewise constant flow field
%
% Note: this function is ready for LDDMM-style flow fields, even
% though the none of the official Dartel tools can generate them
% yet.
% _______________________________________________________________________
%  Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_integrate.m 5506 2013-05-14 17:13:43Z john $

if isa(U,'char'), U = nifti(U); end;
if isa(U,'nifti'), U = U.dat; end;

% Figure out which bits of flow field to use, the numbers of
% Euler time steps, and the scales.
nt  = size(U,4);
K   = ceil(log2(2^K/nt));
big = 2^12;
t   = min(max(t,0),1);
t0  = t(1);
t1  = t(2);
tt0 = t0*nt;
tt1 = t1*nt;
sc  = zeros(1,nt);
if t1>t0,
    for i=1:nt,
        sc(i) = max(round((min(tt1,i)-max(tt0,i-1))*big)/big,0);
    end
    ind = find(sc~=0);
else
    for i=1:nt,
        sc(i) = max(round((min(tt0,i)-max(tt1,i-1))*big)/big,0);
    end
    ind = fliplr(find(sc~=0));
end
if isempty(ind), % t0==t1
    ind = 1;
    sc  = 0;
    ts  = 0;
else
    sc = sc(ind);
    ts = zeros(1,numel(ind));
    for i=1:numel(ind),
        ts(i) = ceil(log2((2^K)*sc(i))-1/big);
        sc(i) = sc(i)*2^(K-ts(i));
    end
    if t0>t1, sc = -sc; end
end


% Do the integrations
if nargout==1,
    % Deformation field only
    u   = squeeze(single(U(:,:,:,ind(1),:)));
    Phi = dartel3('Exp',u,[ts(1), sc(1)]);
    for i=2:numel(ind),
        u   = squeeze(single(U(:,:,:,ind(i),:)));
        Phi = dartel3('comp',Phi,dartel3('Exp',u,[ts(i), sc(i)]));
    end
else
    % Deformation and Jacobian determinant fields
    u          = squeeze(single(U(:,:,:,ind(1),:)));
    [Phi,DPhi] = dartel3('Exp',u,[ts(1), sc(1), 1]);
    for i=2:numel(ind),
        u            = squeeze(single(U(:,:,:,ind(i),:)));
        [Phi1,DPhi1] = dartel3('Exp',u,[ts(i), sc(i), 1]);
        [Phi,DPhi]   = dartel3('comp',Phi,Phi1,DPhi,DPhi1);
        clear Phi1 DPhi1
    end
end
