function T = spm_bias_estimate(V,flags)
% Estimate image nonuniformity.
%
% FORMAT T = spm_bias_estimate(V,flags)
%   V     - filename or vol struct of image
%   flags - a structure containing the following fields
%     nbins  - number of bins in histogram (1024)
%     reg    - amount of regularisation (1)
%     cutoff - cutoff (in mm) of basis functions (35)
%   T     - DCT of bias field.
%
% The objective function is related to minimising the entropy of
% the image histogram, but is modified slightly.
% This fixes the problem with the SPM99 non-uniformity correction
% algorithm, which tends to try to reduce the image intensities. As
% the field was constrainded to have an average value of one, then
% this caused the field to bend upwards in regions not included in
% computations of image non-uniformity.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_bias_estimate.m 1143 2008-02-07 19:33:33Z spm $


def_flags = struct('nbins',256,'reg',0.01,'cutoff',30,'niter',128);
if nargin < 2,
        flags = def_flags;
else
        fnms = fieldnames(def_flags);
        for i=1:length(fnms),
                if ~isfield(flags,fnms{i}),
                    flags.(fnms{i}) = def_flags.(fnms{i});
                end;
        end;
end;

reg     = flags.reg;     % Regularisation
co      = flags.cutoff;  % Highest wavelength of DCT
nh      = flags.nbins;
niter   = flags.niter;

if ischar(V), V = spm_vol(V); end;
mx      = 1.1*get_max(V); % Maximum value in histogram

tmp     = sqrt(sum(V(1).mat(1:3,1:3).^2));
nbas    = max(round((V(1).dim(1:3).*tmp)/co),[1 1 1]);
B1      = spm_dctmtx(V(1).dim(1),nbas(1));
B2      = spm_dctmtx(V(1).dim(2),nbas(2));
B3      = spm_dctmtx(V(1).dim(3),nbas(3));

[T,IC0] = get_priors(V,nbas,reg,4);
IC0     = IC0(2:end,2:end);
%tmp     = diag(IC0);
%tmp     = tmp(tmp>0);
%offset = 0.5*(log(2*pi)*length(tmp) + sum(log(tmp)));
offset  = 0;
fprintf('*** %s ***\nmax = %g, Bases = %dx%dx%d\n', V.fname, mx, nbas);

olpp     = Inf;
x        = (0:(nh-1))'*(mx/(nh-1));
plt      = zeros(64,2);

for iter = 1:flags.niter,

    [Alpha,Beta,ll, h, n] = spm_bias_mex(V,B1,B2,B3,T,[mx nh]);
    T     = T(:);
    T     = T(2:end);
    Alpha = Alpha(2:end,2:end)/n;
    Beta  = Beta(2:end)/n;
    ll    = ll/n;
    lp    = offset + 0.5*(T'*IC0*T);
    lpp   = lp+ll;
    T     = (Alpha + IC0)\(Alpha*T - Beta);
    T     = reshape([0 ; T],nbas);

    [pth,nm] = fileparts(deblank(V.fname));
    S        = fullfile(pth,['bias_' nm '.mat']);
    %S       = ['bias_' nm '.mat'];
    save(S,'V','T','h');
    fprintf('%g %g\n', ll, lp);

    if iter==1,
        fg = spm_figure('FindWin','Interactive');
        if ~isempty(fg),
            spm_figure('Clear',fg);
            ax1 = axes('Position', [0.15 0.1 0.8 0.35],...
                'Box', 'on','Parent',fg);
            ax2 = axes('Position', [0.15 0.6 0.8 0.35],...
                'Box', 'on','Parent',fg);
        end;
        h0 = h;
    end;
    if ~isempty(fg),
        plt(iter,1) = ll;
        plt(iter,2) = lpp;
        plot((1:iter)',plt(1:iter,:),'Parent',ax1,'LineWidth',2);
        set(get(ax1,'Xlabel'),'string','Iteration','FontSize',10);
        set(get(ax1,'Ylabel'),'string','Negative Log-Likelihood','FontSize',10);
        set(ax1,'Xlim',[1 iter+1]);

        plot(x,h0,'r', x,h,'b', 'Parent',ax2);
        set(ax2,'Xlim',[0 x(end)]);
        set(get(ax2,'Xlabel'),'string','Intensity','FontSize',10);
        set(get(ax2,'Ylabel'),'string','Frequency','FontSize',10);
        drawnow;
    end;

    if olpp-lpp < 1e-5 && ~isempty(fg), delete([ax1 ax2]); break; end;
    olpp = lpp;
end;
return;
%=======================================================================

%=======================================================================
function mx = get_max(V)
mx = 0;
for i=1:V.dim(3),
    img = spm_slice_vol(V,spm_matrix([0 0 i]),V.dim(1:2),0);
    mx  = max([mx max(img(:))]);
end;
return;
%=======================================================================

%=======================================================================
function [T,IC0] = get_priors(VF,nbas,reg,o)
% Set up a priori covariance matrix
vx    = sqrt(sum(VF(1).mat(1:3,1:3).^2));
kx    = (((1:nbas(1))'-1)*pi/vx(1)/VF(1).dim(1)*10).^2;
ky    = (((1:nbas(2))'-1)*pi/vx(2)/VF(1).dim(2)*10).^2;
kz    = (((1:nbas(3))'-1)*pi/vx(3)/VF(1).dim(3)*10).^2;

switch o,
case 0, % Cost function based on sum of squares
IC0 =    kron(kz.^0,kron(ky.^0,kx.^0))*reg;

case 1, % Cost function based on sum of squared 1st derivatives
IC0 = (  kron(kz.^1,kron(ky.^0,kx.^0)) +...
         kron(kz.^0,kron(ky.^1,kx.^0)) +...
         kron(kz.^0,kron(ky.^0,kx.^1)) )*reg;

case 2, % Cost function based on sum of squared 2nd derivatives
IC0 = (1*kron(kz.^2,kron(ky.^0,kx.^0)) +...
       1*kron(kz.^0,kron(ky.^2,kx.^0)) +...
       1*kron(kz.^0,kron(ky.^0,kx.^2)) +...
       2*kron(kz.^1,kron(ky.^1,kx.^0)) +...
       2*kron(kz.^1,kron(ky.^0,kx.^1)) +...
       2*kron(kz.^0,kron(ky.^1,kx.^1)) )*reg;

case 3, % Cost function based on sum of squared 3rd derivatives
IC0 = (1*kron(kz.^3,kron(ky.^0,kx.^0)) +...
       1*kron(kz.^0,kron(ky.^3,kx.^0)) +...
       1*kron(kz.^0,kron(ky.^0,kx.^3)) +...
       3*kron(kz.^2,kron(ky.^1,kx.^0)) +...
       3*kron(kz.^2,kron(ky.^0,kx.^1)) +...
       3*kron(kz.^1,kron(ky.^2,kx.^0)) +...
       3*kron(kz.^0,kron(ky.^2,kx.^1)) +...
       3*kron(kz.^1,kron(ky.^0,kx.^2)) +...
       3*kron(kz.^0,kron(ky.^1,kx.^2)) +...
       6*kron(kz.^1,kron(ky.^1,kx.^1)) )*reg;

case 4, % Cost function based on sum of squares of 4th derivatives
IC0 =  (1*kron(kz.^4,kron(ky.^0,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^4,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^0,kx.^4)) +...
        4*kron(kz.^3,kron(ky.^1,kx.^0)) +...
        4*kron(kz.^3,kron(ky.^0,kx.^1)) +...
        4*kron(kz.^1,kron(ky.^3,kx.^0)) +...
        4*kron(kz.^0,kron(ky.^3,kx.^1)) +...
        4*kron(kz.^1,kron(ky.^0,kx.^3)) +...
        4*kron(kz.^0,kron(ky.^1,kx.^3)) +...
        6*kron(kz.^2,kron(ky.^2,kx.^0)) +...
        6*kron(kz.^2,kron(ky.^0,kx.^2)) +...
        6*kron(kz.^0,kron(ky.^2,kx.^2)) +...
       12*kron(kz.^2,kron(ky.^1,kx.^1)) +...
       12*kron(kz.^1,kron(ky.^2,kx.^1)) +...
       12*kron(kz.^1,kron(ky.^1,kx.^2)) )*reg;
otherwise,
error('Unknown regularisation form');
end;

IC0(1) = max([max(IC0) 1e4]);
IC0    = diag(IC0);

% Initial estimate for intensity modulation field
T      = zeros(nbas(1),nbas(2),nbas(3),1);
return;
%=======================================================================
