function x = spm_coreg(varargin)
% Between modality coregistration using information theory
% FORMAT x = spm_coreg(VG,VF,flags)
% VG    - handle for reference image (see spm_vol).
% VF    - handle for source (moved) image.
% flags - a structure containing the following elements:
%          sep      - optimisation sampling steps (mm)
%                     default: [4 2]
%          params   - starting estimates (6 elements)
%                     default: [0 0 0  0 0 0]
%          cost_fun - cost function string:
%                       'mi'  - Mutual Information
%                       'nmi' - Normalised Mutual Information
%                       'ecc' - Entropy Correlation Coefficient
%                       'ncc' - Normalised Cross Correlation
%                     default: 'nmi'
%          tol      - tolerences for accuracy of each param
%                     default: [0.02 0.02 0.02 0.001 0.001 0.001]
%          fwhm     - smoothing to apply to 256x256 joint histogram
%                     default: [7 7]
%          graphics - display coregistration outputs
%                     default: ~spm('CmdLine')
%
% x     - the parameters describing the rigid body rotation, such that a
%         mapping from voxels in G to voxels in F is attained by:
%         VF.mat\spm_matrix(x(:)')*VG.mat
%
% At the end, the voxel-to-voxel affine transformation matrix is
% displayed, along with the histograms for the images in the original
% orientations, and the final orientations. The registered images are
% displayed at the bottom.
%__________________________________________________________________________
%
% The registration method used here is based on the work described in:
% A Collignon, F Maes, D Delaere, D Vandermeulen, P Suetens & G Marchal
% (1995) "Automated Multi-modality Image Registration Based On
% Information Theory". In the proceedings of Information Processing in
% Medical Imaging (1995).  Y. Bizais et al. (eds.).  Kluwer Academic
% Publishers.
%
% The original interpolation method described in this paper has been
% changed in order to give a smoother cost function.  The images are
% also smoothed slightly, as is the histogram.  This is all in order to
% make the cost function as smooth as possible, to give faster convergence
% and less chance of local minima.
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_coreg.m 6435 2015-05-14 09:59:54Z guillaume $

%--------------------------------------------------------------------------
% References
%==========================================================================
%
% Mutual Information
% -------------------------------------------------------------------------
% Collignon, Maes, Delaere, Vandermeulen, Suetens & Marchal (1995).
% "Automated multi-modality image registration based on information theory".
% In Bizais, Barillot & Di Paola, editors, Proc. Information Processing
% in Medical Imaging, pages 263--274, Dordrecht, The Netherlands, 1995.
% Kluwer Academic Publishers.
%
% Wells III, Viola, Atsumi, Nakajima & Kikinis (1996).
% "Multi-modal volume registration by maximisation of mutual information".
% Medical Image Analysis, 1(1):35-51, 1996. 
%
% Entropy Correlation Coefficient
% -------------------------------------------------------------------------
% Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
% "Multimodality image registration by maximisation of mutual
% information". IEEE Transactions on Medical Imaging 16(2):187-198
%
% Normalised Mutual Information
% -------------------------------------------------------------------------
% Studholme, Hill & Hawkes (1998).
% "A normalized entropy measure of 3-D medical image alignment".
% in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.             
%
% Optimisation
% -------------------------------------------------------------------------
% Press, Teukolsky, Vetterling & Flannery (1992).
% "Numerical Recipes in C (Second Edition)".
% Published by Cambridge.
%--------------------------------------------------------------------------

SVNid = '$Rev: 6435 $';

if nargin >= 4
    x = optfun(varargin{:});
    return;
end

%-Say hello
%--------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SVNid);

def_flags          = spm_get_defaults('coreg.estimate');
def_flags.params   = [0 0 0  0 0 0];
def_flags.graphics = ~spm('CmdLine');
if nargin < 3
    flags = def_flags;
else
    flags = varargin{3};
    fnms  = fieldnames(def_flags);
    for i=1:length(fnms)
        if ~isfield(flags,fnms{i})
            flags.(fnms{i}) = def_flags.(fnms{i});
        end
    end
end

if nargin < 1
    VG = spm_vol(spm_select(1,'image','Select reference image'));
else
    VG = varargin{1};
    if ischar(VG), VG = spm_vol(VG); end
end
if nargin < 2
    VF = spm_vol(spm_select(Inf,'image','Select moved image(s)'));
else
    VF = varargin{2};
    if ischar(VF) || iscellstr(VF), VF = spm_vol(char(VF)); end;
end

if ~isfield(VG, 'uint8')
    VG.uint8 = loaduint8(VG);
    vxg      = sqrt(sum(VG.mat(1:3,1:3).^2));
    fwhmg    = sqrt(max([1 1 1]*flags.sep(end)^2 - vxg.^2, [0 0 0]))./vxg;
    VG       = smooth_uint8(VG,fwhmg); % Note side effects
end

sc = flags.tol(:)'; % Required accuracy
sc = sc(1:length(flags.params));
xi = diag(sc*20);
x = zeros(numel(VF),numel(flags.params));

for k=1:numel(VF)
    VFk = VF(k);
    if ~isfield(VFk, 'uint8')
        VFk.uint8 = loaduint8(VFk);
        vxf       = sqrt(sum(VFk.mat(1:3,1:3).^2));
        fwhmf     = sqrt(max([1 1 1]*flags.sep(end)^2 - vxf.^2, [0 0 0]))./vxf;
        VFk       = smooth_uint8(VFk,fwhmf); % Note side effects
    end

    xk  = flags.params(:);
    for samp=flags.sep(:)'
        xk     = spm_powell(xk(:), xi,sc,mfilename,VG,VFk,samp,flags.cost_fun,flags.fwhm);
        x(k,:) = xk(:)';
    end
    if flags.graphics
        display_results(VG(1),VFk(1),xk(:)',flags);
    end
end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#

%==========================================================================
% function o = optfun(x,VG,VF,s,cf,fwhm)
%==========================================================================
function o = optfun(x,VG,VF,s,cf,fwhm)
% The function that is minimised.
if nargin<6, fwhm = [7 7];   end
if nargin<5, cf   = 'mi';    end
if nargin<4, s    = [1 1 1]; end

% Voxel sizes
vxg = sqrt(sum(VG.mat(1:3,1:3).^2));sg = s./vxg;

% Create the joint histogram
H = spm_hist2(VG.uint8,VF.uint8, VF.mat\spm_matrix(x(:)')*VG.mat ,sg);

% Smooth the histogram
lim  = ceil(2*fwhm);
krn1 = spm_smoothkern(fwhm(1),-lim(1):lim(1)) ; krn1 = krn1/sum(krn1); H = conv2(H,krn1);
krn2 = spm_smoothkern(fwhm(2),-lim(2):lim(2))'; krn2 = krn2/sum(krn2); H = conv2(H,krn2);

% Compute cost function from histogram
H  = H+eps;
sh = sum(H(:));
H  = H/sh;
s1 = sum(H,1);
s2 = sum(H,2);

switch lower(cf)
    case 'mi'
        % Mutual Information:
        H   = H.*log2(H./(s2*s1));
        mi  = sum(H(:));
        o   = -mi;
    case 'ecc'
        % Entropy Correlation Coefficient of:
        % Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
        % "Multimodality image registration by maximisation of mutual
        % information". IEEE Transactions on Medical Imaging 16(2):187-198
        H   = H.*log2(H./(s2*s1));
        mi  = sum(H(:));
        ecc = -2*mi/(sum(s1.*log2(s1))+sum(s2.*log2(s2)));
        o   = -ecc;
    case 'nmi'
        % Normalised Mutual Information of:
        % Studholme,  Hill & Hawkes (1998).
        % "A normalized entropy measure of 3-D medical image alignment".
        % in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
        nmi = (sum(s1.*log2(s1))+sum(s2.*log2(s2)))/sum(sum(H.*log2(H)));
        o   = -nmi;
    case 'ncc'
        % Normalised Cross Correlation
        i     = 1:size(H,1);
        j     = 1:size(H,2);
        m1    = sum(s2.*i');
        m2    = sum(s1.*j);
        sig1  = sqrt(sum(s2.*(i'-m1).^2));
        sig2  = sqrt(sum(s1.*(j -m2).^2));
        [i,j] = ndgrid(i-m1,j-m2);
        ncc   = sum(sum(H.*i.*j))/(sig1*sig2);
        o     = -ncc;
    otherwise
        error('Invalid cost function specified');
end


%==========================================================================
% function udat = loaduint8(V)
%==========================================================================
function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 && V.pinfo(1) == 2
    mx = 255*V.pinfo(1) + V.pinfo(2);
    mn = V.pinfo(2);
else
    spm_progress_bar('Init',V.dim(3),...
        ['Computing max/min of ' spm_file(V.fname,'filename')],...
        'Planes complete');
    mx = -Inf; mn =  Inf;
    for p=1:V.dim(3)
        img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
        img = img(isfinite(img));
        mx  = max([max(img(:))+paccuracy(V,p) mx]);
        mn  = min([min(img(:)) mn]);
        spm_progress_bar('Set',p);
    end
end

% Another pass to find a maximum that allows a few hot-spots in the data.
spm_progress_bar('Init',V.dim(3),...
        ['2nd pass max/min of ' spm_file(V.fname,'filename')],...
        'Planes complete');
nh = 2048;
h  = zeros(nh,1);
for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    img = img(isfinite(img));
    img = round((img+((mx-mn)/(nh-1)-mn))*((nh-1)/(mx-mn)));
    h   = h + accumarray(img,1,[nh 1]);
    spm_progress_bar('Set',p);
end
tmp = [find(cumsum(h)/sum(h)>0.9999); nh];
mx  = (mn*nh-mx+tmp(1)*(mx-mn))/(nh-1);

% Load data from file indicated by V into an array of unsigned bytes.
spm_progress_bar('Init',V.dim(3),...
    ['Loading ' spm_file(V.fname,'filename')],...
    'Planes loaded');
udat = zeros(V.dim,'uint8');
st = rand('state'); % st = rng;
rand('state',100); % rng(100,'v5uniform'); % rng('defaults');
for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    acc = paccuracy(V,p);
    if acc==0
        udat(:,:,p) = uint8(max(min(round((img-mn)*(255/(mx-mn))),255),0));
    else
        % Add random numbers before rounding to reduce aliasing artifact
        r = rand(size(img))*acc;
        udat(:,:,p) = uint8(max(min(round((img+r-mn)*(255/(mx-mn))),255),0));
    end
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');
rand('state',st); % rng(st);


%==========================================================================
% function acc = paccuracy(V,p)
%==========================================================================
function acc = paccuracy(V,p)
if ~spm_type(V.dt(1),'intt')
    acc = 0;
else
    if size(V.pinfo,2)==1
        acc = abs(V.pinfo(1,1));
    else
        acc = abs(V.pinfo(1,p));
    end
end


%==========================================================================
% function V = smooth_uint8(V,fwhm)
%==========================================================================
function V = smooth_uint8(V,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(V.uint8,V.uint8,x,y,z,-[i j k]);


%==========================================================================
% function display_results(VG,VF,x,flags)
%==========================================================================
function display_results(VG,VF,x,flags)
fig = spm_figure('FindWin','Graphics');
if isempty(fig), return; end;
set(0,'CurrentFigure',fig);
spm_figure('Clear','Graphics');

%txt = 'Information Theoretic Coregistration';
switch lower(flags.cost_fun)
    case 'mi',  txt = 'Mutual Information Coregistration';
    case 'ecc', txt = 'Entropy Correlation Coefficient Registration';
    case 'nmi', txt = 'Normalised Mutual Information Coregistration';
    case 'ncc', txt = 'Normalised Cross Correlation';
    otherwise, error('Invalid cost function specified');
end

% Display text
%--------------------------------------------------------------------------
ax = axes('Position',[0.1 0.8 0.8 0.15],'Visible','off','Parent',fig);
text(0.5,0.7, txt,'FontSize',16,...
    'FontWeight','Bold','HorizontalAlignment','center','Parent',ax);

Q = inv(VF.mat\spm_matrix(x(:)')*VG.mat);
text(0,0.5, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),'Parent',ax);
text(0,0.3, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),'Parent',ax);
text(0,0.1, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),'Parent',ax);

% Display joint histograms
%--------------------------------------------------------------------------
ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\VG.mat,[1 1 1]);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
    'XTick',[],'YTick',[]);
title('Original Joint Histogram','Parent',ax);
xlabel(spm_file(VG.fname,'short22'),'Parent',ax,'Interpreter','none');
ylabel(spm_file(VF.fname,'short22'),'Parent',ax,'Interpreter','none');

H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\spm_matrix(x(:)')*VG.mat,[1 1 1]);
ax  = axes('Position',[0.6 0.5 0.35 0.3],'Visible','off','Parent',fig);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
    'XTick',[],'YTick',[]);
title('Final Joint Histogram','Parent',ax);
xlabel(spm_file(VG.fname,'short22'),'Parent',ax,'Interpreter','none');
ylabel(spm_file(VF.fname,'short22'),'Parent',ax,'Interpreter','none');

% Display ortho-views
%--------------------------------------------------------------------------
spm_orthviews('Reset');
     spm_orthviews('Image',VG,[0.01 0.01 .48 .49]);
h2 = spm_orthviews('Image',VF,[.51 0.01 .48 .49]);
global st
st.vols{h2}.premul = inv(spm_matrix(x(:)'));
spm_orthviews('Space');

spm_print;
