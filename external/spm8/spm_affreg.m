function [M,scal] = spm_affreg(VG,VF,flags,M,scal)
% Affine registration using least squares.
% FORMAT [M,scal] = spm_affreg(VG,VF,flags,M0,scal0)
%
% VG        - Vector of template volumes.
% VF        - Source volume.
% flags     - a structure containing various options.  The fields are:
%             WG       - Weighting volume for template image(s).
%             WF       - Weighting volume for source image
%                        Default to [].
%             sep      - Approximate spacing between sampled points (mm).
%                        Defaults to 5.
%             regtype  - regularisation type.  Options are:
%                        'none'  - no regularisation
%                        'rigid' - almost rigid body
%                        'subj'  - inter-subject registration (default).
%                        'mni'   - registration to ICBM templates
%             globnorm - Global normalisation flag (1)
% M0        - (optional) starting estimate. Defaults to eye(4).
% scal0     - (optional) starting estimate.
%
% M         - affine transform, such that voxels in VF map to those in
%             VG by   VG.mat\M*VF.mat
% scal      - scaling factors for VG
%
% When only one template is used, then the cost function is approximately
% symmetric, although a linear combination of templates can be used.
% Regularisation is based on assuming a multi-normal distribution for the
% elements of the Henckey Tensor. See:
% "Non-linear Elastic Deformations". R. W. Ogden (Dover), 1984.
% Weighting for the regularisation is determined approximately according
% to:
% "Incorporating Prior Knowledge into Image Registration"
% J. Ashburner, P. Neelin, D. L. Collins, A. C. Evans & K. J. Friston.
% NeuroImage 6:344-352 (1997).
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


if nargin<5, scal = ones(length(VG),1); end;
if nargin<4, M    = eye(4);             end;

def_flags = struct('sep',5, 'regtype','subj','WG',[],'WF',[],'globnorm',1,'debug',0);
if nargin < 2 || ~isstruct(flags),
    flags = def_flags;
else
    fnms = fieldnames(def_flags);
    for i=1:length(fnms),
        if ~isfield(flags,fnms{i}),
            flags.(fnms{i}) = def_flags.(fnms{i});
        end;
    end;
end;

% Check to ensure inputs are valid...
% ---------------------------------------------------------------
if length(VF)>1, error('Can not use more than one source image'); end;
if ~isempty(flags.WF),
    if length(flags.WF)>1,
        error('Can only use one source weighting image');
    end;
    if any(any((VF.mat-flags.WF.mat).^2>1e-8)),
        error('Source and its weighting image must have same orientation');
    end;
    if any(any(VF.dim(1:3)-flags.WF.dim(1:3))),
        error('Source and its weighting image must have same dimensions');
    end;
end;
if ~isempty(flags.WG),
    if length(flags.WG)>1,
        error('Can only use one template weighting image');
    end;
    tmp = reshape(cat(3,VG(:).mat,flags.WG.mat),16,length(VG)+length(flags.WG));
else
    tmp = reshape(cat(3,VG(:).mat),16,length(VG));
end;
if any(any(diff(tmp,1,2).^2>1e-8)),
    error('Reference images must all have the same orientation');
end;
if ~isempty(flags.WG),
    tmp = cat(1,VG(:).dim,flags.WG.dim);
else
    tmp = cat(1,VG(:).dim);
end;
if any(any(diff(tmp(:,1:3),1,1))),
    error('Reference images must all have the same dimensions');
end;
% ---------------------------------------------------------------

% Generate points to sample from, adding some jitter in order to
% make the cost function smoother.
% ---------------------------------------------------------------
rand('state',0); % want the results to be consistant.
dg   = VG(1).dim(1:3);
df   = VF(1).dim(1:3);

if length(VG)==1,
    skip = sqrt(sum(VG(1).mat(1:3,1:3).^2)).^(-1)*flags.sep;
    [x1,x2,x3]=ndgrid(1:skip(1):dg(1)-.5, 1:skip(2):dg(2)-.5, 1:skip(3):dg(3)-.5);
    x1   = x1 + rand(size(x1))*0.5; x1 = x1(:);
    x2   = x2 + rand(size(x2))*0.5; x2 = x2(:);
    x3   = x3 + rand(size(x3))*0.5; x3 = x3(:);
end;

skip = sqrt(sum(VF(1).mat(1:3,1:3).^2)).^(-1)*flags.sep;
[y1,y2,y3]=ndgrid(1:skip(1):df(1)-.5, 1:skip(2):df(2)-.5, 1:skip(3):df(3)-.5);
y1   = y1 + rand(size(y1))*0.5; y1 = y1(:);
y2   = y2 + rand(size(y2))*0.5; y2 = y2(:);
y3   = y3 + rand(size(y3))*0.5; y3 = y3(:);
% ---------------------------------------------------------------

if flags.globnorm,
    % Scale all images approximately equally
    % ---------------------------------------------------------------
    for i=1:length(VG),
        VG(i).pinfo(1:2,:) = VG(i).pinfo(1:2,:)/spm_global(VG(i));
    end;
    VF(1).pinfo(1:2,:) = VF(1).pinfo(1:2,:)/spm_global(VF(1));
end;
% ---------------------------------------------------------------

if length(VG)==1,
    [G,dG1,dG2,dG3]  = spm_sample_vol(VG(1),x1,x2,x3,1);
    if ~isempty(flags.WG),
        WG = abs(spm_sample_vol(flags.WG,x1,x2,x3,1))+eps;
        WG(~isfinite(WG)) = 1;
    end;
end;

[F,dF1,dF2,dF3]  = spm_sample_vol(VF(1),y1,y2,y3,1);
if ~isempty(flags.WF),
    WF = abs(spm_sample_vol(flags.WF,y1,y2,y3,1))+eps;
    WF(~isfinite(WF)) = 1;
end;
% ---------------------------------------------------------------
n_main_its = 0;
ss         = Inf;
W          = [Inf Inf Inf];
est_smo    = 1;
% ---------------------------------------------------------------

for iter=1:256,
    pss   = ss;
    p0    = [0 0 0  0 0 0  1 1 1  0 0 0];

    % Initialise the cost function and its 1st and second derivatives
    % ---------------------------------------------------------------
    n     = 0;
    ss    = 0;
    Beta  = zeros(12+length(VG),1);
    Alpha = zeros(12+length(VG));

    if length(VG)==1,
        % Make the cost function symmetric
        % ---------------------------------------------------------------

        % Build a matrix to rotate the derivatives by, converting from
        % derivatives w.r.t. changes in the overall affine transformation
        % matrix, to derivatives w.r.t. the parameters p.
        % ---------------------------------------------------------------
        dt  = 0.0001;
        R   = eye(13);
        MM0 = inv(VG.mat)*inv(spm_matrix(p0))*VG.mat;
        for i1=1:12,
            p1          = p0;
            p1(i1)      = p1(i1)+dt;
            MM1         = (inv(VG.mat)*inv(spm_matrix(p1))*(VG.mat));
            R(1:12,i1)  = reshape((MM1(1:3,:)-MM0(1:3,:))/dt,12,1);
        end;
        % ---------------------------------------------------------------
        [t1,t2,t3] = coords((M*VF(1).mat)\VG(1).mat,x1,x2,x3);
        msk        = find((t1>=1 & t1<=df(1) & t2>=1 & t2<=df(2) & t3>=1 & t3<=df(3)));
        if length(msk)<32, error_message; end;
        t1         = t1(msk);
        t2         = t2(msk);
        t3         = t3(msk);
        t          = spm_sample_vol(VF(1), t1,t2,t3,1);

        % Get weights
        % ---------------------------------------------------------------
        if ~isempty(flags.WF) || ~isempty(flags.WG),
            if isempty(flags.WF),
                wt = WG(msk);
            else
                wt = spm_sample_vol(flags.WF(1), t1,t2,t3,1)+eps;
                wt(~isfinite(wt)) = 1;
                if ~isempty(flags.WG), wt = 1./(1./wt + 1./WG(msk)); end;
            end;
            wt = sparse(1:length(wt),1:length(wt),wt);
        else
            % wt = speye(length(msk));
            wt = [];
        end;
        % ---------------------------------------------------------------
        clear t1 t2 t3

        % Update the cost function and its 1st and second derivatives.
        % ---------------------------------------------------------------
        [AA,Ab,ss1,n1] = costfun(x1,x2,x3,dG1,dG2,dG3,msk,scal^(-2)*t,G(msk)-(1/scal)*t,wt);
        Alpha = Alpha + R'*AA*R;
        Beta  = Beta  + R'*Ab;
        ss    = ss    + ss1;
        n     = n     + n1;
        % t     = G(msk) - (1/scal)*t;
    end;

    if 1,
        % Build a matrix to rotate the derivatives by, converting from
        % derivatives w.r.t. changes in the overall affine transformation
        % matrix, to derivatives w.r.t. the parameters p.
        % ---------------------------------------------------------------
        dt = 0.0001;
        R  = eye(12+length(VG));
        MM0 = inv(M*VF.mat)*spm_matrix(p0)*M*VF.mat;
        for i1=1:12,
            p1          = p0;
            p1(i1)      = p1(i1)+dt;
            MM1         = (inv(M*VF.mat)*spm_matrix(p1)*M*VF.mat);
            R(1:12,i1)  = reshape((MM1(1:3,:)-MM0(1:3,:))/dt,12,1);
        end;
        % ---------------------------------------------------------------
        [t1,t2,t3] = coords(VG(1).mat\M*VF(1).mat,y1,y2,y3);
        msk        = find((t1>=1 & t1<=dg(1) & t2>=1 & t2<=dg(2) & t3>=1 & t3<=dg(3)));
        if length(msk)<32, error_message; end;

        if length(msk)<32, error_message; end;
        t1 = t1(msk);
        t2 = t2(msk);
        t3 = t3(msk);
        t  = zeros(length(t1),length(VG));

        % Get weights
        % ---------------------------------------------------------------
        if ~isempty(flags.WF) || ~isempty(flags.WG),
            if isempty(flags.WG),
                wt = WF(msk);
            else
                wt = spm_sample_vol(flags.WG(1), t1,t2,t3,1)+eps;
                wt(~isfinite(wt)) = 1;
                if ~isempty(flags.WF), wt = 1./(1./wt + 1./WF(msk)); end;
            end;
            wt = sparse(1:length(wt),1:length(wt),wt);
        else
            wt = speye(length(msk));
        end;
        % ---------------------------------------------------------------

        if est_smo,
            % Compute derivatives of residuals in the space of F
            % ---------------------------------------------------------------
            [ds1,ds2,ds3] = transform_derivs(VG(1).mat\M*VF(1).mat,dF1(msk),dF2(msk),dF3(msk));
            for i=1:length(VG),
                [t(:,i),dt1,dt2,dt3] = spm_sample_vol(VG(i), t1,t2,t3,1);
                ds1   = ds1 - dt1*scal(i); clear dt1
                ds2   = ds2 - dt2*scal(i); clear dt2
                ds3   = ds3 - dt3*scal(i); clear dt3
            end;
            dss   = [ds1'*wt*ds1 ds2'*wt*ds2 ds3'*wt*ds3];
            clear ds1 ds2 ds3
        else
            for i=1:length(VG),
                t(:,i)= spm_sample_vol(VG(i), t1,t2,t3,1);
            end;
        end;

        clear t1 t2 t3

        % Update the cost function and its 1st and second derivatives.
        % ---------------------------------------------------------------
        [AA,Ab,ss2,n2] = costfun(y1,y2,y3,dF1,dF2,dF3,msk,-t,F(msk)-t*scal,wt);
        Alpha = Alpha  + R'*AA*R;
        Beta  = Beta   + R'*Ab;
        ss    = ss     + ss2;
        n     = n      + n2;
    end;

    if est_smo,
        % Compute a smoothness correction from the residuals and their
        % derivatives.  This is analagous to the one used in:
        %   "Analysis of fMRI Time Series Revisited"
        %   Friston KJ, Holmes AP, Poline JB, Grasby PJ, Williams SCR,
        %   Frackowiak RSJ, Turner R.  Neuroimage 2:45-53 (1995).
        % ---------------------------------------------------------------
        vx     = sqrt(sum(VG(1).mat(1:3,1:3).^2));
        pW     = W;
        W      = (2*dss/ss2).^(-.5).*vx;
        W      = min(pW,W);
        if flags.debug, fprintf('\nSmoothness FWHM:  %.3g x %.3g x %.3g mm\n', W*sqrt(8*log(2))); end;
        if length(VG)==1, dens=2; else dens=1; end;
        smo    = prod(min(dens*flags.sep/sqrt(2*pi)./W,[1 1 1]));
        est_smo=0;
        n_main_its = n_main_its + 1;
    end;

    % Update the parameter estimates
    % ---------------------------------------------------------------
    nu      = n*smo;
    sig2    = ss/nu;
    [d1,d2] = reg(M,12+length(VG),flags.regtype);

    soln    = (Alpha/sig2+d2)\(Beta/sig2-d1);
    scal    = scal - soln(13:end);
    M       = spm_matrix(p0 + soln(1:12)')*M;

    if flags.debug,
        fprintf('%d\t%g\n', iter, ss/n);
        piccies(VF,VG,M,scal)
    end;

    % If cost function stops decreasing, then re-estimate smoothness
    % and try again.  Repeat a few times.
    % ---------------------------------------------------------------
    ss = ss/n;
    if iter>1, spm_chi2_plot('Set',ss); end;
    if (pss-ss)/pss < 1e-6,
        est_smo = 1;
    end;
    if n_main_its>3, break; end;

end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [X1,Y1,Z1] = transform_derivs(Mat,X,Y,Z)
% Given the derivatives of a scalar function, return those of the
% affine transformed function
%_______________________________________________________________________

t1 = Mat(1:3,1:3);
t2 = eye(3);
if sum((t1(:)-t2(:)).^2) < 1e-12,
        X1 = X;Y1 = Y; Z1 = Z;
else
        X1    = Mat(1,1)*X + Mat(1,2)*Y + Mat(1,3)*Z;
        Y1    = Mat(2,1)*X + Mat(2,2)*Y + Mat(2,3)*Z;
        Z1    = Mat(3,1)*X + Mat(3,2)*Y + Mat(3,3)*Z;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [d1,d2] = reg(M,n,typ)
% Analytically compute the first and second derivatives of a penalty
% function w.r.t. changes in parameters.

if nargin<3, typ = 'subj'; end;
if nargin<2, n   = 13;     end;

[mu,isig] = priors(typ);
ds  = 0.000001;
d1  = zeros(n,1);
d2  = zeros(n);
p0  = [0 0 0  0 0 0  1 1 1  0 0 0];
h0  = penalty(p0,M,mu,isig);
for i=7:12, % derivatives are zero w.r.t. rotations and translations
    p1    = p0;
    p1(i) = p1(i)+ds;
    h1    = penalty(p1,M,mu,isig);
    d1(i) = (h1-h0)/ds; % First derivative
    for j=7:12,
        p2    = p0;
        p2(j) = p2(j)+ds;
        h2    = penalty(p2,M,mu,isig);
        p3    = p1;
        p3(j) = p3(j)+ds;
        h3    = penalty(p3,M,mu,isig);
        d2(i,j) = ((h3-h2)/ds-(h1-h0)/ds)/ds; % Second derivative
    end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function h = penalty(p,M,mu,isig)
% Return a penalty based on the elements of an affine transformation,
% which is given by:
%   spm_matrix(p)*M
%
% The penalty is based on the 6 unique elements of the Hencky tensor
% elements being multinormally distributed.
%_______________________________________________________________________

% Unique elements of symmetric 3x3 matrix.
els = [1 2 3 5 6 9];

T = spm_matrix(p)*M;
T = T(1:3,1:3);
T = 0.5*logm(T'*T);
T = T(els)' - mu;
h = T'*isig*T;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [mu,isig] = priors(typ)
% The parameters for this distribution were derived empirically from 227
% scans, that were matched to the ICBM space.
%_______________________________________________________________________

mu   = zeros(6,1);
isig = zeros(6);
switch deblank(lower(typ)),

case 'mni', % For registering with MNI templates...
    mu   = [0.0667 0.0039 0.0008 0.0333 0.0071 0.1071]';
    isig = 1e4 * [
        0.0902   -0.0345   -0.0106   -0.0025   -0.0005   -0.0163
       -0.0345    0.7901    0.3883    0.0041   -0.0103   -0.0116
       -0.0106    0.3883    2.2599    0.0113    0.0396   -0.0060
       -0.0025    0.0041    0.0113    0.0925    0.0471   -0.0440
       -0.0005   -0.0103    0.0396    0.0471    0.2964   -0.0062
       -0.0163   -0.0116   -0.0060   -0.0440   -0.0062    0.1144];

case 'rigid', % Constrained to be almost rigid...
    mu   = zeros(6,1);
    isig = eye(6)*1e9;

case 'isochoric', % Volume preserving...
    error('Not implemented');

case 'isotropic', % Isotropic zoom in all directions...
    error('Not implemented');

case 'subj', % For inter-subject registration...
    mu   = zeros(6,1);
    isig = 1e3 * [
        0.8876    0.0784    0.0784   -0.1749    0.0784   -0.1749
        0.0784    5.3894    0.2655    0.0784    0.2655    0.0784
        0.0784    0.2655    5.3894    0.0784    0.2655    0.0784
       -0.1749    0.0784    0.0784    0.8876    0.0784   -0.1749
        0.0784    0.2655    0.2655    0.0784    5.3894    0.0784
       -0.1749    0.0784    0.0784   -0.1749    0.0784    0.8876];

case 'none', % No regularisation...
    mu   = zeros(6,1);
    isig = zeros(6);

otherwise,
    error(['"' typ '" not recognised as type of regularisation.']);
end;
return;

%_______________________________________________________________________
function [y1,y2,y3]=coords(M,x1,x2,x3)
% Affine transformation of a set of coordinates.
%_______________________________________________________________________

y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function A = make_A(x1,x2,x3,dG1,dG2,dG3,t)
% Generate part of a design matrix using the chain rule...
% df/dm = df/dy * dy/dm
% where
%   df/dm is the rate of change of intensity w.r.t. affine parameters
%   df/dy is the gradient of the image f
%   dy/dm crange of position w.r.t. change of parameters
%_______________________________________________________________________

A  = [x1.*dG1 x1.*dG2 x1.*dG3 ...
      x2.*dG1 x2.*dG2 x2.*dG3 ...
      x3.*dG1 x3.*dG2 x3.*dG3 ...
          dG1     dG2     dG3    t];
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [AA,Ab,ss,n] = costfun(x1,x2,x3,dG1,dG2,dG3,msk,lastcols,b,wt)
chunk = 10240;
lm    = length(msk);
AA    = zeros(12+size(lastcols,2));
Ab    = zeros(12+size(lastcols,2),1);
ss    = 0;
n     = 0;

for i=1:ceil(lm/chunk),
    ind  = (((i-1)*chunk+1):min(i*chunk,lm))';
    msk1 = msk(ind);

    A1   = make_A(x1(msk1),x2(msk1),x3(msk1),dG1(msk1),dG2(msk1),dG3(msk1),lastcols(ind,:));
    b1   = b(ind);
    if ~isempty(wt),
        wt1   = wt(ind,ind);
        AA    = AA  + A1'*wt1*A1;
        %Ab   = Ab  + A1'*wt1*b1;
        Ab    = Ab  + (b1'*wt1*A1)';
        ss    = ss  + b1'*wt1*b1;
        n     = n   + trace(wt1);
        clear wt1
    else
        AA    = AA  + A1'*A1;
        %Ab   = Ab  + A1'*b1;
        Ab    = Ab  + (b1'*A1)';
        ss    = ss  + b1'*b1;
        n     = n   + length(msk1);
    end;
    clear A1 b1 msk1 ind
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function error_message
% Display an error message for when things go wrong.
str = { 'There is not enough overlap in the images',...
    'to obtain a solution.',...
    ' ',...
    'Please check that your header information is OK.',...
    'The Check Reg utility will show you the initial',...
    'alignment between the images, which must be',...
    'within about 4cm and about 15 degrees in order',...
    'for SPM to find the optimal solution.'};
spm('alert*',str,mfilename,sqrt(-1));
error('insufficient image overlap')
%_______________________________________________________________________

%_______________________________________________________________________
function piccies(VF,VG,M,scal)
% This is for debugging purposes.
% It shows the linear combination of template images, the affine
% transformed source image, the residual image and a histogram of the
% residuals.
%_______________________________________________________________________

figure(2);
Mt = spm_matrix([0 0 (VG(1).dim(3)+1)/2]);
M  = (M*VF(1).mat)\VG(1).mat;
t  = zeros(VG(1).dim(1:2));
for i=1:length(VG);
    t  = t + spm_slice_vol(VG(i),  Mt,VG(1).dim(1:2),1)*scal(i);
end;
u  = spm_slice_vol(VF(1),M*Mt,VG(1).dim(1:2),1);
subplot(2,2,1);imagesc(t');axis image xy off
subplot(2,2,2);imagesc(u');axis image xy off
subplot(2,2,3);imagesc(u'-t');axis image xy off
%subplot(2,2,4);hist(b,50); % Entropy of residuals may be a nice cost function?
drawnow;
return;
%_______________________________________________________________________
