function VO = spm_write_sn(V,prm,flags,extras)
% Write out warped images
% FORMAT VO = spm_write_sn(V,prm,flags,msk)
% V         - Images to transform (filenames or volume structure).
% prm       - Transformation information (filename or structure).
% flags     - flags structure, with fields...
%           interp   - interpolation method (0-7)
%           wrap     - wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%           vox      - voxel sizes (3 element vector - in mm)
%                      Non-finite values mean use template vox.
%           bb       - bounding box (2x3 matrix - in mm)
%                      Non-finite values mean use template bb.
%           preserve - either 0 or 1.  A value of 1 will "modulate"
%                      the spatially normalised images so that total
%                      units are preserved, rather than just
%                      concentrations.
%           prefix   - Prefix for normalised images. Defaults to 'w'.
% msk       - An optional cell array for masking the spatially
%             normalised images (see below).
%
% Warped images are written prefixed by "w".
%
% Non-finite vox or bounding box suggests that values should be derived
% from the template image.
%
% Don't use interpolation methods greater than one for data containing
% NaNs.
%__________________________________________________________________________
%
% FORMAT msk = spm_write_sn(V,prm,flags,'mask')
% V          - Images to transform (filenames or volume structure).
% prm        - Transformation information (filename or structure).
% flags      - flags structure, with fields...
%            wrap - wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%            vox  - voxel sizes (3 element vector - in mm)
%                   Non-finite values mean use template vox.
%            bb   - bounding box (2x3 matrix - in mm)
%                   Non-finite values mean use template bb.
% msk        - a cell array for masking a series of spatially normalised
%              images.
%
%
%_________________________________________________________________________
%
% FORMAT VO = spm_write_sn(V,prm,'modulate')
% V         - Spatially normalised images to modulate (filenames or
%             volume structure).
% prm       - Transformation information (filename or structure).
%
%  After nonlinear spatial normalization, the relative volumes of some
%  brain structures will have decreased, whereas others will increase.
%  The resampling of the images preserves the concentration of pixel
%  units in the images, so the total counts from structures that have
%  reduced volumes after spatial normalization will be reduced by an
%  amount proportional to the volume reduction.
%
%  This routine rescales images after spatial normalization, so that
%  the total counts from any structure are preserved.  It was written
%  as an optional step in performing voxel based morphometry.
%
%__________________________________________________________________________
% Copyright (C) 1996-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_write_sn.m 4878 2012-08-31 12:07:30Z guillaume $


SVNid = '$Rev: 4878 $';
 
%-Say hello
%--------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SVNid);

%-Get images
%--------------------------------------------------------------------------
if isempty(V), return; end

if ischar(prm), prm = load(prm);  end
if ischar(V),   V   = spm_vol(V); end

%-Modulate
%--------------------------------------------------------------------------
if nargin==3 && ischar(flags) && strcmpi(flags,'modulate')
    if nargout==0
        modulate(V,prm);
    else
        VO = modulate(V,prm);
    end
    return;
end

%-Get parameters
%--------------------------------------------------------------------------
def_flags        = spm_get_defaults('old.normalise.write');

if nargin < 3
    flags = def_flags;
else
    fnms = fieldnames(def_flags);
    for i=1:length(fnms)
        if ~isfield(flags,fnms{i})
            flags.(fnms{i}) = def_flags.(fnms{i});
        end
    end
end

%-Mask
%--------------------------------------------------------------------------
[x,y,z,mat] = get_xyzmat(prm,flags.bb,flags.vox);

if nargin==4
    if ischar(extras) && strcmpi(extras,'mask')
        VO = get_snmask(V,prm,x,y,z,flags.wrap);
        return;
    end
    if iscell(extras)
        msk = extras;
    end
end

if nargout>0 && numel(V)>8
    error('Too many images to save in memory');
end

%-Warp images
%--------------------------------------------------------------------------
if ~exist('msk','var')
    msk = get_snmask(V,prm,x,y,z,flags.wrap);
end

if nargout==0
    if isempty(prm.Tr)
        affine_transform(V,prm,x,y,z,mat,flags,msk);
    else
        nonlin_transform(V,prm,x,y,z,mat,flags,msk);
    end
else
    if isempty(prm.Tr)
        VO = affine_transform(V,prm,x,y,z,mat,flags,msk);
    else
        VO = nonlin_transform(V,prm,x,y,z,mat,flags,msk);
    end
end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#


%==========================================================================
%-function VO = affine_transform(V,prm,x,y,z,mat,flags,msk)
%==========================================================================
function VO = affine_transform(V,prm,x,y,z,mat,flags,msk)

[X,Y] = ndgrid(x,y);
d     = [flags.interp*[1 1 1]' flags.wrap(:)];

spm_progress_bar('Init',numel(V),'Resampling','volumes/slices completed');
for i=1:numel(V)
    VO     = make_hdr_struct(V(i),x,y,z,mat, flags.prefix);
    if flags.preserve
        VO.fname = spm_file(VO.fname, 'prefix','m');
    end
    detAff = det(prm.VF(1).mat*prm.Affine/prm.VG(1).mat);
    if flags.preserve, VO.pinfo(1:2,:) = VO.pinfo(1:2,:)/detAff; end

    %Dat= zeros(VO.dim(1:3));
    Dat = single(0);
    Dat(VO.dim(1),VO.dim(2),VO.dim(3)) = 0;

    C = spm_bsplinc(V(i),d);

    for j=1:length(z)   % Cycle over planes
        [X2,Y2,Z2]  = mmult(X,Y,z(j),V(i).mat\prm.VF(1).mat*prm.Affine);
        dat         = spm_bsplins(C,X2,Y2,Z2,d);
        if flags.preserve, dat = dat*detAff; end
        dat(msk{j}) = NaN;

        Dat(:,:,j) = single(dat);

        if numel(V)<5, spm_progress_bar('Set',i-1+j/length(z)); end
    end
    if nargout~=0
        VO.pinfo  = [1 0]';
        VO.dt     = [spm_type('float32') spm_platform('bigend')];
        VO.dat    = Dat;
    else
        spm_write_vol(VO, Dat);
    end
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


%==========================================================================
%-function VO = nonlin_transform(V,prm,x,y,z,mat,flags,msk)
%==========================================================================
function VO = nonlin_transform(V,prm,x,y,z,mat,flags,msk)

[X,Y] = ndgrid(x,y);
Tr = prm.Tr;
BX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1);
BY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1);
BZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1);
if flags.preserve
    DX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1,'diff');
    DY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1,'diff');
    DZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1,'diff');
end
d  = [flags.interp*[1 1 1]' flags.wrap(:)];

spm_progress_bar('Init',numel(V),'Resampling','volumes completed');
for i=1:numel(V)
    VO     = make_hdr_struct(V(i),x,y,z,mat, flags.prefix);
    if flags.preserve
        VO.fname = spm_file(VO.fname, 'prefix','m');
    end
    detAff = det(prm.VF(1).mat*prm.Affine/prm.VG(1).mat);

    % Accumulate data
    %Dat= zeros(VO.dim(1:3));
    Dat = single(0);
    Dat(VO.dim(1),VO.dim(2),VO.dim(3)) = 0;

    C = spm_bsplinc(V(i),d);

    for j=1:length(z)   % Cycle over planes
        % Nonlinear deformations
        %------------------------------------------------------------------
        tx = get_2Dtrans(Tr(:,:,:,1),BZ,j);
        ty = get_2Dtrans(Tr(:,:,:,2),BZ,j);
        tz = get_2Dtrans(Tr(:,:,:,3),BZ,j);
        X1 = X    + BX*tx*BY';
        Y1 = Y    + BX*ty*BY';
        Z1 = z(j) + BX*tz*BY';

        [X2,Y2,Z2]  = mmult(X1,Y1,Z1,V(i).mat\prm.VF(1).mat*prm.Affine);
        dat         = spm_bsplins(C,X2,Y2,Z2,d);
        dat(msk{j}) = NaN;

        if ~flags.preserve
            Dat(:,:,j) = single(dat);
        else
            j11 = DX*tx*BY' + 1; j12 = BX*tx*DY';     j13 = BX*get_2Dtrans(Tr(:,:,:,1),DZ,j)*BY';
            j21 = DX*ty*BY';     j22 = BX*ty*DY' + 1; j23 = BX*get_2Dtrans(Tr(:,:,:,2),DZ,j)*BY';
            j31 = DX*tz*BY';     j32 = BX*tz*DY';     j33 = BX*get_2Dtrans(Tr(:,:,:,3),DZ,j)*BY' + 1;

            % The determinant of the Jacobian reflects relative volume changes
            %-----------------------------------------------------------------
            dat       = dat .* (j11.*(j22.*j33-j23.*j32) - j21.*(j12.*j33-j13.*j32) + j31.*(j12.*j23-j13.*j22)) * detAff;
            Dat(:,:,j) = single(dat);
        end
        if numel(V)<5, spm_progress_bar('Set',i-1+j/length(z)); end
    end
    if nargout==0
        if flags.preserve, VO = rmfield(VO,'pinfo'); end
        VO = spm_write_vol(VO,Dat);
    else
        VO.pinfo  = [1 0]';
        VO.dt     = [spm_type('float32') spm_platform('bigend')];
        VO.dat    = Dat;
    end
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


%==========================================================================
%-function VO = modulate(V,prm)
%==========================================================================
function VO = modulate(V,prm)

spm_progress_bar('Init',numel(V),'Modulating','volumes completed');
for i=1:numel(V)
    VO          = V(i);
    VO          = rmfield(VO,'pinfo');
    VO.fname    = spm_file(VO.fname, 'prefix','m');
    detAff      = det(prm.VF(1).mat*prm.Affine/prm.VG(1).mat);
    %Dat        = zeros(VO.dim(1:3));
    Dat         = single(0);
    Dat(VO.dim(1),VO.dim(2),VO.dim(3)) = 0;
    [x,y,z,mat] = get_xyzmat(prm,NaN,NaN,VO);

    if sum((mat(:)-VO.mat(:)).^2)>1e-7, error('Orientations not compatible'); end

    Tr = prm.Tr;

    if isempty(Tr)
        for j=1:length(z),   % Cycle over planes
            dat        = spm_slice_vol(V(i),spm_matrix([0 0 j]),V(i).dim(1:2),0);
            Dat(:,:,j) = single(dat);
            if numel(V)<5, spm_progress_bar('Set',i-1+j/length(z)); end
        end
    else
        BX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1);
        BY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1);
        BZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1);
        DX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1,'diff');
        DY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1,'diff');
        DZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1,'diff');

        for j=1:length(z)   % Cycle over planes

            tx = get_2Dtrans(Tr(:,:,:,1),BZ,j);
            ty = get_2Dtrans(Tr(:,:,:,2),BZ,j);
            tz = get_2Dtrans(Tr(:,:,:,3),BZ,j);

            j11 = DX*tx*BY' + 1; j12 = BX*tx*DY';     j13 = BX*get_2Dtrans(Tr(:,:,:,1),DZ,j)*BY';
            j21 = DX*ty*BY';     j22 = BX*ty*DY' + 1; j23 = BX*get_2Dtrans(Tr(:,:,:,2),DZ,j)*BY';
            j31 = DX*tz*BY';     j32 = BX*tz*DY';     j33 = BX*get_2Dtrans(Tr(:,:,:,3),DZ,j)*BY' + 1;

            % The determinant of the Jacobian reflects relative volume changes
            %-----------------------------------------------------------------
            dat        = spm_slice_vol(V(i),spm_matrix([0 0 j]),V(i).dim(1:2),0);
            dat        = dat .* (j11.*(j22.*j33-j23.*j32) - j21.*(j12.*j33-j13.*j32) + j31.*(j12.*j23-j13.*j22)) * detAff;
            Dat(:,:,j) = single(dat);

            if numel(V)<5, spm_progress_bar('Set',i-1+j/length(z)); end
        end
    end

    if nargout==0
        VO = spm_write_vol(VO,Dat);
    else
        VO.pinfo  = [1 0]';
        VO.dt     = [spm_type('float32') spm_platform('bigend')];
        VO.dat    = Dat;
    end
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


%==========================================================================
%-function VO   = make_hdr_struct(V,x,y,z,mat,prefix)
%==========================================================================
function VO   = make_hdr_struct(V,x,y,z,mat,prefix)
VO            = V;
VO.fname      = spm_file(V.fname, 'prefix',prefix);
VO.mat        = mat;
VO.dim(1:3)   = [length(x) length(y) length(z)];
VO.pinfo      = V.pinfo;
VO.descrip    = 'spm - 3D normalized';


%==========================================================================
%-function T2 = get_2Dtrans(T3,B,j)
%==========================================================================
function T2 = get_2Dtrans(T3,B,j)
d   = [size(T3) 1 1 1];
tmp = reshape(T3,d(1)*d(2),d(3));
T2  = reshape(tmp*B(j,:)',d(1),d(2));


%==========================================================================
%-function Mask = getmask(X,Y,Z,dim,wrp)
%==========================================================================
function Mask = getmask(X,Y,Z,dim,wrp)
% Find range of slice
tiny = 5e-2;
Mask = true(size(X));
if ~wrp(1), Mask = Mask & (X >= (1-tiny) & X <= (dim(1)+tiny)); end
if ~wrp(2), Mask = Mask & (Y >= (1-tiny) & Y <= (dim(2)+tiny)); end
if ~wrp(3), Mask = Mask & (Z >= (1-tiny) & Z <= (dim(3)+tiny)); end


%==========================================================================
%-function [X2,Y2,Z2] = mmult(X1,Y1,Z1,Mult)
%==========================================================================
function [X2,Y2,Z2] = mmult(X1,Y1,Z1,Mult)
if length(Z1) == 1
    X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + (Mult(1,3)*Z1 + Mult(1,4));
    Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + (Mult(2,3)*Z1 + Mult(2,4));
    Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + (Mult(3,3)*Z1 + Mult(3,4));
else
    X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
    Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
    Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
end


%==========================================================================
%-function msk = get_snmask(V,prm,x,y,z,wrap)
%==========================================================================
function msk = get_snmask(V,prm,x,y,z,wrap)
% Generate a mask for where there is data for all images
%--------------------------------------------------------------------------
msk = cell(length(z),1);
t1 = cat(3,V.mat);
t2 = cat(1,V.dim);
t  = [reshape(t1,[16 length(V)])' t2(:,1:3)];
Tr = prm.Tr;
[X,Y] = ndgrid(x,y);

BX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1);
BY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1);
BZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1);

if numel(V)>1 && any(any(diff(t,1,1)))
    spm_progress_bar('Init',length(z),'Computing available voxels','planes completed');
    for j=1:length(z)   % Cycle over planes
        Count = zeros(length(x),length(y));
        if isempty(Tr)
            % Generate a mask for where there is data for all images
            %--------------------------------------------------------------
            for i=1:numel(V),
                [X2,Y2,Z2] = mmult(X,Y,z(j),V(i).mat\prm.VF(1).mat*prm.Affine);
                Count      = Count + getmask(X2,Y2,Z2,V(i).dim(1:3),wrap);
            end
        else
            % Nonlinear deformations
            %--------------------------------------------------------------
            X1 = X    + BX*get_2Dtrans(Tr(:,:,:,1),BZ,j)*BY';
            Y1 = Y    + BX*get_2Dtrans(Tr(:,:,:,2),BZ,j)*BY';
            Z1 = z(j) + BX*get_2Dtrans(Tr(:,:,:,3),BZ,j)*BY';

            % Generate a mask for where there is data for all images
            %--------------------------------------------------------------
            for i=1:numel(V)
                [X2,Y2,Z2] = mmult(X1,Y1,Z1,V(i).mat\prm.VF(1).mat*prm.Affine);
                Count      = Count + getmask(X2,Y2,Z2,V(i).dim(1:3),wrap);
            end
        end
        msk{j} = uint32(find(Count ~= numel(V)));
        spm_progress_bar('Set',j);
    end
    spm_progress_bar('Clear');
else
    for j=1:length(z), msk{j} = uint32([]); end
end


%==========================================================================
%-function [x,y,z,mat] = get_xyzmat(prm,bb,vox,VG)
%==========================================================================
function [x,y,z,mat] = get_xyzmat(prm,bb,vox,VG)
% The old voxel size and origin notation is used here.
% This requires that the position and orientation
% of the template is transverse.  It would not be
% straitforward to account for templates that are
% in different orientations because the basis functions
% would no longer be seperable.  The seperable basis
% functions mean that computing the deformation field
% from the parameters is much faster.

% bb  = sort(bb);
% vox = abs(vox);

if nargin<4
    VG = prm.VG(1);

    if all(~isfinite(bb(:))) && all(~isfinite(vox(:)))
        x   = 1:VG.dim(1);
        y   = 1:VG.dim(2);
        z   = 1:VG.dim(3);
        mat = VG.mat;
        return;
    end
end

[bb0 vox0] = spm_get_bbox(VG, 'old');
if ~all(isfinite(vox(:))), vox = vox0; end
if ~all(isfinite(bb(:))),  bb  = bb0;  end

msk       = find(vox<0);
bb        = sort(bb);
bb(:,msk) = flipud(bb(:,msk));

% Adjust bounding box slightly - so it rounds to closest voxel.
% Comment out if not needed.
%bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
%bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
%bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

M   = prm.VG(1).mat;
vxg = sqrt(sum(M(1:3,1:3).^2));
if det(M(1:3,1:3))<0, vxg(1) = -vxg(1); end
ogn = M\[0 0 0 1]';
ogn = ogn(1:3)';

% Convert range into range of voxels within template image
x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);

og  = -vxg.*ogn;

% Again, chose whether to round to closest voxel.
%of  = -vox.*(round(-bb(1,:)./vox)+1);
of = bb(1,:)-vox;

M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
mat = prm.VG(1).mat*inv(M1)*M2;

LEFTHANDED = true;
if (LEFTHANDED && det(mat(1:3,1:3))>0) || (~LEFTHANDED && det(mat(1:3,1:3))<0),
    Flp = [-1 0 0 (length(x)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
    mat = mat*Flp;
    x   = flipud(x(:))';
end


%==========================================================================
%-function VO = write_dets(P,bb,vox)
%==========================================================================
function VO = write_dets(P,bb,vox)
if nargin==1
    job = P;
    P   = job.P;
    bb  = job.bb;
    vox = job.vox;
end

spm_progress_bar('Init',numel(P),'Writing','volumes completed');

for i=1:numel(V)
    prm = load(deblank(P{i}));
    [x,y,z,mat] = get_xyzmat(prm,bb,vox);
    Tr = prm.Tr;
    BX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1);
    BY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1);
    BZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1);
    DX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1,'diff');
    DY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1,'diff');
    DZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1,'diff');

    VO = struct(...
        'fname',   spm_file(P{i}, 'prefix','jy_'),...
        'dim',     [numel(x),numel(y),numel(z)],...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'pinfo',   [1 0 0]',...
        'mat',     mat,...
        'n',       1,...
        'descrip', 'Jacobian determinants');
    VO     = spm_create_vol(VO);
    detAff = det(prm.VF(1).mat*prm.Affine/prm.VG(1).mat);
    Dat    = single(0);
    Dat(VO.dim(1),VO.dim(2),VO.dim(3)) = 0;

    for j=1:length(z)   % Cycle over planes
        % Nonlinear deformations
        tx = get_2Dtrans(Tr(:,:,:,1),BZ,j);
        ty = get_2Dtrans(Tr(:,:,:,2),BZ,j);
        tz = get_2Dtrans(Tr(:,:,:,3),BZ,j);

        %------------------------------------------------------------------
        j11 = DX*tx*BY' + 1; j12 = BX*tx*DY';     j13 = BX*get_2Dtrans(Tr(:,:,:,1),DZ,j)*BY';
        j21 = DX*ty*BY';     j22 = BX*ty*DY' + 1; j23 = BX*get_2Dtrans(Tr(:,:,:,2),DZ,j)*BY';
        j31 = DX*tz*BY';     j32 = BX*tz*DY';     j33 = BX*get_2Dtrans(Tr(:,:,:,3),DZ,j)*BY' + 1;
        % The determinant of the Jacobian reflects relative volume changes.
        %------------------------------------------------------------------
        dat       = (j11.*(j22.*j33-j23.*j32) - j21.*(j12.*j33-j13.*j32) + j31.*(j12.*j23-j13.*j22)) * detAff;
        Dat(:,:,j) = single(dat);
        if numel(P)<5, spm_progress_bar('Set',i-1+j/length(z)); end
    end
    VO = spm_write_vol(VO,Dat);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');
