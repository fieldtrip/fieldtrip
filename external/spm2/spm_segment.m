function VO = spm_segment(VF,PG,flags)
% Segment an MR image into Gray, White & CSF.
%
% FORMAT VO = spm_segment(PF,PG,flags)
% PF    - name(s) of image(s) to segment (must have same dimensions).
% PG    - name(s) of template image(s) for realignment.
%       - or a 4x4 transformation matrix which maps from the image to
%         the set of templates.
% flags - a structure normally based on defaults.segment
% VO    - optional output volume
%
%                      The algorithm is four step:
%
% 1) Determine the affine transform which best matches the image with a
%    template image. If the name of more than one image is passed, then
%    the first image is used in this step. This step is not performed if
%    no template images are specified.
%
% 2) Perform Cluster Analysis with a modified Mixture Model and a-priori
%    information about the likelihoods of each voxel being one of a
%    number of different tissue types. If more than one image is passed,
%    then they they are all assumed to be in register, and the voxel
%    values are fitted to multi-normal distributions.
%
% 3) Perform morphometric operations on the grey and white partitions
%    in order to more accurately identify brain tissue. This is then used
%    to clean up the grey and white matter segments. 
%
% 4) If no output argument is specified, then the segmented images are
%    written to disk. The names of these images have "_seg1", "_seg2"
%    & "_seg3" appended to the name of the first image passed.
%
%_______________________________________________________________________
% Refs:
%
% Ashburner J & Friston KJ (1997) Multimodal Image Coregistration and
% Partitioning - a Unified Framework. NeuroImage 6:209-217
%
%_______________________________________________________________________
%
% The template image, and a-priori likelihood images are modified
% versions of those kindly supplied by Alan Evans, MNI, Canada
% (ICBM, NIH P-20 project, Principal Investigator John Mazziotta).
%_______________________________________________________________________
% @(#)spm_segment.m	2.28a John Ashburner 04/04/01

% Create some suitable default values
%-----------------------------------------------------------------------

def_flags.estimate.priors = str2mat(...
        fullfile(spm('Dir'),'apriori','gray.mnc'),...
        fullfile(spm('Dir'),'apriori','white.mnc'),...
        fullfile(spm('Dir'),'apriori','csf.mnc'));
def_flags.estimate.reg    = 0.01;
def_flags.estimate.cutoff = 30;
def_flags.estimate.samp   = 3;
def_flags.estimate.bb     =  [[-88 88]' [-122 86]' [-60 95]'];
def_flags.estimate.affreg.smosrc = 8;
def_flags.estimate.affreg.regtype = 'mni';
def_flags.estimate.affreg.weight = '';
def_flags.write.cleanup   = 1;
def_flags.write.wrt_cor   = 1;
def_flags.graphics        = 1;

if nargin<3, flags = def_flags; end;
if ~isfield(flags,'estimate'),        flags.estimate        = def_flags.estimate;        end;
if ~isfield(flags.estimate,'priors'), flags.estimate.priors = def_flags.estimate.priors; end;
if ~isfield(flags.estimate,'reg'),    flags.estimate.reg    = def_flags.estimate.reg;    end;
if ~isfield(flags.estimate,'cutoff'), flags.estimate.cutoff = def_flags.estimate.cutoff; end;
if ~isfield(flags.estimate,'samp'),   flags.estimate.samp   = def_flags.estimate.samp;   end;
if ~isfield(flags.estimate,'bb'),     flags.estimate.bb     = def_flags.estimate.bb;     end;
if ~isfield(flags.estimate,'affreg'), flags.estimate.affreg = def_flags.estimate.affreg; end;
if ~isfield(flags.estimate.affreg,'smosrc'),
	flags.estimate.affreg.smosrc = def_flags.estimate.affreg.smosrc;
end;
if ~isfield(flags.estimate.affreg,'regtype'),
	flags.estimate.affreg.regtype = def_flags.estimate.affreg.regtype;
end;
if ~isfield(flags.estimate.affreg,'weight'),
	flags.estimate.affreg.weight = def_flags.estimate.affreg.weight;
end;
if ~isfield(flags,'write'),         flags.write         = def_flags.write;         end;
if ~isfield(flags.write,'cleanup'), flags.write.cleanup = def_flags.write.cleanup; end;
if ~isfield(flags.write,'wrt_cor'), flags.write.wrt_cor = def_flags.write.wrt_cor; end;
if ~isfield(flags,'graphics'),      flags.graphics      = def_flags.graphics;      end;

%-----------------------------------------------------------------------

if ischar(VF), VF= spm_vol(VF); end;

SP         = init_sp(flags.estimate,VF,PG);
[x1,x2,x3] = get_sampling(SP.MM,VF,flags.estimate.samp,flags.estimate.bb);
BP         = init_bp(VF, flags.estimate.cutoff, flags.estimate.reg);
CP         = init_cp(VF,x3);
sums       = zeros(8,1);

for pp=1:length(x3),
	[raw,msk] = get_raw(VF,x1,x2,x3(pp));
	s         = get_sp(SP,x1,x2,x3(pp));
	CP        = update_cp_est(CP,s,raw,msk,pp);
	sums      = sums + reshape(sum(sum(s,1),2),8,1);
end;
sums = sums/sum(sums);
CP   = shake_cp(CP);

[CP,BP,SP] = run_segment(CP,BP,SP,VF,sums,x1,x2,x3);

%save segmentation_results.mat CP BP SP VF sums

[g,w,c] = get_gwc(VF,BP,SP,CP,sums,flags.write.wrt_cor);
if flags.write.cleanup, [g,w,c] = clean_gwc(g,w,c); end;

% Create the segmented images.
%-----------------------------------------------------------------------
%offs  = cumsum(repmat(prod(VF(1).dim(1:2)),1,VF(1).dim(3)))-prod(VF(1).dim(1:2));
%pinfo = [repmat([1/255 0]',1,VF(1).dim(3)) ; offs];
[pth,nm,xt] = fileparts(deblank(VF(1).fname));
for j=1:3,
	tmp   = fullfile(pth,[nm '_seg' num2str(j) xt]);
	VO(j) = struct(...
		'fname',tmp,...
		'dim',    [VF(1).dim(1:3) 2],...
		'mat',    VF(1).mat,...
		'pinfo',  [1/255 0 0]',...
		'descrip','Segmented image');
end;

if nargout==0,
	VO = spm_create_vol(VO);

	spm_progress_bar('Init',VF(1).dim(3),'Writing Segmented','planes completed');
	for pp=1:VF(1).dim(3),
		VO(1) = spm_write_plane(VO(1),double(g(:,:,pp))/255,pp);
		VO(2) = spm_write_plane(VO(2),double(w(:,:,pp))/255,pp);
		VO(3) = spm_write_plane(VO(3),double(c(:,:,pp))/255,pp);
		spm_progress_bar('Set',pp);
	end;
	VO = spm_close_vol(VO);
	spm_progress_bar('Clear');
end;

VO(1).dat = g; VO(1).pinfo = VO(1).pinfo(1:2,:);
VO(2).dat = w; VO(2).pinfo = VO(2).pinfo(1:2,:);
VO(3).dat = c; VO(3).pinfo = VO(3).pinfo(1:2,:);

if flags.graphics, display_graphics(VF,VO,CP.mn,CP.cv,CP.mg); end;

return;
%=======================================================================

%=======================================================================
function [y1,y2,y3] = affine_transform(x1,x2,x3,M)
	y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
	y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
	y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);
return;
%=======================================================================

%=======================================================================
function display_graphics(VF,VS,mn,cv,mg)
% Do the graphics
nb = 3;
spm_figure('Clear','Graphics');
fg = spm_figure('FindWin','Graphics');
if ~isempty(fg),
	% Show some text
	%-----------------------------------------------------------------------
	ax = axes('Position',[0.05 0.8 0.9 0.2],'Visible','off','Parent',fg);
	text(0.5,0.80, 'Segmentation','FontSize',16,'FontWeight','Bold',...
		'HorizontalAlignment','center','Parent',ax);

	text(0,0.65, ['Image:  ' spm_str_manip(VF(1).fname,'k50d')],...
		'FontSize',14,'FontWeight','Bold','Parent',ax);

	text(0,0.40, 'Means:','FontSize',12,'FontWeight','Bold','Parent',ax);
	text(0,0.30, 'Std devs:' ,'FontSize',12,'FontWeight','Bold','Parent',ax);
	text(0,0.20, 'N vox:','FontSize',12,'FontWeight','Bold','Parent',ax);
	for j=1:nb,
		text((j+0.5)/(nb+1),0.40, num2str(mn(1,j)),...
			'FontSize',12,'FontWeight','Bold',...
			'HorizontalAlignment','center','Parent',ax);
		text((j+0.5)/(nb+1),0.30, num2str(sqrt(cv(1,1,j))),...
			'FontSize',12,'FontWeight','Bold',...
			'HorizontalAlignment','center','Parent',ax);
		text((j+0.5)/(nb+1),0.20, num2str(mg(1,j)/sum(mg(1,:))),...
			'FontSize',12,'FontWeight','Bold',...
			'HorizontalAlignment','center','Parent',ax);
	end;
	if length(VF) > 1,
		text(0,0.10,...
		'Note: only means and variances for the first image are shown',...
		'Parent',ax,'FontSize',12);
	end;

	M1 = VS(1).mat;
	M2 = VF(1).mat;
	for i=1:5,
		M   = spm_matrix([0 0 i*VF(1).dim(3)/6]);
		img = spm_slice_vol(VF(1),M,VF(1).dim(1:2),1);
		img(1,1) = eps;
		ax = axes('Position',...
			[0.05 0.75*(1-i/5)+0.05 0.9/(nb+1) 0.75/5],...
			'Visible','off','Parent',fg);
		imagesc(rot90(img), 'Parent', ax);
		set(ax,'Visible','off','DataAspectRatio',[1 1 1]);

		for j=1:3,
			img = spm_slice_vol(VS(j),M2\M1*M,VF(1).dim(1:2),1);
			ax  = axes('Position',...
				[0.05+j*0.9/(nb+1) 0.75*(1-i/5)+0.05 0.9/(nb+1) 0.75/5],...
				'Visible','off','Parent',fg);
			image(rot90(img*64), 'Parent', ax);
			set(ax,'Visible','off','DataAspectRatio',[1 1 1]);
		end;
	end;

	spm_print;
	drawnow;
end;
return;
%=======================================================================

%=======================================================================
function M = get_affine_mapping(VF,VG,aflags)

if ~isempty(VG) & ischar(VG), VG = spm_vol(VG); end;

if ~isempty(VG) & isstruct(VG),
	% Affine registration so that a priori images match the image to
	% be segmented.
	%-----------------------------------------------------------------------

	VFS = spm_smoothto8bit(VF(1),aflags.smosrc);

	% Scale all images approximately equally
	% ---------------------------------------------------------------
	for i=1:length(VG),
		VG(i).pinfo(1:2,:) = VG(i).pinfo(1:2,:)/spm_global(VG(i));
	end;
	VFS(1).pinfo(1:2,:) = VFS(1).pinfo(1:2,:)/spm_global(VFS(1));

	spm_chi2_plot('Init','Affine Registration','Mean squared difference','Iteration');
	flags     = struct('sep',aflags.smosrc, 'regtype',aflags.regtype,'WG',[],'globnorm',0,'debug',0);
	M         = eye(4);
	[M,scal]  = spm_affreg(VG, VFS, flags, M);

	if ~isempty(aflags.weight), flags.WG = spm_vol(aflags.weight); end;

	flags.sep = aflags.smosrc/2;
	M         = spm_affreg(VG, VFS, flags, M,scal);
	spm_chi2_plot('Clear');

elseif all(size(VG) == [4 4])
	% Assume that second argument is a matrix that will do the job
	%-----------------------------------------------------------------------
	M = VG;
else
	% Assume that image is normalized
	%-----------------------------------------------------------------------
	M = eye(4);
end
return;
%=======================================================================

%=======================================================================
function [x1,x2,x3] = get_sampling(MM,VF,samp,bb1)
% Voxels to sample during the cluster analysis
%-----------------------------------------------------------------------

% A bounding box for the brain in Talairach space.
%bb = [ [-88 88]' [-122 86]' [-60 95]'];
%c = [bb(1,1) bb(1,2) bb(1,3) 1
%     bb(1,1) bb(1,2) bb(2,3) 1
%     bb(1,1) bb(2,2) bb(1,3) 1
%     bb(1,1) bb(2,2) bb(2,3) 1
%     bb(2,1) bb(1,2) bb(1,3) 1
%     bb(2,1) bb(1,2) bb(2,3) 1
%     bb(2,1) bb(2,2) bb(1,3) 1
%     bb(2,1) bb(2,2) bb(2,3) 1]';
%tc = MM\c;
%tc = tc(1:3,:)';
%mx = max(tc);
%mn = min(tc);
%bb = [mn ; mx];
%vx = sqrt(sum(VF(1).mat(1:3,1:3).^2));
%samp = round(max(abs([4 4 4]./vx), [1 1 1]));
%x1 = bb(1,1):samp(1):bb(2,1);
%x2 = bb(1,2):samp(2):bb(2,2);
%x3 = bb(1,3):samp(3):bb(2,3);
%return;

% A bounding box for the brain in Talairach space.
if nargin<4, bb1 = [ [-88 88]' [-122 86]' [-60 95]']; end;

% A mapping from a unit radius sphere to a hyper-ellipse
% that is just enclosed by the bounding box in Talairach
% space.
M0 = [diag(diff(bb1)/2) mean(bb1)';[0 0 0 1]];

% The mapping from voxels to Talairach space is MM,
% so the ellipse in the space of the image becomes:
M0 = MM\M0;

% So to work out the bounding box in the space of the
% image that just encloses the hyper-ellipse.
tmp = M0(1:3,1:3);
tmp = diag(tmp*tmp'/diag(sqrt(diag(tmp*tmp'))));
bb  = round([M0(1:3,4)-tmp M0(1:3,4)+tmp])';
bb  = min(max(bb,[1 1 1 ; 1 1 1]),[VF(1).dim(1:3) ; VF(1).dim(1:3)]);

% Want to sample about every 3mm
tmp  = sqrt(sum(VF(1).mat(1:3,1:3).^2))';
samp = round(max(abs(tmp.^(-1)*samp), [1 1 1]'));

x1 = bb(1,1):samp(1):bb(2,1);
x2 = bb(1,2):samp(2):bb(2,2);
x3 = bb(1,3):samp(3):bb(2,3);

return;
%=======================================================================

%=======================================================================
function [CP,BP,SP] = run_segment(CP,BP,SP,VF,sums,x1,x2,x3)
oll = -Inf;
spm_chi2_plot('Init','Segmenting','Log-likelihood','Iteration #');

for iter = 1:64,
	ll= 0;
	for pp = 1:length(x3), % Loop over planes
		bf        = get_bp(BP,x1,x2,x3(pp));
		[raw,msk] = get_raw(VF,x1,x2,x3(pp));
		s         = get_sp(SP,x1,x2,x3(pp));
		cor       = bf.*raw;
		[P,ll0]   = get_p(cor,msk,s,sums,CP,bf);
		ll        = ll + ll0;
		CP        = update_cp_est(CP,P,cor,msk,pp);
		BP        = update_bp_est(BP,P,cor,CP,msk,x1,x2,x3(pp));
	end;

	BP = update_bp(BP);
	if iter>1, spm_chi2_plot('Set',ll); end;
	%fprintf('\t%g\n', ll);

	% Stopping criterion
	%-----------------------------------------------------------------------
	if iter == 2,
		ll2 = ll;
	elseif iter > 2 & abs((ll-oll)/(ll-ll2)) < 0.0001
		break;
	end;
	oll = ll;
end;
spm_chi2_plot('Clear');
return;
%=======================================================================

%=======================================================================
function BP = init_bp(VF,co,reg)
m        = length(VF);
tmp      = sqrt(sum(VF(1).mat(1:3,1:3).^2));
BP.nbas  = max(round((VF(1).dim(1:3).*tmp)/co),[1 1 1]);
BP.B1    = spm_dctmtx(VF(1).dim(1),BP.nbas(1));
BP.B2    = spm_dctmtx(VF(1).dim(2),BP.nbas(2));
BP.B3    = spm_dctmtx(VF(1).dim(3),BP.nbas(3));

nbas     = BP.nbas;
if prod(BP.nbas)>1,
	% Set up a priori covariance matrix
	vx = sqrt(sum(VF(1).mat(1:3,1:3).^2));
	kx=(pi*((1:nbas(1))'-1)*pi/vx(1)/VF(1).dim(1)*10).^2;
	ky=(pi*((1:nbas(2))'-1)*pi/vx(2)/VF(1).dim(2)*10).^2;
	kz=(pi*((1:nbas(3))'-1)*pi/vx(3)/VF(1).dim(3)*10).^2;

	% Cost function based on sum of squares of 4th derivatives
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

	%IC0(1) = max(IC0);
	BP.IC0 = diag(IC0(2:end));

	% Initial estimate for intensity modulation field
	BP.T   = zeros(nbas(1),nbas(2),nbas(3),length(VF));
	%-----------------------------------------------------------------------
else
	BP.T   = zeros([1 1 1 length(VF)]);
	BP.IC0 = [];
end;
BP.Alpha = zeros(prod(BP.nbas(1:3)),prod(BP.nbas(1:3)),m);
BP.Beta  = zeros(prod(BP.nbas(1:3)),m);
return;
%=======================================================================

%=======================================================================
function BP = update_bp_est(BP,p,cor,CP,msk,x1,x2,x3)
if prod(BP.nbas)<=1, return; end;
B1 = BP.B1(x1,:);
B2 = BP.B2(x2,:);
B3 = BP.B3(x3,:);
for j=1:size(BP.Alpha,3),
	cr  = cor(:,:,j);
	w1 = zeros(size(cr));
	w2 = zeros(size(cr));
	for i=[1 2 3 4 5 6 7 8],
		tmp = p(:,:,i)*CP.cv(j,j,i)^(-1);
		w1  = w1 + tmp.*(CP.mn(j,i) - cr);
		w2  = w2 + tmp;
	end;
	wt1       = 1 + cr.*w1;
	wt2       = cr.*(cr.*w2 - w1);
	wt1(~msk) = 0;
	wt2(~msk) = 0;

	BP.Beta(:,j)    = BP.Beta(:,j)    + kron(B3',spm_krutil(wt1,B1,B2,0));
	BP.Alpha(:,:,j) = BP.Alpha(:,:,j) + kron(B3'*B3,spm_krutil(wt2,B1,B2,1));
end;
return;
%=======================================================================
 
%=======================================================================
function BP = update_bp(BP)
if prod(BP.nbas)<=1, return; end;
for j=1:size(BP.Alpha,3),
	x     = BP.T(:,:,:,j);
	x     = x(:);
	x     = x(2:end);
	Alpha = BP.Alpha(2:end,2:end,j);
	Beta  = BP.Beta(2:end,j);
	x     = (Alpha + BP.IC0)\(Alpha*x + Beta);

	BP.T(:,:,:,j) = reshape([0 ; x],BP.nbas(1:3));
	BP.Alpha      = zeros(size(BP.Alpha));
	BP.Beta       = zeros(size(BP.Beta));
end;
return;
%=======================================================================
 
%=======================================================================
function bf = get_bp(BP,x1,x2,x3)
bf = ones(length(x1),length(x2),size(BP.Alpha,3));
if prod(BP.nbas)<=1, return; end;
B1 = BP.B1(x1,:);
B2 = BP.B2(x2,:);
B3 = BP.B3(x3,:);
for i=1:size(BP.Alpha,3),
	t = reshape(reshape(BP.T(:,:,:,i),...
		BP.nbas(1)*BP.nbas(2),BP.nbas(3))*B3', BP.nbas(1), BP.nbas(2));
	bf(:,:,i) = exp(B1*t*B2');
end;
return;
%=======================================================================
 
%=======================================================================
function [dat,msk] = get_raw(VF,x1,x2,x3)
[X1,X2,X3] = ndgrid(x1,x2,x3);
for i=1:length(VF),
	[Y1,Y2,Y3] = affine_transform(X1,X2,X3,VF(i).mat\VF(1).mat);
	dat(:,:,i) = spm_sample_vol(VF(i),Y1,Y2,Y3,1);
end;
msk = all(dat,3) & all(isfinite(double(dat)),3);
return;
%=======================================================================
 
%=======================================================================
function CP = init_cp(VF,x3)
n = 8;
m = length(VF);
p = length(x3);
CP.mom0 = zeros(1,n,p)+eps;
CP.mom1 = zeros(m,n,p);
CP.mom2 = zeros(m,m,n,p)+eps;

% Occasionally the dynamic range of the images is such that many voxels
% all have the same intensity.  Adding cv0 is an attempt to improve the
% stability of the algorithm if this occurs. The value 0.083 was obtained
% from var(rand(1000000,1)).  It prbably isn't the best way of doing
% things, but it appears to work.
CP.cv0 = zeros(m,m);
for i=1:m,
	if spm_type(VF(i).dim(4),'intt'),
		CP.cv0(i,i)=0.083*mean(VF(i).pinfo(1,:));
	end;
end;
return;
%=======================================================================
 
%=======================================================================
function CP = shake_cp(CP)
CP.mom0(:,5,:)   = CP.mom0(:,1,:);
CP.mom0(:,6,:)   = CP.mom0(:,2,:);
CP.mom0(:,7,:)   = CP.mom0(:,3,:);
CP.mom1(:,5,:)   = CP.mom1(:,1,:);
CP.mom1(:,6,:)   = CP.mom1(:,2,:);
CP.mom1(:,7,:)   = CP.mom1(:,3,:);
CP.mom1(:,8,:)   = 0;
CP.mom2(:,:,5,:) = CP.mom2(:,:,1,:);
CP.mom2(:,:,6,:) = CP.mom2(:,:,2,:);
CP.mom2(:,:,7,:) = CP.mom2(:,:,3,:);
return;
%=======================================================================
 
%=======================================================================
function CP = update_cp_est(CP,P,dat,msk,p)
m   = size(dat,3);
d   = size(P);
P   = reshape(P,[d(1)*d(2),d(3)]);
dat = reshape(dat,[d(1)*d(2),m]);
P(~msk(:),:)   = [];
dat(~msk(:),:) = [];
for i=1:size(CP.mom0,2),
	CP.mom0(1,i,p)   = sum(P(:,i));
	CP.mom1(:,i,p)   = sum((P(:,i)*ones(1,m)).*dat)';
	CP.mom2(:,:,i,p) = ((P(:,i)*ones(1,m)).*dat)'*dat;
end;

for i=1:size(CP.mom0,2),
	CP.mg(1,i)   = sum(CP.mom0(1,i,:),3);
	CP.mn(:,i)   = sum(CP.mom1(:,i,:),3)/CP.mg(1,i);

	tmp          = (CP.mg(1,i).*CP.mn(:,i))*CP.mn(:,i)';
	tmp          = tmp-eye(size(tmp))*eps*10000;
	CP.cv(:,:,i) = (sum(CP.mom2(:,:,i,:),4) - tmp)/CP.mg(1,i) + CP.cv0;
end;
CP.mg   = CP.mg/sum(CP.mg);
return;
%=======================================================================
 
%=======================================================================
function [p,ll] = get_p(cor,msk,s,sums,CP,bf)
d   = [size(cor) 1 1];
n   = size(CP.mg,2);
cor = reshape(cor,d(1)*d(2),d(3));
cor = cor(msk,:);
p   = zeros(d(1)*d(2),n);
if ~any(msk), p  = reshape(p,d(1),d(2),n); ll=0; return; end;

for i=1:n,
	amp       = 1/sqrt((2*pi)^d(3) * det(CP.cv(:,:,i)));
	dst       = (cor-ones(size(cor,1),1)*CP.mn(:,i)')/sqrtm(CP.cv(:,:,i));
	dst       = sum(dst.*dst,2);
	tmp       = s(:,:,i);
	p(msk,i)  = (amp*CP.mg(1,i)/sums(i))*exp(-0.5*dst).*tmp(msk) +eps;
end;
sp = sum(p,2);
ll = sum(log(sp(msk).*bf(msk)+eps));
sp(~msk) = Inf;
for i=1:n, p(:,i) = p(:,i)./sp; end;
p  = reshape(p,d(1),d(2),n);
return;
%=======================================================================
 
%=======================================================================
function SP = init_sp(flags,VF,PG)
SP.VB       = spm_vol(flags.priors);
MM          = get_affine_mapping(VF,PG,flags.affreg);
%VF          = spm_vol(PF);
SP.MM       = MM*VF(1).mat;
SP.w        = 0.98;
return;
%=======================================================================
 
%=======================================================================
function s = get_sp(SP,x1,x2,x3)
[X1,X2,X3] = ndgrid(x1,x2,x3);
[Y1,Y2,Y3] = affine_transform(X1,X2,X3,SP.VB(1).mat\SP.MM);
w1  = SP.w;
w2  = (1-w1)/2;
s   = zeros([size(Y1),4]);
for i=1:3,
	s(:,:,i) = spm_sample_vol(SP.VB(i),Y1,Y2,Y3,1)*w1+w2;
end;
s(:,:,4:8)   = repmat(abs(1-sum(s(:,:,1:3),3))/5,[1 1 5]);
return;
%=======================================================================
 
%=======================================================================
function [g,w,c] = get_gwc(VF,BP,SP,CP,sums,wc)

if wc,
	VC = VF;
	for j=1:length(VF),
		[pth,nm,xt]    = fileparts(deblank(VF(j).fname));
		VC(j).fname    = fullfile(pth,['m' nm xt]);
		VC(j).descrip  = 'Bias corrected image';
	end;
	VC = spm_create_vol(VC);
end;

spm_progress_bar('Init',VF(1).dim(3),'Creating Segmented','planes completed');
x1 = 1:VF(1).dim(1);
x2 = 1:VF(1).dim(2);
x3 = 1:VF(1).dim(3);

g = uint8(0); g(VF(1).dim(1),VF(1).dim(2),VF(1).dim(3)) = 0;
w = uint8(0); w(VF(1).dim(1),VF(1).dim(2),VF(1).dim(3)) = 0;
c = uint8(0); c(VF(1).dim(1),VF(1).dim(2),VF(1).dim(3)) = 0;

for pp=1:length(x3),
	bf        = get_bp(BP,x1,x2,x3(pp));
	[raw,msk] = get_raw(VF,x1,x2,x3(pp));
	cor       = raw.*bf;
	if wc,
		for j=1:length(VC),
			VC(j) = spm_write_plane(VC(j),cor(:,:,j),pp);
		end;
	end;
	s         = get_sp(SP,x1,x2,x3(pp));
	p         = get_p(cor,msk,s,sums,CP,bf);
	g(:,:,pp) = uint8(round(p(:,:,1)*255));
	w(:,:,pp) = uint8(round(p(:,:,2)*255));
	c(:,:,pp) = uint8(round(p(:,:,3)*255));

	spm_progress_bar('Set',pp);
end;
spm_progress_bar('Clear');

if wc, spm_close_vol(VC); end;
return;
%=======================================================================
 
%=======================================================================
function [g,w,c] = clean_gwc(g,w,c)
b    = w;
b(1) = w(1);

% Build a 3x3x3 seperable smoothing kernel
%-----------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Erosions and conditional dilations
%-----------------------------------------------------------------------
niter = 32;
spm_progress_bar('Init',niter,'Extracting Brain','Iterations completed');
for j=1:niter,
	if j>2, th=0.15; else th=0.6; end; % Dilate after two its of erosion.
	for i=1:size(b,3),
		gp = double(g(:,:,i));
		wp = double(w(:,:,i));
		bp = double(b(:,:,i))/255;
		bp = (bp>th).*(wp+gp);
		b(:,:,i) = uint8(round(bp));
	end;
	spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
	spm_progress_bar('Set',j);
end;
th = 0.05;
for i=1:size(b,3),
	gp       = double(g(:,:,i))/255;
	wp       = double(w(:,:,i))/255;
	cp       = double(c(:,:,i))/255;
	bp       = double(b(:,:,i))/255;
	bp       = ((bp>th).*(wp+gp))>th;
	g(:,:,i) = uint8(round(255*gp.*bp./(gp+wp+cp+eps)));
	w(:,:,i) = uint8(round(255*wp.*bp./(gp+wp+cp+eps)));
	c(:,:,i) = uint8(round(255*(cp.*bp./(gp+wp+cp+eps)+cp.*(1-bp))));
end;
spm_progress_bar('Clear');
return;
