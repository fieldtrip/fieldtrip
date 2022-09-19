function [model, mat_r, mat_a, Va] = ft_omri_align_scan(model, Vo)

% function [model, mat_r, mat_a, Va] = ft_omri_align_scan(model, Vo)
%
% Estimate 6-DOF motion parameters and align volume Vo with reference model.
%
% INPUTS
% model - data structure containing reference image and flags
% Vo    - original image (to be registered to the reference model
%
% OUTPUTS
% model - same as input, but with modified voxel mask (for keeping track of "missing voxels")
% Va    - aligned image (possibly rotated + translated version of Vo)
% mat_r - parameters of 6-dof transformation in homogenous matrix form
%         transformation from world coordinates of Vo to world coordinates of reference
% mat_a - pixel to world coordinate matrix of Vo (absolute "position" of Vo)

% Basic code taken from SPM8, adapted for online use 2010, S.Klanke

startTime = tic;

%-----------------------------------------------------------------------
V  = smooth_vol(Vo, model.interp, model.wrap, model.fwhm, model.mat);
d  = [size(V) 1 1];
d  = d(1:3);
ss = Inf;
countdown = -1;

mat_a = model.mat;

for iter=1:64,
	[y1,y2,y3] = coords([0 0 0  0 0 0],model.mat, mat_a, model.x1, model.x2, model.x3);
	msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
	if length(msk)<32
		error 'Too little information left';
	end

	F          = spm_bsplins(V, y1(msk),y2(msk),y3(msk),model.deg);
	if ~isempty(model.wt), F = F.*model.wt(msk); end

	A          = model.A0(msk,:);
	b1         = model.b(msk);
	sc         = sum(b1)/sum(F);
	b1         = b1-F*sc;
	soln       = (A'*A)\(A'*b1);

	p          = [0 0 0  0 0 0  1 1 1  0 0 0];
	p(model.lkp)     = p(model.lkp) + soln';
	mat_a      = inv(spm_matrix(p))*mat_a;

	pss        = ss;
	ss         = sum(b1.^2)/length(b1);
	
	%fprintf(1,'[%3i] %f\n',iter,ss);
	
	if (pss-ss)/pss < 1e-8 && countdown == -1, % Stopped converging.
		countdown = 2;
	end;
	if countdown ~= -1,
		if countdown==0, break; end
		countdown = countdown -1;
	end;
	if toc(startTime) > model.time
		warning 'Timeout for realignment - returning suboptimal solution';
		break;
	end
end

mat_r = mat_a/model.mat;

if nargout > 3
	Mpix  = mat_a\model.mat;
	%Va = kspace3d(Vo, Mpix);
	[Va, mask] = reslice_vol(Vo, Mpix, model.interp);
	model.mask = model.mask & mask;
end
