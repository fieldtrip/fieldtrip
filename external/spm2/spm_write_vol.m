function V = spm_write_vol(V,Y)
% Write an image volume to disk, setting scales and offsets as appropriate
% FORMAT V = spm_write_vol(V,Y)
% V (input)  - a structure containing image volume information (see spm_vol)
% Y          - a one, two or three dimensional matrix containing the image voxels
% V (output) - data structure after modification for writing.
%_______________________________________________________________________
% @(#)spm_write_vol.m	2.9 John Ashburner 03/02/26

if ndims(Y)>3, error('Can only handle a maximum of 3 dimensions.'), end

if ~isfield(V,'pinfo'), V.pinfo = [1,0,0]'; end

dim = [size(Y) 1 1 1];
if ~all(dim(1:3) == V.dim(1:3)) | (size(V.pinfo,2)~=1 & size(V.pinfo,2)~=dim(3)),
	error('Incompatible dimensions.');
end


% Set scalefactors and offsets
%-----------------------------------------------------------------------
dt = V.dim(4); if dt>256, dt = dt/256; end;
if any(dt == [128+2 128+4 128+8]),
	% Convert to a form that Analyze will support
	dt = dt - 128;
end;
s            = find(dt == [2 4 8 128+2 128+4 128+8]);
dmnmx        = [0 -2^15 -2^31 -2^7 0 0 ; 2^8-1 2^15-1 2^31-1 2^7-1 2^16 2^32];
dmnmx        = dmnmx(:,s);
V.pinfo(1,:) = 1;
V.pinfo(2,:) = 0;
mxs          = zeros(dim(3),1)+NaN;
mns          = zeros(dim(3),1)+NaN;
if ~isempty(s),
	for p=1:dim(3),
		tmp    = double(Y(:,:,p));
		tmp    = tmp(isfinite(tmp));
		if ~isempty(tmp),
			mxs(p) = max(tmp);
			mns(p) = min(tmp);
		end;
	end;

	if size(V.pinfo,2) ~= 1,
		for p=1:dim(3),
			mx = mxs(p);
			mn = mns(p);
			if ~isfinite(mx), mx = 0; end;
			if ~isfinite(mn), mn = 0; end;
			if mx~=mn,
				V.pinfo(1,p) = (mx-mn)/(dmnmx(2)-dmnmx(1));
				V.pinfo(2,p) = ...
					(dmnmx(2)*mn-dmnmx(1)*mx)/(dmnmx(2)-dmnmx(1));
			else,
				V.pinfo(1,p) = 0;
				V.pinfo(2,p) = mx;
			end;
		end;
	else,
		mx = max(mxs(isfinite(mxs)));
		mn = min(mns(isfinite(mns)));
		if isempty(mx), mx = 0; end;
		if isempty(mn), mn = 0; end;
		if mx~=mn,
			V.pinfo(1,1) = (mx-mn)/(dmnmx(2)-dmnmx(1));
			V.pinfo(2,1) = (dmnmx(2)*mn-dmnmx(1)*mx)/(dmnmx(2)-dmnmx(1));
		else,
			V.pinfo(1,1) = 0;
			V.pinfo(2,1) = mx;
		end;
	end;
end;

%-Create and write image
%-----------------------------------------------------------------------
V = spm_create_vol(V);
for p=1:V.dim(3),
	V = spm_write_plane(V,Y(:,:,p),p);
end;
V = spm_close_vol(V);
