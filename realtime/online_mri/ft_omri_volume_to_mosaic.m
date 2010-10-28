function M = ft_omri_volume_to_mosaic(V)
% function M = ft_omri_volume_to_mosaic(V)
% 
% Reshuffle [MxNxS] volume to [XxY] mosaic

% (C) 2010 S.Klanke

[h,w,ns] = size(V);
nn = ceil(sqrt(ns));
mm = ceil(ns/nn);

M = zeros(w*mm, h*nn, class(V));
row = 1;
col = 1;
for n = 1:ns
	is = 1 + (row-1)*h;
	ie = row*h;
	js = 1 + (col-1)*w;
	je = col*w;
		
	M(is:ie, js:je) = rot90(V(:,:,n));
	col = col + 1;
	if col>nn
	   row = row + 1;
	   col = 1;
	end
end