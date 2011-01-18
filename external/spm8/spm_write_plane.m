function V = spm_write_plane(V,dat,n)
% Write transverse plane(s) of image data.
% FORMAT VO = spm_write_plane(V,dat,n)
% V   - data structure containing image information.
%       - see spm_vol for a description.
% dat - the two/three dimensional image to write.
% n   - the plane number(s) (beginning from 1). If an entire volume
%       should be written, n should contain the single character ':'
%       instead of plane numbers.
%
% VO  - (possibly) modified data structure containing image information.
%       It is possible that future versions of spm_write_plane may
%       modify scalefactors (for example).
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_write_plane.m 1143 2008-02-07 19:33:33Z spm $

% For performance reasons, on network filesystems one should write
% out as large contiguous blocks data as possible at once. Therefore,
% multiple planes or even entire volumes should be handled here. 
% Dimension checking is left to mat2file. 

if isfield(V,'n'), 
    n1= num2cell(V.n);
    n = {n n1{:}}; 
else 
    n = {n}; 
end;
S     = struct('type','()','subs',{{':',':',n{:}}});
V.private.dat = subsasgn(V.private.dat,S,dat);
