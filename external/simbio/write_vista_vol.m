function write_vista_vol(dim,seg,filename)

% WRITE_VISTA_VOL is substitute matlab version for the mex file write_vista_vol.cpp
% 
% This matlab version might not contain all the features of the cpp file
% and is provided only in case the source .cpp file does not compile.
%
% with:
% dim = size(seg)   dimensions of the volume
% seg               volume, a set of images (like MR or CT)
% filename          the name to be saved (with extension .v)
%
% $Id$

% number of points
N = prod(dim);
nbands  = size(seg,3);
nframes = nbands;
nrows   = size(seg,2); 
ncolumns= size(seg,1);
VFileDelimiter = '\f\n';

% write header
try
  fid = fopen(filename,'w');
  fprintf(fid,['V-data 2 {\n']);
  fprintf(fid,['	image: image {\n']);
  fprintf(fid,['		data: 0\n']);
  fprintf(fid,['		length: ' num2str(N) '\n']);
  fprintf(fid,['		nbands: ' num2str(nbands) '\n']);
  fprintf(fid,['		nframes: ' num2str(nframes) '\n']);
  fprintf(fid,['		nrows: ' num2str(nrows) '\n']);
  fprintf(fid,['		ncolumns: ' num2str(ncolumns) '\n']);
  fprintf(fid,['		repn: ubyte\n']);
  fprintf(fid,['		voxel: "1.000000 1.000000 1.000000"\n']);
  fprintf(fid,['	}\n']);
  fprintf(fid,['}']);
  fprintf(fid,['\n' VFileDelimiter]);
  fwrite(fid,cast(seg,'uint8'),'uint8');
  fclose(fid);
end

