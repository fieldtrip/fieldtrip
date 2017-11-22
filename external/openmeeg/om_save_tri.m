function [] = om_save_tri(filename,points,faces,nrmls)

% OM_SAVE_TRI   Save .tri file
%
%   SYNTAX
%       [] = OM_SAVE_TRI(FILENAME,POINTS,FACES,NORMALS)
%

% Copyright (C) 2010-2017, OpenMEEG developers

if nargin<4 || isempty(nrmls)
    nrmls = om_normals(points,faces);
end

fid = fopen(filename,'w');
npoints = size(points,1);
fprintf(fid,'- %g\n',npoints);
fprintf(fid,'%g %g %g %g %g %g\n',[points , nrmls]');
nfaces = size(faces,1);
faces = faces-1;
fprintf(fid,'- %g %g %g\n', [nfaces nfaces nfaces]);
fprintf(fid,'%g %g %g\n',faces');

fclose(fid);

end %  function
