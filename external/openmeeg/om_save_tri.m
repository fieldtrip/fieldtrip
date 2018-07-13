function [] = om_save_tri(filename,points,faces,nrmls)

% OM_SAVE_TRI   Save .tri file
%
%   SYNTAX
%       [] = OM_SAVE_TRI(FILENAME,POINTS,FACES,NORMALS)
%

% Copyright (C) 2010-2017, OpenMEEG developers

% OpenMEEG vs 2.3 and up internally adjusts the convention for surface
% normals, but OpenMEEG v2.2 expects surface normals to point inwards;
% this checks and corrects if needed
ok = checknormals(points,faces);
if ~ok
    faces = fliplr(faces);
end

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

function ok = checknormals(points,faces)
% FIXME: this method is rigorous only for star shaped surfaces
ok = 0;
% translate to the center
org = mean(points,1);
points(:,1) = points(:,1) - org(1);
points(:,2) = points(:,2) - org(2);
points(:,3) = points(:,3) - org(3);

w = sum(solid_angle(points, faces));

if w<0 && (abs(w)-4*pi)<1000*eps
  ok = 0;
  disp('Your surface normals are outwards oriented, flipping for OpenMEEG')
elseif w>0 && (abs(w)-4*pi)<1000*eps
  ok = 1;
%   ft_warning('your normals are inwards oriented')
else
  ft_error('your surface probably is irregular\n')
  ok = 0;
end
end %   function
