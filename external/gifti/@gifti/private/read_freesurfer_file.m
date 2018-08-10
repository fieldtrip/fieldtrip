function this = read_freesurfer_file(filename)
% Low level reader of FreeSurfer files (ASCII triangle surface file)
% FORMAT this = read_freesurfer_file(filename)
% filename    - FreeSurfer file
%
% See http://wideman-one.com/gw/brain/fs/surfacefileformats.htm
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: read_freesurfer_file.m 5322 2013-03-13 15:04:14Z guillaume $


fid = fopen(filename,'rt');
if fid == -1, error('Cannot open "%s".',filename); end

fgetl(fid); % #!ascii 
N = fscanf(fid,'%d',2);
this.vertices = fscanf(fid,'%f %f %f %d',[4 N(1)])';
this.faces    = fscanf(fid,'%d %d %d %d',[4 N(2)])';

fclose(fid);

this.vertices = this.vertices(:,1:3);
this.faces    = this.faces(:,1:3) + 1;
