function savesmf(v,f,fname)
%
% savesmf(v,f,fname)
%
% save a surface mesh to smf format
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/21
%
% input:
%      v: input, surface node list, dimension (nn,3)
%      f: input, surface face element list, dimension (be,3)
%      fname: output file name
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fid=fopen(fname,'wt');
fprintf(fid,'v %f %f %f\n',v');
fprintf(fid,'f %d %d %d\n',f');
fclose(fid);
