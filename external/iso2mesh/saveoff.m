function saveoff(v,f,fname)
%
% saveoff(v,f,fname)
%
% save a surface mesh to Geomview Object File Format (OFF)
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/03/28
%
% input:
%      v: input, surface node list, dimension (nn,3)
%      f: input, surface face element list, dimension (be,3)
%      fname: output file name
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fid=fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end
fprintf(fid,'OFF\n');
fprintf(fid,'%d\t%d\t%d\n',length(v),length(f),0);
fprintf(fid,'%f\t%f\t%f\n',v');
face=[size(f,2)*ones(size(f,1),1) f-1];
format=[repmat('%d\t',1,size(face,2)-1) '%d\n'];
fprintf(fid,format,face');
fclose(fid);

