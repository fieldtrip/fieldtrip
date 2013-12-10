function dat=readinr(fname)
%
% vol=readinr(fname)
%
% load a volume from an INR file
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2009/05/03
%
% input:
%      fname: input file name
%
% output:
%      dat: output, data read from the inr file
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fid=fopen(fname,'rb');
s=fread(fid,256,'uchar');

s=char(s)';

if(regexp(s,'#INRIMAGE-4')~=1)
   error('INRIMAGE header was not found')
end

nx=regexp(s,'XDIM\s*=\s*([0-9]+)','tokens');
if(length(nx)) 
   nx=str2num(nx{1}{1});
else
   error('no XDIM found');
end

ny=regexp(s,'YDIM\s*=\s*([0-9]+)','tokens');
if(length(ny)) 
   ny=str2num(ny{1}{1});
else
   error('no YDIM found');
end

nz=regexp(s,'ZDIM\s*=\s*([0-9]+)','tokens');
if(length(nz)) 
   nz=str2num(nz{1}{1});
else
   error('no ZDIM found');
end

nv=regexp(s,'VDIM\s*=\s*([0-9]+)','tokens');
if(length(nv))
   nv=str2num(nv{1}{1});
else
   nv=1;
end

type=regexp(s,'TYPE=([a-z ]+)','tokens');
if(length(type))
   type=type{1}{1};
else
   error('no TYPE found');
end

pixel=regexp(s,'PIXSIZE=([0-9]+)','tokens');
if(length(pixel))
   pixel=str2num(pixel{1}{1});
else
   error('no PIXSIZE found');
end

%header=sprintf(['#INRIMAGE-4#{\nXDIM=%d\nYDIM=%d\nZDIM=%d\nVDIM=1\nTYPE=%s\n' ...
%  'PIXSIZE=%d bits\nCPU=decm\nVX=1\nVY=1\nVZ=1\n'],size(vol),btype,bitlen);

if(strcmp(type,'unsigned fixed') & pixel==8)
   dtype='uint8';
elseif(strcmp(type,'unsigned fixed') & pixel==16)
   dtype='uint16';
elseif(strcmp(type,'float') & pixel==32)
   dtype='float';
elseif(strcmp(type,'float') & pixel==64)
   dtype='double';
else
   error('volume format not supported');
end


dat=fread(fid,nx*ny*nz*nv,dtype);
fclose(fid);

if(nv==1)
   dat=reshape(dat,[nx,ny,nz]);
else
   dat=reshape(dat,[nx,ny,nz,nv]);
end
