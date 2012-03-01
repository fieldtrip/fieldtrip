function savebinstl(node,elem,fname,solidname)
%
% savebinstl(node,elem,fname,solidname)
%
% save a tetrahedral mesh to a binary STL (Standard Tessellation Language) file
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2012/02/14
%
% input:
%      node: input, surface node list, dimension Nx3
%      elem: input, tetrahedral element list; if size(elem,2)==3, it is a surface
%      fname: output file name
%      solidname: an optional string for the name of the object
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fid=fopen(fname,'wb');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end

if(nargin<3)
    error('incomplete input');
end
if(isempty(node) || size(node,2)<3)
    error('invalid node input');
end

hddata=struct('Ver',1,'Creator','iso2mesh','Date',datestr(now));
if(nargin<4) 
    solidname='';
else
    hddata.name=solidname;
end

headerstr=savejson('',hddata);
headerstr=regexprep(headerstr,'[\t\n\r]','');
headerstr=regexprep(headerstr,'": ','":');
if(length(headerstr)>80)
    headerstr=headerstr(1:80); % this makes the header an invalid json string
else
    headerstr(81)=0;
end

fwrite(fid,headerstr(1:80),'char');

if(~isempty(elem))
  if(size(elem,2)>=5)
	elem(:,5:end)=[];
  end
  if(size(elem,2)==4)
	elem=meshreorient(node,elem);
	elem=volface(elem);
  end
  ev=surfplane(node,elem);
  ev=ev(:,1:3)./repmat(sqrt(sum(ev(:,1:3).*ev(:,1:3),2)),1,3);

  len=size(elem,1);
  fwrite(fid,len,'uint32');

  for i=1:len
      fwrite(fid,[ev(i,:)',node(elem(i,:),1:3)'],'float32');
      fwrite(fid,0,'uint16');
  end
end

fclose(fid);
