function facecell=finddisconnsurf(f)
%
% facecell=finddisconnsurf(f)
% 
% subroutine to extract disconnected surfaces from a 
% cluster of surfaces
% 
% author: Qianqian Fang (fangq@nmr.mgh.harvard.edu)
% Date: 2008/03/06
%
% input: 
%     f: faces defined by node indices for all surface triangles
%
% output:
%     facecell: separated disconnected surface node indices
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

facecell={};
subset=[];
while(length(f))
	idx=reshape(ismember(f,f(1,:)), size(f));
	[ii,jj]=find(sum(idx')');
	while(length(ii)>0)
		if(isempty(ii)) break; end
		%ii=unique(ii);
		subset(end+1:end+length(ii),:)=f(ii,:);
		f(ii,:)=[];
		idx=reshape(ismember(f,subset), size(f));
	        [ii,jj]=find(sum(idx')');
	end
	if(length(subset))
		facecell{end+1}=subset;
        	subset=[];
	end
end
