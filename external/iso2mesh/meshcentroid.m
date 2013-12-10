function centroid=meshcentroid(v,f)
%
% centroid=meshcentroid(v,f)
% 
% compute the centroids of a mesh defined by nodes and elements
% (surface or tetrahedra) in R^n space
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%      v: surface node list, dimension (nn,3)
%      f: surface face element list, dimension (be,3)
%
% output:
%      centroid: centroid positions, one row for each element
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(~iscell(f))
    ec=reshape(v(f(:,1:size(f,2))',:)', [size(v,2) size(f,2) size(f,1)]);
    centroid=squeeze(mean(ec,2))';
else
    len=length(f);
    centroid=zeros(len,size(v,2));
    try
        for i=1:len
            fc=f{i};
            if(~isempty(fc)) % need to set centroid to NaN if fc is empty?
                vlist=fc{1};
                centroid(i,:)=mean(v(vlist(find(~isnan(vlist))),:));
            end
        end
    catch
        error('malformed face cell array');
    end
end

