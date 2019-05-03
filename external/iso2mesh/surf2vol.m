function [img, v2smap]=surf2vol(node,face,xi,yi,zi,varargin)
%
% [img, v2smap]=surf2vol(node,face,xi,yi,zi,'options',values,...)
%
% convert a triangular surface to a shell of voxels in a 3D image
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 node: node list of the triangular surface, 3 columns for x/y/z
%	 face: triangle node indices, each row is a triangle
%              if face contains the 4th column, it indicates the label of
%              the face triangles (each face componment must be closed); if
%              face contains 5 columns, it stores a tetrahedral mesh with
%              labels, where the first 4 columns are the element list and 
%              the last column is the element label;
%	 xi,yi,zi: x/y/z grid for the resulting volume
%        options: 'fill', if set to 1, the enclosed voxels are labeled by 1
%                 'label', if set to 1, the enclosed voxels are labeled by
%                          the corresponding label of the face or element;
%                          setting 'label' to 1 also implies 'fill'.
%
% output:
%	 img: a volumetric binary image at position of ndgrid(xi,yi,zi)
%        v2smap (optional): a 4x4 matrix denoting the Affine transformation to map
%             the voxel coordinates back to the mesh space.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fprintf(1,'converting a closed surface to a volumetric binary image ...\n');

opt=varargin2struct(varargin{:});
label=jsonopt('label',0,opt);

elabel=[1];
if(size(face,2)>=4)
        elabel=unique(face(:,end));
        if(size(face,2)==5)
            label=1;
            el=face;
            face=[];
            for i=1:length(elabel)
                    fc=volface(el(el(:,5)==elabel(i),1:4));
                    fc(:,4)=elabel(i);
                    face=[face ; fc];
            end
        end
else
        fc=face;
end

img=zeros(length(xi),length(yi),length(zi));

for i=1:length(elabel)
        if(size(face,2)==4)
            fc=face(face(:,4)==elabel(i),1:3);
        end
        im=surf2volz(node(:,1:3),fc(:,1:3),xi,yi,zi);
        im=im | shiftdim(surf2volz(node(:,[3 1 2]),fc(:,1:3),zi,xi,yi),1);
        im=im | shiftdim(surf2volz(node(:,[2 3 1]),fc(:,1:3),yi,zi,xi),2);

        v2smap=[];

        % here we assume the grid is uniform; surf2vol can handle non-uniform grid, 
        % but the affine output is not correct in this case

        if(jsonopt('fill',0,opt) || label)
                im=imfill(im,'holes');
                if(label)
                    im=im*elabel(i);
                end
        end
        img=max(im,img);
end
 
if(nargout>1) 
        dlen=abs([xi(2)-xi(1) yi(2)-yi(1) zi(2)-zi(1)]);
        p0=min(node);
        offset=p0;
        v2smap=diag(abs(dlen));
        v2smap(4,4)=1;
        v2smap(1:3,4)=offset';
end