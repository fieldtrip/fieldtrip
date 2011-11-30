function img=surf2volz(node,face,xi,yi,zi)
%
% img=surf2volz(node,face,xi,yi,zi)
%
% convert a triangular surface to a shell of voxels in a 3D image
% along the z-axis
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 node: node list of the triangular surface, 3 columns for x/y/z
%	 face: triangle node indices, each row is a triangle
%	 xi,yi,zi: x/y/z grid for the resulting volume
%
% output:
%	 img: a volumetric binary image at position of ndgrid(xi,yi,zi)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

ne=size(face,1);
img=zeros(length(xi),length(yi),length(zi),'uint8');

dx0=min(abs(diff(xi)));
dx=dx0/2;
dy0=min(abs(diff(yi)));
dy=dy0/2;
dz0=min(abs(diff(zi)));

dl=sqrt(dx*dx+dy*dy);

minz=min(node(:,3));
maxz=max(node(:,3));
iz=hist([minz,maxz],zi);
hz=find(iz);
iz=hz(1):min(length(zi),hz(end)+1);

for i=1:length(iz);
    plane=[0 100 zi(iz(i)); 100 0 zi(iz(i)); 0 0 zi(iz(i))];
    [bcutpos,bcutvalue,bcutedges]=qmeshcut(face(:,1:3),node,node(:,1),plane);
    if(isempty(bcutpos))
        continue;
    end
    enum=length(bcutedges);
    for j=1:enum
        e0=bcutpos(bcutedges(j,1),1:2);
        e1=bcutpos(bcutedges(j,2),1:2);
        len=ceil(sum(abs(e1-e0))/(abs(dx)+abs(dy)))+1;
        dd=(e1-e0)/len;

        posx= floor(((e0(1)+(0:len)*dd(1)-xi(1)))/dx0)';
        posy= floor(((e0(2)+(0:len)*dd(2)-yi(1)))/dy0)';
        pos=[posx, posy];
        pos(find(posx>length(xi) | posy>length(yi) | posx<=0|posy<=0), :)=[];

        if(length(pos)>0)
            zz=floor(((zi(iz(i))-zi(1)))/dy0);
            for k=1:size(pos,1)
                img(pos(k,1),pos(k,2),zz)=1;
            end
          %img(sub2ind(size(img),pos(:,1),pos(:,2),i*ones(size(pos,1),1)))=1;
        end
    end
end

