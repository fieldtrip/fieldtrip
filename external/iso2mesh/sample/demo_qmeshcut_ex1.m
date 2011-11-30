% qmeshcut demonstration
%
% by Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% to demonstrate how to use qmeshcut to produce cross-sectional plot 
% of an un-structured (tetrahedral) mesh

% run vol2mesh demo 1 to create a 3d mesh

demo_vol2mesh_ex1

% define a plane by 3 points, in this case, z=mean(node(:,3))

z0=mean(node(:,3));

plane=[min(node(:,1)) min(node(:,2)) z0
       min(node(:,1)) max(node(:,2)) z0
       max(node(:,1)) min(node(:,2)) z0];

% run qmeshcut to get the cross-section information at z=mean(node(:,1))
% use the x-coordinates as the nodal values

[cutpos,cutvalue,facedata]=qmeshcut(elem(:,1:4),node,node(:,1),plane);

% plot your results

figure;
hsurf=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'facecolor','none');
hold on;
if(isoctavemesh)
  hcut=patch('Faces',facedata,'Vertices',cutpos);
else
  hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','interp');
end
%set(hcut, 'linestyle','none')
axis equal;

% qmeshcut can also cut a surface

[bcutpos,bcutvalue,bcutedges]=qmeshcut(face(:,1:3),node,node(:,1),plane);
[bcutpos,bcutedges]=removedupnodes(bcutpos,bcutedges);
bcutloop=extractloops(bcutedges);

bcutloop(isnan(bcutloop))=[]; % there can be multiple loops, remove the separators

% plot the plane-surface cuts

plot3(bcutpos(bcutloop,1),bcutpos(bcutloop,2),bcutpos(bcutloop,3),'r','LineWidth',4);

% essencially, this should be the same as you do a removedupnodes(cutpos,facedata)
% and then call extractloop(facedata)
