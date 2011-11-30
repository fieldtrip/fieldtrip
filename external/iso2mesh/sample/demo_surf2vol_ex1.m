%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   demo script to convert a closed surface to a binary image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preparation
% user must add the path of iso2mesh to matlab path list
% addpath('../');

%% load the sample data
load rat_head.mat

% first, generate a surface from the original image
% similar to demo_shortcuts_ex1.m

[node,face,regions,holes]=v2s(volimage,0.5,3);

node=sms(node,face(:,1:3),3,0.5); % apply 3 mesh smoothing

mdim=ceil(max(node)+1);
dstep=0.25;
zslice=15;
xrange=0:dstep:mdim(1);
yrange=0:dstep:mdim(2);
zrange=0:dstep:mdim(3);
img=surf2vol(node,face(:,1:3),xrange,yrange,zrange);

imagesc(squeeze(img(:,:,zslice))); % z=10

hold on

z0=zslice*dstep;
plane=[min(node(:,1)) min(node(:,2)) z0
       min(node(:,1)) max(node(:,2)) z0
       max(node(:,1)) min(node(:,2)) z0];

% run qmeshcut to get the cross-section information at z=mean(node(:,1))
% use the x-coordinates as the nodal values

[bcutpos,bcutvalue,bcutedges]=qmeshcut(face(:,1:3),node,node(:,1),plane);
[bcutpos,bcutedges]=removedupnodes(bcutpos,bcutedges);
bcutloop=extractloops(bcutedges);
bcutloop(isnan(bcutloop))=[]; % there can be multiple loops, remove the separators
plot(bcutpos(bcutloop,2)*(1/dstep),bcutpos(bcutloop,1)*(1/dstep),'w');

if(isoctavemesh)
  if(~exist('bwfill'))
    error('you need to install octave-image toolbox first');
  end
  img2=zeros(size(img),'uint8');
  for i=1:size(img,3)
    img2(:,:,i)=bwfill(img(:,:,i),'holes');
  end
  img2=img2+img;
else
  img2=imfill(img,'holes')+img;
end
figure;
imagesc(squeeze(img2(:,:,zslice))); % z=10
hold on;
plot(bcutpos(bcutloop,2)*(1/dstep),bcutpos(bcutloop,1)*(1/dstep),'y--');
