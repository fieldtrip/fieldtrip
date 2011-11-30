%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   demo script for surface smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preparation
% user must add the path of iso2mesh to matlab path list
% addpath('../');

% user need to add the full path to .../iso2mesh/bin directory
% to windows/Linux/Unix PATH environment variable

%% load the sample data
load rat_head.mat

% volimage is a volumetric image such as an X-ray or MRI image
% A,b are registration matrix and vector, respectively
%% perform mesh generation

[node,face]=v2s(volimage,0.5,2,'cgalmesh');

face=face(:,1:3);

p0=min(node);
p1=max(node);

rownum=3;
colnum=4;
figure;
subplot(rownum,colnum,1);

plotmesh(node,face(:,1:3));
if(~isoctavemesh) 
	title({'Laplacian+HC Smoothing Test','no smoothing'}); 
else
	title('Laplacian+HC - no smoothing');
end
axis equal;
set(gca,'xlim',[p0(1),p1(1)],'ylim',[p0(2),p1(2)],'zlim',[p0(3),p1(3)])

%=========================================================
% apply Laplacian+HC smoothing
%=========================================================

n1=node;
for i=1:rownum*colnum-1
  n1=sms(n1,face(:,1:3),1,0.5); % apply Laplacian+HC mesh smoothing
  subplot(rownum,colnum,i+1);
  plotmesh(n1,face(:,1:3));
  title(['iter=' num2str(i)]);
  axis equal;
  set(gca,'xlim',[p0(1),p1(1)],'ylim',[p0(2),p1(2)],'zlim',[p0(3),p1(3)])
end

%=========================================================
% apply Laplacian smoothing
%=========================================================

figure;
subplot(rownum,colnum,1);

plotmesh(node,face(:,1:3));
if(~isoctavemesh)
        title({'Laplacian Smoothing Test','no smoothing'});
else
        title('Laplacian - no smoothing');
end
axis equal;
set(gca,'xlim',[p0(1),p1(1)],'ylim',[p0(2),p1(2)],'zlim',[p0(3),p1(3)])

conn=meshconn(face(:,1:3),size(node,1));

n1=node;
for i=1:rownum*colnum-1
  n1=smoothsurf(n1,[],conn,1,0.5,'laplacian');
  subplot(rownum,colnum,i+1);
  plotmesh(n1,face(:,1:3));
  title(['iter=' num2str(i)]);
  axis equal;
  set(gca,'xlim',[p0(1),p1(1)],'ylim',[p0(2),p1(2)],'zlim',[p0(3),p1(3)])
end


%=========================================================
% apply Low-pass smoothing
%=========================================================

figure;
subplot(rownum,colnum,1);

plotmesh(node,face(:,1:3));
if(~isoctavemesh)
        title({'Low-pass Smoothing Test','no smoothing'});
else
        title('Low-pass - no smoothing');
end
axis equal;
set(gca,'xlim',[p0(1),p1(1)],'ylim',[p0(2),p1(2)],'zlim',[p0(3),p1(3)])

conn=meshconn(face(:,1:3),size(node,1));

n1=node;
for i=1:rownum*colnum-1
  n1=smoothsurf(n1,[],conn,1,0.5,'lowpass',0.5);
  subplot(rownum,colnum,i+1);
  plotmesh(n1,face(:,1:3));
  title(['iter=' num2str(i)]);
  axis equal;
  set(gca,'xlim',[p0(1),p1(1)],'ylim',[p0(2),p1(2)],'zlim',[p0(3),p1(3)])
end

