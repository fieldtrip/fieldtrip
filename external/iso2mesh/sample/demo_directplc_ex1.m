%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% angular resolution
dt=pi/40;		
t=0:dt:2*pi-dt;

%% parameters of a side-cut fiber
h=100;
r1=20;
r2=25;
a=-1;	 
b=0;
c=1;
d=-h;

%% key nodes of a side-cut fiber

n1=[r1*sin(t(:)) r1*cos(t(:)) zeros(size(t(:)))];
n2=[r2*sin(t(:)) r2*cos(t(:)) zeros(size(t(:)))];

n3=[r1*sin(t(:)) r1*cos(t(:)) -d-(a*r1*sin(t(:))+b*r1*cos(t(:)))/c];
n4=[r2*sin(t(:)) r2*cos(t(:)) -d-(a*r2*sin(t(:))+b*r2*cos(t(:)))/c];

no=[n1;n2;n3;n4];

%% PLCs of the side-cut fiber

clear fc;
count=1;  
for i=1:length(t)-1
   % the last number in each cell is the fc id
   fc{count}={[i+length(t) i+3*length(t) i+3*length(t)+1 i+length(t)+1],1}; count=count+1;
   fc{count}={[i i+2*length(t) i+2*length(t)+1 i+1],2}; count=count+1;
end
i=length(t);
fc{count}={[i+length(t) i+3*length(t) 1+3*length(t) 1+length(t)],1}; count=count+1;
fc{count}={[i i+2*length(t) 1+2*length(t) 1],2}; count=count+1;

fc{count}={1:1+length(t)-1,3};count=count+1;  % bottom inner circle
fc{count}={[1+length(t):1+length(t)*2-1 nan fliplr(1:1+length(t)-1)],4};count=count+1; % button outter circle
fc{count}={1+length(t)*2:1+length(t)*3-1,5};count=count+1;  % top inner circle
fc{count}={[1+length(t)*3:1+length(t)*4-1 nan fliplr(1+length(t)*2:1+length(t)*3-1)],6}; % top outter circle

%% mesh generation of the cladding for the side-cut fiber
%[node,elem,face]=s2m(no,face,1,50);
[node,elem,face]=surf2mesh(no,fc,min(no),max(no),1,50,[0 0 1],[],0);

plotmesh(node,elem,'x>0 | y>0');
axis equal;
