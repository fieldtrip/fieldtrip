function p=smoothsurf(node,mask,conn,iter,useralpha,usermethod,userbeta)
%
% p=smoothsurf(node,mask,conn,iter,useralpha,usermethod,userbeta)
%
% smoothing a surface mesh
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/21
%
% input:
%    node:  node coordinates of a surface mesh
%    mask:  flag whether a node is movable: 0 movable, 1 non-movable
%           if mask=[], it assumes all nodes are movable
%    conn:  input, a cell structure of length size(node), conn{n}
%           contains a list of all neighboring node ID for node n,
%           this can be computed from meshconn function
%    iter:  smoothing iteration number
%    useralpha: scaler, smoothing parameter, v(k+1)=(1-alpha)*v(k)+alpha*mean(neighbors)
%    usermethod: smoothing method, including 'laplacian','laplacianhc' and 'lowpass'
%    userbeta: scaler, smoothing parameter, for 'laplacianhc'
%
% output:
%    p: output, the smoothed node coordinates
%
% recommendations
%    Based on [Bade2006], 'Lowpass' method outperforms 'Laplacian-HC' in volume
%    preserving and both are significantly better than the standard Laplacian method
%
%    [Bade2006]  R. Bade, H. Haase, B. Preim, "Comparison of Fundamental Mesh 
%                Smoothing Algorithms for Medical Surface Models," 
%                Simulation and Visualization, pp. 289-304, 2006. 
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

p=node;
if(isempty(mask))
    nn=size(node,1);
    idx=1:nn;
else
    idx=find(mask==0)';
    nn=length(idx);
end
alpha=0.5;
method='laplacian';
beta=0.5;
if(nargin>4)
    alpha=useralpha;
    if(nargin>5)
        method=usermethod;
        if(nargin>6)
            beta=userbeta;
        end
    end
end
ibeta=1-beta;
ialpha=1-alpha;

for i=1:nn
    if(length(conn{idx(i)})==0)
        idx(i)=0;
    end
end
idx=idx(idx>0);
nn=length(idx);

if(strcmp(method,'laplacian'))
    for j=1:iter
        for i=1:nn
            p(idx(i),:)=ialpha*p(idx(i),:)+alpha*mean(node(conn{idx(i)},:)); 
        end
        node=p;
    end
elseif(strcmp(method,'laplacianhc'))
    for j=1:iter
        q=p;
        for i=1:nn
            p(idx(i),:)=mean(q(conn{idx(i)},:));
        end
        b=p-(alpha*node+ialpha*q);
        for i=1:nn
            p(idx(i),:)=p(idx(i),:)-(beta*b(i,:)+ibeta*mean(b(conn{idx(i)},:))); 
        end
    end
elseif(strcmp(method,'lowpass'))
    beta=-1.02*alpha;
    ibeta=1-beta;
    for j=1:iter
        for i=1:nn
            p(idx(i),:)=ialpha*node(idx(i),:)+alpha*mean(node(conn{idx(i)},:)); 
        end
        node=p;
        for i=1:nn
            p(idx(i),:)=ibeta *node(idx(i),:)+beta *mean(node(conn{idx(i)},:)); 
        end
        node=p;
    end
end
