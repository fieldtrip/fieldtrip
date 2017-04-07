function hm=plotmesh(node,varargin)
%
% hm=plotmesh(node,face,elem,opt)
%
% plot surface and volumetric meshes
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%      node: a node coordinate list, 3 columns for x/y/z; if node has a 
%            4th column, it will be used to set the color at each node.
%      face: a triangular surface face list; if face has a 4th column,
%            it will be used to separate the surface into 
%            sub-surfaces and display them in different colors;
%            face can be a cell array, each element of the array represents
%            a polyhedral facet of the mesh, if an element is an array with
%            two array subelements, the first one is the node index, the
%            second one is a scalar as the group id of the facet.
%      elem: a tetrahedral element list; if elem has a 5th column,
%            it will be used to separate the mesh into 
%            sub-domains and display them in different colors.
%      opt:  additional options for the plotting
%
%            for simple point plotting, opt can be markers
%            or color options, such as 'r.', or opt can be 
%            a logic statement to select a subset of the mesh,
%            such as 'x>0 & y+z<1', or an equation defining
%            a plane at which a mesh cross-section is plotted, for
%            example 'y=2*x'; opt can have more than one
%            items to combine these options, for example: 
%            plotmesh(...,'x>0','r.'); the range selector must
%            appear before the color/marker specifier
%
% in the event where all of the above inputs have extra settings related to 
% the color of the plot, the priorities are given in the following order:
%
%          opt > node(:,4) > elem(:,5) > face(:,4)
%
% output:
%   hm: handle or handles (vector) to the plotted surfaces
%
% example:
%
%   h=plotmesh(node,'r.');
%   h=plotmesh(node,'x<20','r.');
%   h=plotmesh(node,face);
%   h=plotmesh(node,face,'y>10');
%   h=plotmesh(node,face,'facecolor','r');
%   h=plotmesh(node,elem,'x<20');
%   h=plotmesh(node,elem,'x<20 & y>0');
%   h=plotmesh(node,face,elem);
%   h=plotmesh(node,face,elem,'linestyle','--');
%   h=plotmesh(node,elem,'z=20');
% 
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

selector=[];
opt=[];
face=[];
elem=[];

if(nargin>1)
   hasopt=0;
   for i=1:length(varargin)
   	if(ischar(varargin{i}))
		if(~isempty(regexp(varargin{i},'[x-zX-Z]')) && ~isempty(regexp(varargin{i},'[><=&|]')))
			selector=varargin{i};
			if(nargin>=i+1) opt=varargin(i+1:end); end
		else
			opt=varargin(i:end);
		end
		if(i==1)
			face=[];elem=[];
		elseif(i==2)
			if(iscell(varargin{1}) | size(varargin{1},2)<4)
				face=varargin{1}; elem=[];
			elseif(size(varargin{1},2)==4)
                faceid=unique(varargin{1}(:,4));
                if(length(faceid)==1)
                    face=varargin{1}; elem=[];
                elseif(any(hist(varargin{1}(:,4),unique(varargin{1}(:,4)))>50))
                    face=varargin{1}; elem=[];
                else
                    elem=varargin{1}; face=[];
                end
			else
				elem=varargin{1}; face=[];
			end
		elseif(i==3)
			face=varargin{1};
			elem=varargin{2};
		end
		hasopt=1;
		break;
	end
   end
   if(hasopt==0)
   	if(length(varargin)>=2)
		face=varargin{1};
		elem=varargin{2};
		if(length(varargin)>2) opt=varargin(3:end); end
	elseif(iscell(varargin{1}) | size(varargin{1},2)<4)
		face=varargin{1}; elem=[];
	elseif(size(varargin{1},2)==4)
	    faceid=unique(varargin{1}(:,4));
            if(length(faceid)==1)
	        face=varargin{1}; elem=[];
	    elseif(any(hist(varargin{1}(:,4),unique(varargin{1}(:,4)))>50))
                face=varargin{1}; elem=[];
	    else
                elem=varargin{1}; face=[];
	    end
	else
		elem=varargin{1}; face=[];
	end
   end
end

holdstate=ishold;
if(~holdstate)
    cla;
end
if(size(node,2)==4 && size(elem,2)==5)
    warning(['You have specified the node colors by both the 4th ' ...
            'and 5th columns of node and face inputs, respectively. ' ...
            'The node input takes priority']);
end
if(isempty(face) && isempty(elem))
   if(isempty(selector))
        if(isempty(opt))
   		h=plot3(node(:,1),node(:,2),node(:,3),'o');
	else
   		h=plot3(node(:,1),node(:,2),node(:,3),opt{:});
	end
   else
	x=node(:,1);
	y=node(:,2);
	z=node(:,3);
	idx=eval(['find(' selector ')']);
    if(~isempty(idx))
	    if(isempty(opt))
		h=plot3(node(idx,1),node(idx,2),node(idx,3),'o');
	    else
		h=plot3(node(idx,1),node(idx,2),node(idx,3),opt{:});
        end
    else
        warning('nothing to plot');
	end
   end
end

if(~isempty(face))
   hold on;
   if(isempty(selector))
        if(isempty(opt))
   		h=plotsurf(node,face);
	else
   		h=plotsurf(node,face,opt{:});
	end
   else
    if(iscell(face))
       cent=meshcentroid(node,face);
    else
       cent=meshcentroid(node,face(:,1:3));
    end
	x=cent(:,1);
    y=cent(:,2);
	z=cent(:,3);
    idx=eval(['find(' selector ')']);
    if(~isempty(idx))
        if(iscell(face))
            h=plotsurf(node,face(idx),opt{:});
        else
    		h=plotsurf(node,face(idx,:),opt{:});
        end
    else
        warning('no surface to plot');
	end
   end
end

if(~isempty(elem))
   hold on;
   if(isempty(selector))
        if(isempty(opt))
   		h=plottetra(node,elem);
	else
   		h=plottetra(node,elem,opt{:});
	end
   else
   cent=meshcentroid(node,elem(:,1:4));
   x=cent(:,1);
   y=cent(:,2);
   z=cent(:,3);
   if(regexp(selector,'='))
      if(size(node,2)==4)
          [cutpos,cutvalue,facedata]=qmeshcut(elem,node(:,1:3),node(:,4),selector);  
      elseif(size(node,2)==3)
          [cutpos,cutvalue,facedata]=qmeshcut(elem,node,node(:,3),selector);
      else
          error('plotmesh can only plot 3D tetrahedral meshes');
      end
      h=patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,'facecolor','interp',opt{:});
   else
      idx=eval(['find(' selector ')']);
      if(~isempty(idx))
	    if(isempty(opt))
		h=plottetra(node,elem(idx,:));
	    else
		h=plottetra(node,elem(idx,:),opt{:});
        end
      else
        warning('no tetrahedral element to plot');
	end
     end
   end
end

if(exist('h','var') & ~holdstate)
  hold off;
end
if(exist('h','var'))
  if(any(get(gca,'dataaspectratio')>1e8))
     view(3);
  end
  axis equal;
end
if(exist('h','var') & nargout>=1)
  hm=h;
end
