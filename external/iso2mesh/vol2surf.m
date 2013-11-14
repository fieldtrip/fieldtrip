function [no,el,regions,holes]=vol2surf(img,ix,iy,iz,opt,dofix,method,isovalues)
%
% [no,el,regions,holes]=vol2surf(img,ix,iy,iz,opt,dofix,method,isovalues)
%
% converting a 3D volumetric image to surfaces
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 img: a volumetric binary image; if img is empty, vol2surf will
%	      return user defined surfaces via opt.surf if it exists
%	 ix,iy,iz: subvolume selection indices in x,y,z directions
%	 opt: function parameters
%	   if method is 'cgalsurf' or 'cgalpoly':
%	     opt=a float number>1: max radius of the Delaunay sphere(element size) 
%	     opt.radbound: same as above, max radius of the Delaunay sphere
%	     opt.distbound: maximum deviation from the specified isosurfaces
%	     opt(1,2,...).radbound: same as above, for each levelset
%	   if method is 'simplify':
%	     opt=a float number<1: compression rate for surf. simplification
%	     opt.keeyratio=a float less than 1: same as above, same for all surf.
%	     opt(1,2,..).keeyratio: setting compression rate for each levelset
%	   opt(1,2,..).maxsurf: 1 - only use the largest disjointed surface
%				0 - use all surfaces for that levelset
%          opt(1,2,..).side: - 'upper': threshold at upper interface
%                              'lower': threshold at lower interface
%	   opt(1,2,..).maxnode: - the maximum number of surface node per levelset
%	   opt(1,2,..).holes: user specified holes interior pt list
%	   opt(1,2,..).regions: user specified regions interior pt list
%	   opt(1,2,..).surf.{node,elem}: add additional surfaces
%	   opt(1,2,..).{A,B}: linear transformation for each surface
%	   opt.autoregion: if set to 1, vol2surf will try to determine 
%              the interior points for each closed surface automatically
%	 dofix: 1: perform mesh validation&repair, 0: skip repairing
%	 method: - if method is 'simplify', iso2mesh will first call
%		   binsurface to generate a voxel-based surface mesh and then
%		   use meshresample/meshcheckrepair to create a coarser mesh;
%		 - if method is 'cgalsurf', iso2mesh will call the surface
%		   extraction program from CGAL to make surface mesh
%		 - if method is not specified, 'cgalsurf' is assumed by default
%	 isovalues: a list of isovalues where the levelset is defined
%
% output: 
%	 no:  list of nodes on the resulting suface mesh, 3 columns for x,y,z
%	 el:  list of trianglular elements on the surface, [n1,n2,n3,region_id]
%	 regions: list of interior points for all sub-region, [x,y,z]
%	 holes:   list of interior points for all holes, [x,y,z]
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fprintf(1,'extracting surfaces from a volume ...\n');

el=[];
no=[];

if(isstruct(opt) & isfield(opt,'holes')) 
    holes=opt.holes;
else
    holes=[];
end
if(isstruct(opt) & isfield(opt,'regions')) 
    regions=opt.regions;
else
    regions=[];
end
maxlevel=0;

if(~isempty(img))

    img=img(ix,iy,iz);
    dim=size(img);
    newdim=dim+[2 2 2];
    newimg=zeros(newdim);
    newimg(2:end-1,2:end-1,2:end-1)=img;

    if(nargin<8)
        maxlevel=max(newimg(:));
        isovalues=1:maxlevel;
    else
        isovalues=unique(sort(isovalues));
        maxlevel=length(isovalues);
    end

    for i=1:maxlevel
      if(i<maxlevel)
          levelmask=int8(newimg>=isovalues(i) & newimg<isovalues(i+1));
      else
          levelmask=int8(newimg>=isovalues(i));
      end
      [levelno,levelel]=binsurface(levelmask);
      if(~isempty(levelel))
          if(isstruct(opt) & isfield(opt,'autoregion'))
              if(opt.autoregion)
                  seeds=surfseeds(levelno,levelel);
              else
                  seeds=surfinterior(levelno,levelel);
              end
          else
              seeds=surfinterior(levelno,levelel);
          end
          if(~isempty(seeds))
              disp([sprintf('region %d centroid :',i) sprintf('\t%f %f %f\n', seeds')]);
              if(~isempty(regions))
                  regions(end+1:end+size(seeds,1),:)=seeds;
              else
                  regions=seeds;
              end
          end
      end
    end

    for i=1:maxlevel
        fprintf(1,'processing threshold level %d...\n',i);

        if(nargin>=7 & strcmp(method,'simplify'))

          [v0,f0]=binsurface(newimg>=isovalues(i)); % not sure if binsurface works for multi-value arrays
          % with binsurface, I think the following line is not needed anymore
          %  v0(:,[1 2])=v0(:,[2 1]); % isosurface(V,th) assumes x/y transposed
          if(dofix)  [v0,f0]=meshcheckrepair(v0,f0);  end  

          if(isstruct(opt) & length(opt)==maxlevel) keepratio=opt(i).keepratio;
          elseif (isstruct(opt) & length(opt)==1) keepratio=opt.keepratio;
          else keepratio=opt;  end;

          % first, resample the surface mesh with cgal
          fprintf(1,'resampling surface mesh for level %d...\n',i);
          [v0,f0]=meshresample(v0,f0,keepratio);

          % iso2mesh is not stable for meshing small islands,remove them (max 3x3x3 voxels)
          f0=removeisolatedsurf(v0,f0,3);

          if(dofix) [v0,f0]=meshcheckrepair(v0,f0); end

        elseif(nargin<7 | strcmp(method,'cgalsurf') | strcmp(method,'cgalpoly'))
          if(isstruct(opt) & length(opt)==maxlevel) radbound=opt(i).radbound;
          elseif (isstruct(opt) & length(opt)==1) radbound=opt.radbound;
          else radbound=opt;  end;

          maxsurfnode=40000;  % maximum node numbers for each level
          if(isstruct(opt) & length(opt)==maxlevel) 
              if(isfield(opt(i),'maxnode')) maxsurfnode=opt(i).maxnode; end
          elseif (isstruct(opt) & length(opt)==1 )
              if(isfield(opt(1),'maxnode')) 
                 maxsurfnode=opt.maxnode; 
              end
          end

	  distbound=radbound;
          if(isstruct(opt) & length(opt)==maxlevel)
              if(isfield(opt(i),'distbound')) distbound=opt(i).distbound; end
          elseif (isstruct(opt) & length(opt)==1 )
              if(isfield(opt(1),'distbound')) distbound=opt.distbound; end
          end
	  surfside='';
          if(isstruct(opt) & length(opt)==maxlevel)
              if(isfield(opt(i),'side')) surfside=opt(i).side; end
          elseif (isstruct(opt) & length(opt)==1)
              if(isfield(opt(1),'side')) surfside=opt(1).side; end
          end
	  if(~isempty(surfside))
	     newimg0=newimg;
	  end
          if(strcmp(surfside,'upper'))
	      newimg(find(newimg<=isovalues(i)-1e-9))=isovalues(i)-1e-9;
	  elseif(strcmp(surfside,'lower'))
	      newimg(find(newimg>=isovalues(i)+1e-9))=isovalues(i)+1e-9;
	  end
          perturb=1e-4*abs(max(isovalues));
          if(all(newimg>isovalues(i)-perturb)) perturb=-perturb;  end
          [v0,f0]=vol2restrictedtri(newimg,isovalues(i)-perturb,regions(i,:),...
                     sum(newdim.*newdim)*2,30,radbound,distbound,maxsurfnode);

	  if(~isempty(surfside))
            newimg=newimg0;
	    clear newimg0;
	  end
        else
            error('method can only be one of "cgalsurf", "cgalpoly" or "simplify".');
        end

        % if use defines maxsurf=1, take only the largest closed surface
        if(isstruct(opt))
            if( (isfield(opt,'maxsurf') && length(opt)==1 && opt.maxsurf==1) | ...
                    (length(opt)==maxlevel && isfield(opt(i),'maxsurf') && opt(i).maxsurf==1))
                    f0=maxsurf(finddisconnsurf(f0));
            end
        end

        % if a transformation matrix/offset vector supplied, apply them

        if(isstruct(opt) & length(opt)==maxlevel) 
          if(isfield(opt(i),'A') & isfield(opt(i),'B'))
            v0=(opt(i).A*v0'+repmat(opt(i).B(:),1,size(v0,1)))';
          end
        elseif (isstruct(opt) & length(opt)==1) 
          if(isfield(opt,'A') & isfield(opt,'B'))
            v0=(opt.A*v0'+repmat(opt.B(:),1,size(v0,1)))';
          end
        end

        % if user specified holelist and regionlist, append them
        if(isstruct(opt)  & length(opt)==maxlevel)
        if(isfield(opt(i),'hole'))
                holes=[holes;opt(i).hole]
        end
        if(isfield(opt(i),'region'))
                regions=[regions;opt(i).region]
        end
        end

        if(i==0)
        el=[f0 (i+1)*ones(size(f0,1),1)];
        no=v0;
        else
        el=[el;f0+length(no) (i+1)*ones(size(f0,1),1)];
        no=[no;v0];
        end
    end

    %some final fix and scaling
    no(:,1:3)=no(:,1:3)-1; % because we padded the image with a 1 voxel thick null layer in newimg

    no(:,1)=no(:,1)*(max(ix)-min(ix)+1)/dim(1)+(min(ix)-1);
    no(:,2)=no(:,2)*(max(iy)-min(iy)+1)/dim(2)+(min(iy)-1);
    no(:,3)=no(:,3)*(max(iz)-min(iz)+1)/dim(3)+(min(iz)-1);

end  % if not isempty(img)

if(isstruct(opt) & isfield(opt,'surf'))
   for i=1:length(opt.surf)
	opt.surf(i).elem(:,4)=maxlevel+i;
        el=[el;opt.surf(i).elem+length(no)];
        no=[no;opt.surf(i).node];
   end
end

fprintf(1,'surface mesh generation is complete\n');

