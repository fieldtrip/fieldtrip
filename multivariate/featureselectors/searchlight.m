classdef searchlight < featureselector
%SEARCHLIGHT searchlight analysis approach to feature selection
%
% The searchlight uses a sphere with a certain radius and a certain
% stepsize in order to scan a volume and classify based on elements inside
% the sphere.
%
% PARAMETERS:
% subset      % the best performing feature subset
% centers     % centers of the spheres
% spheres     % the feature subset per sphere
% value       % performance estimates per sphere
%
% NOTE:
% A searchlight only returns the feature subset. Subsequently a predictor
% must be applied to it.
%
% Copyright (c) 2010, Marcel van Gerven

  properties
  
    radius                % radius of the hypersphere in terms of array elements
    
    step                  % stepsize in terms of array elements    
   
    mask                  % optional mask
        
  end
  
  methods
    
    function obj = searchlight(varargin)
      
      obj = obj@featureselector(varargin{:});
      
      assert(~isempty(obj.procedure));
      
      if isempty(obj.validator)
        obj.validator = crossvalidator('verbose',true);
      end
      
      obj.validator.procedure = mva(obj.procedure);

    end
    
    function p = estimate(obj,X,Y)
      
      if isempty(obj.radius)
        obj.radius =  max(1,floor(min(obj.indims)./4));
      end
      rad = obj.radius;
      
      if isempty(obj.step)
        obj.step = max(1,floor(min(obj.indims)./4));
      end
      stepsize = obj.step;
      
      % identify centers
      dd = cell(1,numel(obj.indims));
      for c=1:numel(obj.indims)
        dd{c} = 1:stepsize:obj.indims(c);
      end
      
      p.centers = cartprod(dd{:});
      
      % identify centers which are inside the mask
      if ~isempty(obj.mask)
     
        cidx = subv2ind(obj.indims,p.centers);
        p.centers = p.centers(logical(obj.mask(cidx)),:);
        
      end
      
      % identify subsets which fall in each sphere
      
      if obj.verbose
        fprintf('estimating %d spheres with radius %g with %g steps for a volume of size [%s]\n',...
          size(p.centers,1),rad,stepsize,num2str(obj.indims));
      end
      
      if ~isempty(obj.mask)
        midx = find(obj.mask(:));
      end
      
      p.spheres = cell(size(p.centers,1),1); 
      n=0;
      for c=1:length(p.spheres)
        
        if obj.verbose
          fprintf('estimating sphere %d of %d\n',c,size(p.centers,1));
        end
        
        % generate potential points (the cube)
        v = cell(1,numel(obj.indims));
        for j=1:numel(obj.indims)
          v{j} = floor(p.centers(c,j)-rad):ceil(p.centers(c,j)+rad);
        end
        
        cube = cartprod(v{:});
        
        % reduce to elements inside the cube
        cube = cube(all(cube>0,2) & all(bsxfun(@minus,cube,obj.indims)<=0,2),:);
        
        % reduce to index elements inside the sphere
        p.spheres{c} = subv2ind(obj.indims,cube(sqrt(sum(bsxfun(@minus,cube,p.centers(c,:)).^2,2)) <= rad,:));
        
        if ~isempty(obj.mask)
          
          % reduce to elements inside the mask
          p.spheres{c} = p.spheres{c}(logical(obj.mask(p.spheres{c})));
          
          % transform to mask indices
          [a,p.spheres{c}] = ismember(p.spheres{c},midx);
          
        end
        
        n = n+length(p.spheres{c});
        
      end
      
      if obj.verbose
        fprintf('average sphere volume: %g\n',n/length(p.spheres));
      end
      
      % test each sphere
      
      p.value = zeros(length(p.spheres),1);
      for c=1:length(p.spheres)

        vld = obj.validator;
        
        vld = vld.validate(X(:,p.spheres{c}),Y);
           
        p.value(c) = vld.evaluate('metric',obj.metric);        

        if obj.verbose
          fprintf('performance for sphere %d of %d using metric %s: %g\n',c,length(p.spheres),obj.metric,p.value(c));
        end

      end
      
      [a,b] = max(p.value);
      p.subset = p.spheres{b};
      
    end
    
    function [m,desc] = getmodel(obj)
      
      s = obj.params.spheres;
      v = obj.params.value;
      
      if ~isempty(obj.mask)
        midx = find(obj.mask(:));
      end
      
      m = zeros(obj.indims);
      n = zeros(obj.indims);
      for c=1:length(s)
        
        if ~isempty(obj.mask)
          m(midx(s{c})) = m(midx(s{c})) + v(c);
          n(midx(s{c})) = n(midx(s{c})) + 1;
        else
          m(s{c}) = m(s{c}) + v(c);
          n(s{c}) = n(s{c}) + 1;
        end
      end
      
      m(n(:)~=0) = m(n(:)~=0) ./ n(n(:)~=0);
      
      m = {m};
      
      desc = {sprintf('average searchlight performance per feature using %s',obj.metric)};
      
    end
    
  end
  
end