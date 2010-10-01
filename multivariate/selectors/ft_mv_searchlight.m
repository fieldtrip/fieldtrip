classdef ft_mv_searchlight < ft_mv_selector
%FT_MV_SEARCHLIGHT searchlight analysis approach to feature selection
%
% The searchlight uses a sphere with a certain radius and a certain
% stepsize in order to scan a volume and classify based on elements inside
% the sphere.
%
% NOTE:
% A searchlight only returns the feature subset. Subsequently a predictor
% must be applied to it.
%
% Copyright (c) 2010, Marcel van Gerven

  properties
  
    indims                % dimensions of the input data
    radius                % radius of the hypersphere in terms of array elements
    step                  % stepsize in terms of array elements    
    mask                  % optional logical mask of size indims (input features are only those in the mask)
    neighbours            % a sparse adjacency matrix specifying the neighbourhood structure for irregular data (don't use in conjunction with mask)
    
    center = false;       % only selects center features if true
    nspheres = 1;         % select the union of the best nspheres spheres
    compact = true;       % save validator if compact = false
    
    centers               % centers of each sphere
    spheres               % the features belonging to each sphere
    original              % in original space
    value                 % evaluation metric
    pvalue                % significance value
    
    vld                   % saved validators if compact = false
    
  end
  
  methods
    
    function obj = ft_mv_searchlight(varargin)
      
      obj = obj@ft_mv_selector(varargin{:});
            
    end
    
    function obj = train(obj,X,Y)

      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % estimate spheres
      
      if isempty(obj.spheres)
        
        [obj.centers,obj.spheres,obj.original] = obj.estimate_spheres();
      
      end
      
      % test each sphere
      nsp = length(obj.spheres);
      
      obj.value = zeros(nsp,1);
      obj.pvalue = zeros(nsp,1);
      obj.vld = cell(nsp,1);
      for c=1:nsp

        vld = obj.validator;
        
        if iscell(X)
          XX = cell(size(X));
          for cc=1:length(X)
            XX{cc} = X{c}(:,obj.spheres{c});
          end
          vld = vld.train(XX,Y);
        else
          vld = vld.train(X(:,obj.spheres{c}),Y);
        end
        
        % report performance
        obj.value(c) = vld.performance();
        
        % report significance
        obj.pvalue(c) = vld.significance();
        
        if obj.verbose
          fprintf('performance for sphere %d of %d: %g (p-value: %g)\n',c,length(obj.spheres),obj.value(c),obj.pvalue(c));
        end

        % save validator
        if ~obj.compact
          obj.vld{c} = vld;
        end
        
      end
      
      [a,b] = sort(obj.value,'descend');
      if obj.center
        obj.subset = unique(subv2ind(obj.indims,obj.centers(b(1:obj.nspheres),:)));
      else
        obj.subset = unique(cell2mat(obj.spheres(b(1:obj.nspheres))));
      end
      
      
    end
    
    function [centers,spheres,original] = estimate_spheres(obj)
      
      % dimensions are retrieved from mask
      if ~isempty(obj.mask)
        obj.indims = size(obj.mask);
        midx = find(obj.mask(:));
      end
      
      % set default radius
      if isempty(obj.radius)
        obj.radius =  max(1,floor(min(obj.indims)./4));
      end
      rad = obj.radius;
      
      % set default step size
      if isempty(obj.step)
        obj.step = max(1,floor(min(obj.indims)./4));
      end
      stepsize = obj.step;
      
      % identify centers
      dd = cell(1,numel(obj.indims));
      for c=1:numel(obj.indims)
        dd{c} = 1:stepsize:obj.indims(c);
      end
      
      centers = cartprod(dd{:});
      
      % identify centers which are inside the mask
      if ~isempty(obj.mask)
     
        cidx = subv2ind(obj.indims,centers);
        centers = centers(logical(obj.mask(cidx)),:);
        
      end
      
      % identify subsets which fall in each sphere
      
      if obj.verbose
        fprintf('estimating %d spheres with radius %g with %g steps for a volume of size [%s]\n',...
          size(centers,1),rad,stepsize,num2str(obj.indims));
      end
      
      spheres = cell(size(centers,1),1);
      n=0;
      
      if ~isempty(obj.neighbours)
        
        if obj.verbose
          fprintf('building neighbourhood\n');
        end
        
        % centers as variable indices
        centers = subv2ind(obj.indims,centers);
        
        nidx = cell(size(obj.neighbours,1),1);
        for c=1:length(nidx)
          nidx{c} = find(obj.neighbours(c,:));
        end
            
        for c=1:length(p.spheres)
          
          if obj.verbose
            fprintf('estimating sphere %d of %d\n',c,size(centers,1));
          end
            
          % explicit neighbourhood structure for irregular data
          % radius becomes path length in the adjacency matrix
          
          spheres{c} = centers(c);
          for j=1:obj.radius
            spheres{c} = unique([spheres{c} cell2mat(nidx(spheres{c})')]);
          end
          
        end
        
      else
        
        
        for c=1:length(spheres)
          
          if obj.verbose
            fprintf('estimating sphere %d of %d\n',c,size(centers,1));
          end
          
          
          % generate potential points (the cube)
          v = cell(1,numel(obj.indims));
          for j=1:numel(obj.indims)
            v{j} = floor(centers(c,j)-rad):ceil(centers(c,j)+rad);
          end
          
          cube = cartprod(v{:});
          
          % reduce to elements inside the cube
          cube = cube(all(cube>0,2) & all(bsxfun(@minus,cube,obj.indims)<=0,2),:);
          
          % reduce to index elements inside the sphere
          spheres{c} = subv2ind(obj.indims,cube(sqrt(sum(bsxfun(@minus,cube,centers(c,:)).^2,2)) <= rad,:));
          
        end
        
        % centers as variable indices
        centers = subv2ind(obj.indims,centers);
        
      end
      
      % make spheres consistent with mask
      if ~isempty(obj.mask)

        tmp = spheres;
        for c=1:length(spheres)
          
          % reduce to elements inside the mask
          tmp{c} = spheres{c}(logical(obj.mask(spheres{c})));
          
          % transform to mask indices
          [a,spheres{c}] = ismember(tmp{c},midx);
          
        end
        
        % original sphere indices
        original = tmp;
        
      end
      
      if obj.verbose
        fprintf('average sphere volume: %g\n',sum(cellfun(@(x)(numel(x)),spheres))/length(spheres));
      end
      
    end
    
    function [m,desc] = model(obj)
      
      s = obj.spheres;
      v = obj.value;
      
      bmask = ~isempty(obj.mask);
      msk = find(obj.mask);
      
      m = zeros(obj.indims);
      n = zeros(obj.indims);
      for c=1:length(s)
                
        if bmask
          sidx = msk(s{c});
        else
          sidx = s{c};
        end
        
        m(sidx) = m(sidx) + v(c);
        n(sidx) = n(sidx) + 1;
        
      end
      
      m(n(:)~=0) = m(n(:)~=0) ./ n(n(:)~=0);
      
      m = {m};
      
      desc = {'average searchlight performance per feature'};
      
    end
    
  end
  
end