classdef ft_mv_searchlight < ft_mv_selector
%FT_MV_SEARCHLIGHT searchlight analysis approach to feature selection
%
% The searchlight uses a sphere with a certain radius and a certain
% stepsize in order to scan a volume and classify based on elements inside
% the sphere. Subsequently it selects the best n spheres and returns the
% optimal feature subset. The number nspheres is computed using an inner
% cross-validation with the specified validator.
%
% EXAMPLE:
% [a,b,c] = ft_mv_test('mva',{ft_mv_searchlight('step',5,'radius',5,'mask',masks{1},'validator',ft_mv_crossvalidator('sigtest','binomial','metric','accuracy','mva',ft_mv_svm,'nfolds',5),'verbose',true) ft_mv_svm},'nfolds',5)
%
% NOTE:
% A searchlight only returns the feature subset. Subsequently a predictor
% must be applied to it.
%
% Copyright (c) 2010, Marcel van Gerven

  properties
  
    indims                % dimensions of the input data
    radius = 3            % radius of the hypersphere in terms of array elements (diameter will be 1 + 2 * radius).
    step = 1              % stepsize in terms of array elements    
    mask                  % optional logical mask of size indims (input features are only those in the mask)
    neighbours            % a sparse adjacency matrix specifying the neighbourhood structure for irregular data (don't use in conjunction with mask)
    
    center = false;       % only selects the feature in the centre of the sphere if true
    compact = true;       % save validator if compact = false

    nspheres = [];        % select the union of the best nspheres spheres; 
                          % if empty then use the validator to compute the
                          % best n; if array then test all values and
                          % select best n. Set to 0 if you are only
                          % interested in getting the performance estimates
                          % per sphere

    performance = [];     % performance as function of number of selected spheres (or sphere centers)
    
    centers               % centers of each sphere
    spheres               % the features belonging to each sphere
    original              % in original space
    value                 % evaluation metric
    pvalue                % significance value (uncorrected for the number of spheres)
    
    vld                   % saved validators if compact = false
    
    exclude = true;       % if false then only the center within a sphere is required to be part of the optional mask
    
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
      
      % dimensions are retrieved from optional mask
      if ~isempty(obj.mask)
        obj.indims = size(obj.mask);
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
      
      if obj.nspheres~=0
        
        % sort spheres according to performance
        [a,b] = sort(obj.value,'descend');
        
        % use validator to compute best spheres
        if isempty(obj.nspheres)
          obj.nspheres = 1:length(obj.spheres);
        end
        
        if length(obj.nspheres) > 1
          
          if obj.verbose
            fprintf('determining optimal number of spheres\n');
          end
          
          obj.performance = zeros(1,length(obj.nspheres));
          for nn=1:length(obj.nspheres)
            
            n = obj.nspheres(nn);
            
            if obj.center
              obj.subset = unique(subv2ind(obj.indims,obj.centers(b(1:n),:)));
            else
              obj.subset = unique(cell2mat(obj.spheres(b(1:n))));
            end
            
            if iscell(X)
              XX = cell(size(X));
              for cc=1:length(X)
                XX{cc} = X{c}(:,obj.subset);
              end
              v = obj.validator.train(XX,Y);
            else
              v = obj.validator.train(X(:,obj.subset),Y);
            end
            
            obj.performance(n) = v.performance();
            
            if obj.verbose
              fprintf('%d spheres : %f performance\n',n,obj.performance(nn));
            end
            
          end
          
          [tmp,optn] = max(obj.performance);
          
          optn = obj.nspheres(optn);
          
        else
          optn = obj.nspheres;
        end
        
        % get final subset
        if obj.center
          obj.subset = unique(subv2ind(obj.indims,obj.centers(b(1:optn),:)));
        else
          obj.subset = unique(cell2mat(obj.spheres(b(1:optn))));
        end
        
        if obj.verbose
          fprintf('selected %d spheres\n',optn);
        end
        
      end
      
    end
    
    function [centers,spheres,original] = estimate_spheres(obj)
      
      assert(~isempty(obj.indims));
      
      % dimensions are retrieved from mask
      if ~isempty(obj.mask)
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
      if obj.exclude && ~isempty(obj.mask)

        tmp = spheres;
        for c=1:length(spheres)
          
          % reduce to elements inside the mask
          tmp{c} = spheres{c}(logical(obj.mask(spheres{c})));
          
          % transform to mask indices
          [a,spheres{c}] = ismember(tmp{c},midx);
          
        end
        
        % original sphere indices
        original = tmp;
        
      else
        
        original = spheres;
       
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