classdef searchlight < dml.method
% SEARCHLIGHT searchlight analysis.
%
%   DESCRIPTION
%   The searchlight uses a sphere with a certain radius and stepsize in order 
%   to scan a volume and classify based on elements inside the sphere. 
%   Subsequently one may either:
%   * Use m.model to return the performance estimates throughout the volume
%   * Select the best n spheres and retrieve this subset for the input data. 
%     The number nspheres is computed using a nested cross-validation with 
%     the specified validator.
%
%   EXAMPLE:
%
%   X = rand(10,1000); Y = [1 1 1 1 1 2 2 2 2 2]';
%
%   Perform searchlight and return classification accuracy computed with a 
%   default crossvalidator and mapped back to the standard volume. 
%
%   indims = [10 10 10];
%   m = dml.searchlight('step',2,'radius',2,'indims',indims,'verbose',true,'stats',{'accuracy'});
%   m = m.train(X,Y); r = m.model;
%
%   Alternatively, a multidimensional logical mask can be used to specify
%   the input dimensions and restrict to a subset of the original volume; the
%   original volume will still be used to estimate the searchlight
%   spheres:
%
%   mymask = true(10,10,10); mymask(:,:,1:5) = false;
%   m = dml.searchlight('step',2,'radius',2,'mask',mymask,'verbose',true,'stats',{'accuracy'});
%   m = m.train(X(:,find(mymask(:))),Y); r = m.model;
%
%   We may also specify a different, irregular, neighbourhood structure
%   in conjunction with the input dimensions or a particular mask:
%
%   mymask = true(10,10,10); mymask(:,:,1:5) = false;
%   nb = sparse(1000,1000); prm = randperm(1e6); nb(prm(1:1000)) = 1; nb = (nb + nb') ~= 0;
%   m = dml.searchlight('step',2,'radius',2,'neighbours',nb,'mask',mymask,'verbose',true,'stats',{'accuracy'});
%   m = m.train(X(:,find(mymask(:))),Y); r = m.model;
%
%   The output of this searchlight analysis can also be used for feature selection.
%   In that case, one needs to specify nspheres, containing the numbers of
%   spheres to test. The optimal subset of spheres will be determined using
%   the specified crossvalidator:
%
%   mymask = true(10,10,10);
%   m = dml.searchlight('nspheres',[1 2 3 4 5],'step',2,'radius',2,'mask',mymask,'verbose',true,'stats',{'accuracy'});
%   m = m.train(X,Y); r = m.model;
%
%   Instead of using a cross-validator, one may also use a permutation
%   object to determine sphere performance:
%
%   indims = [10 10 10];
%   p = dml.permutation('stat','accuracy','validator',dml.crossvalidator('mva',dml.svm),'nperm',10,'verbose',true);
%   m = dml.searchlight('step',3,'radius',2,'indims',indims,'verbose',true,'validator',p);
%   m = m.train(X,Y); r = m.model;
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
  
    radius = 3            % radius of the hypersphere in terms of array elements (diameter will be 1 + 2 * radius).
    step = 1              % stepsize in terms of array elements    
    
    mask                  % optional logical mask of size indims (input features are only those in the mask)
    neighbours            % a sparse adjacency matrix specifying the neighbourhood structure for irregular data (don't use in conjunction with mask)
    
    stats = {}            % the statistics to save in value (e.g. {'accuracy','binomial'}); empty stats field will just call validator.statistic
    
    center = false;       % only selects the feature in the centre of the sphere if true
    
    exclude = true;       % if false then only the center within a sphere is required to be part of the optional mask
     
    validator = dml.crossvalidator('mva',{dml.standardizer dml.svm}) % the validator used to determine the final used feature subset
    compact = true;       % save validators if compact = false
    
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
    
    value                 % evaluation metric per sphere
    
    subset                % the used feature subset
    
  end
  
  methods
    
    function obj = searchlight(varargin)
      
      obj = obj@dml.method(varargin{:});
            
    end
    
    function obj = train(obj,X,Y)
      
      if isempty(obj.mask) && isempty(obj.indims)
        error('please specify input dimensions or logical mask'); 
      end
      
      if isempty(obj.mask)
        obj.mask = true(obj.indims);
      end
      
      if ~obj.exclude && ((iscell(X) && size(X{1},2) ~= nnz(obj.mask)) || (~iscell(X) && size(X,2) ~= nnz(obj.mask)))
          error('number of features should match nonzero elements in mask');
      end
      obj.indims = size(obj.mask);
       
      if ~iscell(obj.stats), obj.stats = {obj.stats}; end
      
      % estimate spheres      
      if obj.restart || isempty(obj.spheres)
        [obj.centers,obj.spheres,obj.original] = obj.estimate_spheres();
      end
      
      % test each sphere
      nsp = length(obj.spheres);

      V =  obj.validator;
      if ~obj.compact, obj.validator = cell(nsp,1); end
      
      nstats = max(1,length(obj.stats));
      obj.value = zeros(nsp,nstats);
      for c=1:nsp

        if iscell(X)
          XX = cell(size(X));
          for cc=1:length(X)
            XX{cc} = X{cc}(:,obj.spheres{c});
          end
          vald = V.train(XX,Y);
        else
          vald = V.train(X(:,obj.spheres{c}),Y);
        end
        
        % report performance
        if isempty(obj.stats)
          obj.value(c) = vald.statistic;
        else
          for j=1:nstats
            obj.value(c,j) = vald.statistic(obj.stats{j});
          end
        end
        
        if obj.verbose
          fprintf('performance for sphere %d of %d: %g\n',c,length(obj.spheres),obj.value(c,1));
        end
        
        % save validator
        if ~obj.compact, obj.validator{c} = vald; end
        
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
              v = V.train(XX,Y);
            else
              v = V.train(X(:,obj.subset),Y);
            end
            
            obj.performance(n) = v.stats(obj.stats{1});
            
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
          fprintf('selected %d sphere(s)\n',optn);
        end
        
      end
      
    end
    
    function Y = test(obj,X)
            
      Y = X(:,obj.subset);
      
    end
    
    function [centers,spheres,original] = estimate_spheres(obj)
      
      assert(~isempty(obj.indims));
            
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
      if ~isempty(obj.mask) && ~all(obj.mask(:))
     
        cidx = subv2ind(obj.indims,centers);
        centers = centers(logical(obj.mask(cidx)),:);
        
      end
      
      % identify subsets which fall in each sphere
      
      if obj.verbose
        fprintf('estimating %d spheres with radius %g with %g steps for a volume of size [%s]\n',...
          size(centers,1),rad,stepsize,num2str(obj.indims));
      end
      
      spheres = cell(size(centers,1),1);
      
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
            
        for c=1:length(spheres)
          
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
      
      % make spheres fully consistent with mask
      if obj.exclude && ~all(obj.mask(:))

        midx = find(obj.mask(:));

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
    
    function m = model(obj)
      % return performance mapped back to original space
      % we return all quantities specified in the stats field
      
      for j=1:length(obj.stats)
      
        s = obj.spheres;
        v = obj.value(:,j);
        
        bmask = ~all(obj.mask(:));
        msk = find(obj.mask);
        
        mt = zeros(obj.indims);
        n = zeros(obj.indims);
        for c=1:length(s)
          
          if bmask
            sidx = msk(s{c});
          else
            sidx = s{c};
          end
          
          mt(sidx) = mt(sidx) + v(c);
          n(sidx) = n(sidx) + 1;
          
        end
        
        nidx = n(:)~=0;
        mt(nidx) = mt(nidx) ./ n(nidx);
        mt(n(:)==0) = nan;
        
        m.(obj.stats{j}) = mt;
        
      end
      
    end
    
  end
  
end