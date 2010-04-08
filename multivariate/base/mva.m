classdef mva
%MVA multivariate analysis
%   
%   A multivariate analysisis just a sequence of multivariate methods 
%   {method1 method2 method3 ...} that are
%   called in this order and where the output of the previous method acts
%   as input to the next method.
%
%   An mvmethod always has a train and test function for the
%   two different modes of operation. An mva is
%   just a sequence of methods {method1 method2 method3 ...} that are
%   called in this order and where the output of the previous method acts
%   as input to the next method. Note that method derives from handle
%   so any changes are permanent. Methods:
%
%   PREPROCESSOR: produces transformed data
%   FEATURESELECTOR: selects and/or extracts features and produces data
%   CLASSIFIER: produces posterior probabilities
%   REGRESSOR: produces real-valued predictions
%   OPTIMIZER: wrapper class that optimizes parameters of another method
%   
%   the output of a classifier can be used as input to:
%
%   evaluate: takes posteriors and labels and produces an evaluation measure
%
%   parameters of a method should be set by means of a constructor call
%   method(obj,varargin)
%   
%   data is given by an N x M matrix with N examples/times and M features
%   design is an N x 1 matrix with the N class labels 1,2,3,...
%
%   multiple datasets are entered using cell arrays of data/design
%
%   Options:
%   'verbose' : verbose output [false]
%
%   Copyright (c) 2008, Marcel van Gerven

    properties

        mvmethods; % the methods that specify the classification procedure
        nmethods;   % the number of methods
        
        verbose = false;
    end
    
    methods
      
       function obj = mva(mvmethods,varargin)
       % constructor expects mva methods

        if ~nargin
            error('methods not specified');
        end
        
        if isa(mvmethods,'mva')
          obj = mvmethods;
          return;
        end
        
        % cast to cell array if only one method is specified
        if ~iscell(mvmethods), mvmethods = { mvmethods }; end

        % check if this is a valid classification procedure:
        if ~(all(cellfun(@(x)(isa(x,'mvmethod') || iscell(x)),mvmethods)))
          error('invalid multivariate analysis; check contents of mva');
        end
        
        obj.mvmethods = mvmethods;
        obj.nmethods = length(mvmethods);
       
       end
              
       function obj = train(obj,data,design)
         % train just calls the methods' train functions in order to produce a posterior
         
         if nargin<3
           design = [];
         end
         
         for c=1:obj.nmethods
           
           if iscell(obj.mvmethods{c})
            % deal with ensemble methods (i.e., nested cell arrays)
             
             if iscell(data) && ~obj.mvmethods{c}{1}.istransfer()
               % apply each submethod to each of the datasets
               
               if length(data) == length(obj.mvmethods{c})
                 for d=1:length(data)
                   obj.mvmethods{c}{d} = cellfun(@(x)(x.train(data{d},design{d})),obj.mvmethods{c}{d},'UniformOutput',false);
                   if c<obj.nmethods
                     data{d} = cellfun(@(x)(x.test(data{d})),obj.mvmethods{c}{d},'UniformOutput',false);
                   end
                 end
               else
                 error('cannot handle multiple datasets');
               end
               
             else
               % apply each submethod to the replicated dataset
               
               obj.mvmethods{c} = cellfun(@(x)(x.train(data,design)),obj.mvmethods{c},'UniformOutput',false);
               if c<obj.nmethods
                 data = cellfun(@(x)(x.test(data)),obj.mvmethods{c},'UniformOutput',false);
                 design = repmat({design},size(data));
               end
             end
             
           else
             
             if iscell(data) && ~obj.mvmethods{c}.istransfer()
               % if the method is not a transfer learner then we apply
               % the method to each dataset separately and convert
               % the method to a cell array
               
               m = cell(1,length(data));
               for d=1:length(data)
                 
                 m{d} = obj.mvmethods{c}.train(data{d},design{d});
                 if c<obj.nmethods
                   data{d} = m{d}.test(data{d});
                 end
               end
               
               obj.mvmethods{c} = m;
               
             else
               
               obj.mvmethods{c} = obj.mvmethods{c}.train(data,design);
               if c<obj.nmethods
                 data = obj.mvmethods{c}.test(data);
               end
             end
             
           end
         end
         
       end
       
       function data = test(obj,data)
         % test just calls the methods' test functions in order to produce a posterior
         
         if isempty(data)
           data = [];
           return;
         end
         
         for c=1:obj.nmethods
           
           if iscell(obj.mvmethods{c})
             % deal with ensemble methods
             
             if iscell(data) && ~obj.mvmethods{c}{1}.istransfer()
               
               if length(data) == length(obj.mvmethods{c})
                 for d=1:length(data)
                   data{d} = obj.mvmethods{c}{d}.test(data{d});
                 end
               else
                 error('cannot handle multiple datasets');
               end
               
             else
               
               data = cellfun(@(x)(x.test(data)),obj.mvmethods{c},'UniformOutput',false);
             end
             
           else
             
             data = obj.mvmethods{c}.test(data);
             
           end
           
         end
         
       end
       
       function data = untest(obj,data)
         % test just calls the methods' untest functions in order to
         % produce an output; note that method calling is in reversed
         % order
         
         if isempty(data)
           data = [];
           return;
         end
         
         for c=obj.nmethods:-1:1
           
           if iscell(obj.mvmethods{c})
             % deal with ensemble methods
             
             if iscell(data) && ~obj.mvmethods{c}{1}.istransfer()
               
               if length(data) == length(obj.mvmethods{c})
                 for d=1:length(data)
                   data{d} = obj.mvmethods{c}{d}.untest(data{d});
                 end
               else
                 error('cannot handle multiple datasets');
               end
               
             else
               
               data = cellfun(@(x)(x.untest(data)),obj.mvmethods{c},'UniformOutput',false);
             end
             
           else
             
             data = obj.mvmethods{c}.untest(data);
             
           end
           
         end
         
       end
       
       function p = predict(obj,data)
         % calls test functions and converts posterior into classifications
         
         for c=1:(length(obj.mvmethods)-1)
           data = obj.mvmethods{c}.test(data);
         end
         
         p = obj.mvmethods{end}.predict(data);
       end
       
       function s = name(obj)
         % returns multivariate analysis as a string
         
         s = subname(obj.mvmethods);
         
         function ss = subname(m)
            
           if iscell(m)
             ss = ['{ ' cell2mat(cellfun(@(x)([subname(x) ' ']),m,'UniformOutput',false)) '}'];
           else
             ss = class(m);
           end
           
         end
         
       end
       
       function [m,desc] = getmodel(obj)
         % return the model implied by this classification procedure
         
         % get the final model
         if iscell(obj.mvmethods{end})
           
           ntasks = length(obj.mvmethods{end});
           
           mm = cell(1,ntasks);
           for c=1:ntasks
             
             [mm{c},desc] = obj.mvmethods{end}{c}.getmodel();
             
           end
           
           m = cell(size(mm{1},1),ntasks);
           for c=1:ntasks
             for i=1:size(mm{1},1)
               m{i,c} = mm{c}{i};
             end
           end
           
         else
           [m,desc] = obj.mvmethods{end}.getmodel();
         end
         
       end
       
    end
    
end
