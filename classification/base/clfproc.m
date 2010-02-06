classdef clfproc
%CLFPROC classification procedure class
%   
%   A classification procedure is just a sequence of classification methods 
%   {method1 method2 method3 ...} that are
%   called in this order and where the output of the previous method acts
%   as input to the next method.
%
%   A classification method always has a train and test function for the
%   two different modes of operation. A classification procedure is
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
%
%   $Log: clfproc.m,v $
%

    properties

        clfmethods; % the methods that specify the classification procedure
        nmethods;   % the number of methods
        
        verbose = false;
    end
    
    methods
      
       function obj = clfproc(clfmethods,varargin)
       % constructor expects classification methods

        if ~nargin
            error('methods not specified');
        end
        
        if isa(clfmethods,'clfproc')
          obj = clfmethods;
          return;
        end
        
        % cast to cell array if only one method is specified
        if ~iscell(clfmethods), clfmethods = { clfmethods }; end

        % check if this is a valid classification procedure:
        if ~(all(cellfun(@(x)(isa(x,'clfmethod') || iscell(x)),clfmethods)))
          error('invalid classification procedure; check contents of clfproc');
        end
        
        obj.clfmethods = clfmethods;
        obj.nmethods = length(clfmethods);
       
       end
              
       function obj = train(obj,data,design)
         % train just calls the methods' train functions in order to produce a posterior
         
         if nargin<3
           design = [];
         end
         
         % cast to datasets if necessary
         
         if iscell(data)
           for c=1:length(data)
             if ~isa(data{c},'dataset')
               data{c} = dataset(data{c});
             end
           end
         elseif ~isa(data,'dataset')
           data = dataset(data);
         end
         
         if iscell(design)
           for c=1:length(design)
             if ~isa(design{c},'dataset')
               design{c} = dataset(design{c});
             end
           end
         elseif ~isa(design,'dataset')
           design = dataset(design);
         end
         
         for c=1:obj.nmethods
           
           if iscell(obj.clfmethods{c})
             
             if iscell(data) && ~isa(obj.clfmethods{c}{1},'transfer_learner')
               
               if length(data) == length(obj.clfmethods{c})
                 for d=1:length(data)
                   obj.clfmethods{c}{d} = cellfun(@(x)(x.train(data{d},design{d})),obj.clfmethods{c}{d},'UniformOutput',false);
                   if c<obj.nmethods
                     data{d} = cellfun(@(x)(x.test(data{d})),obj.clfmethods{c}{d},'UniformOutput',false);
                   end
                 end
               else
                 error('cannot handle multiple datasets');
               end
               
             else
               
               % deal with ensemble methods (i.e., nested cell arrays)
               obj.clfmethods{c} = cellfun(@(x)(x.train(data,design)),obj.clfmethods{c},'UniformOutput',false);
               if c<obj.nmethods
                 data = cellfun(@(x)(x.test(data)),obj.clfmethods{c},'UniformOutput',false);
               end
             end
             
           else
             
             if iscell(data) && ~isa(obj.clfmethods{c},'transfer_learner')
               % if the method is not a transfer learner then we apply
               % the method to each dataset separately and convert
               % the method to a cell array
               
               m = cell(1,length(data));
               for d=1:length(data)
                 m{d} = obj.clfmethods{c}.train(data{d},design{d});
                 if c<obj.nmethods
                   data{d} = m{d}.test(data{d});
                 end
               end
               
               obj.clfmethods{c} = m;
               
             else
               
               obj.clfmethods{c} = obj.clfmethods{c}.train(data,design);
               if c<obj.nmethods
                 data = obj.clfmethods{c}.test(data);
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
         
         % cast to datasets if necessary
         
         if iscell(data)
           for c=1:length(data)
             if ~isa(data{c},'dataset')
               data{c} = dataset(data{c});
             end
           end
         elseif ~isa(data,'dataset')
           data = dataset(data);
         end
         
         for c=1:obj.nmethods
           
           if iscell(obj.clfmethods{c})
             % deal with ensemble methods
             
             if iscell(data) && ~isa(obj.clfmethods{c}{1},'transfer_learner')
               
               if length(data) == length(obj.clfmethods{c})
                 for d=1:length(data)
                   data{d} = obj.clfmethods{c}{d}.test(data{d});
                 end
               else
                 error('cannot handle multiple datasets');
               end
               
             else
               
               data = cellfun(@(x)(x.test(data)),obj.clfmethods{c},'UniformOutput',false);
             end
             
           else
             
             data = obj.clfmethods{c}.test(data);
             
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
         
         % cast to datasets if necessary
         
         if iscell(data)
           for c=1:length(data)
             if ~isa(data{c},'dataset')
               data{c} = dataset(data{c});
             end
           end
         elseif ~isa(data,'dataset')
           data = dataset(data);
         end
         
         for c=obj.nmethods:-1:1
           
           if iscell(obj.clfmethods{c})
             % deal with ensemble methods
             
             if iscell(data) && ~isa(obj.clfmethods{c}{1},'transfer_learner')
               
               if length(data) == length(obj.clfmethods{c})
                 for d=1:length(data)
                   data{d} = obj.clfmethods{c}{d}.untest(data{d});
                 end
               else
                 error('cannot handle multiple datasets');
               end
               
             else
               
               data = cellfun(@(x)(x.untest(data)),obj.clfmethods{c},'UniformOutput',false);
             end
             
           else
             
             data = obj.clfmethods{c}.untest(data);
             
           end
           
         end
         
       end
       
       function p = predict(obj,data)
         % calls test functions and converts posterior into classifications
         
         for c=1:(length(obj.clfmethods)-1)
           data = obj.clfmethods{c}.test(data);
         end
         
         p = obj.clfmethods{end}.predict(data);
       end
       
       function s = name(obj)
         % returns classification procedure as a string
         
         s = '{ ';
         for j=1:length(obj.clfmethods)
           s = [s class(obj.clfmethods{j}) ' '];
         end
         s = strcat(s,' }');
         
       end
       
       function [m,desc] = getmodel(obj)
         % return the model implied by this classification procedure
         
         % get the final model
         if iscell(obj.clfmethods{end})
           
           ntasks = length(obj.clfmethods{end});
           
           mm = cell(1,ntasks);
           for c=1:ntasks
             
             [mm{c},desc] = obj.clfmethods{end}{c}.getmodel();
             
           end
           
           m = cell(size(mm{1},1),ntasks);
           for c=1:ntasks
             for i=1:size(mm{1},1)
               m{i,c} = mm{c}{i};
             end
           end
           
         else
           [m,desc] = obj.clfmethods{end}.getmodel();
         end
         
       end
       
    end
    
end
