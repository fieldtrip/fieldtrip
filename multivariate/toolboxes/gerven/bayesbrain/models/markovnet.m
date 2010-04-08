classdef markovnet < graphicalmodel
%MARKOVNET Markov network class
%
%   A Markov network is constructed using
%
%   mn = mnet(factors,varargin)
%
%   where the factors are potentials
%
%   SEE ALSO:
%   graphicalmodel.m
%   potential.m
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: markovnet.m,v $
%

   methods
       function obj = markovnet(factors,varargin)
              
           assert(iscell(factors));
           
           assert(all(cellfun(@(x)(isa(x,'pot')),factors)));
           
           obj = obj@graphicalmodel(factors,varargin{:});

           % construct the graph from the factors
           
           g = sparse(obj.length(),obj.length());
           for i=1:obj.length()

               clqdom = obj.factors{i}.domain;
               if length(clqdom) > 1
                   g(clqdom,clqdom) = true;
               end
           end
           for i=1:obj.length(), g(i,i) = 0; end
           obj.g = graph(g);
           
       end
   end
end 
