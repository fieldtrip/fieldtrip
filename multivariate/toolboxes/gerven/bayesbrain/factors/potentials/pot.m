classdef pot
%POT abstract potential class 
%
%   Each potential is specified by its continuous domain, its discrete domain,
%   and its parameters. The parameters are specified by a cell array where
%   each element specifies a parameter. The element is a multidimensional
%   array which is indexed by the discrete parent configurations. Parents
%   are indexed by the order with which they occur in a graphical model.
%
%   SEE ALSO:
%       continuous_pot
%       discrete_pot
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: pot.m,v $
%

   properties
    
       cdomain = [];    % continuous domain
       ddomain = [];    % discrete domain
       dsize   = 0;     % size of discrete domain
   end

   methods
       function obj = pot(cdomain,ddomain,dsize)
           % constructor
                     
           obj.cdomain = cdomain;
           obj.ddomain = ddomain;
           obj.dsize   = dsize;
           
       end
       function dom = domain(obj)
           % return variables in this potential
           dom = [obj.cdomain obj.ddomain];
       end       
       function len = length(obj)
           % return number of variables in this potential
           len = length([obj.cdomain obj.ddomain]);
       end        
   end
   
   methods (Abstract)
       pot = observe(pot, evidence)  % observe some values of the potential
       pot = mtimes(a,b)
       n = states(obj) % number of states
   end
   
end 
