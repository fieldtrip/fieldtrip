classdef cpd
%CPD abstract conditional probability distribution class 
%
%   Each cpd is specified by its continuous parents, its discrete parents,
%   and its parameters. 
%
%   SEE ALSO
%    discrete_cpd
%    continuous_cpd
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: cpd.m,v $
%

    properties

        name;             % optional name of this cpd
        
        child = 0;        % index of this cpd
        cparents = [];    % continuous parents
        dparents = [];    % discrete parents      
        
        ess;     % expected sufficient statistics
 
    end
    
    methods
        function dom = domain(obj)
           dom = [obj.child obj.cparents obj.dparents]; 
        end
        function par = parents(obj)
            par = [obj.cparents obj.dparents];
        end
    end

    methods (Abstract)
        pot = cpd2pot(obj) % transforms cpd to potential
        sz  = dsize(obj) % size of the discrete part
        n   = states(obj) % number of states
    end
end