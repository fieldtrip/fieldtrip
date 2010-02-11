classdef discrete_cpd < cpd
%DISCRETE_CPD abstract discrete conditional probability distribution class 
%
%   Each discrete cpd is specified by its discrete parents,
%   and its parameters. 
%
%   SEE ALSO:
%       multinomial_cpd
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: discrete_cpd.m,v $
%

    properties
       
        statenames; % optional state names
        
    end

    methods
        function obj = discrete_cpd(child,dparents)
            % constructor
            
            obj.child = child;
            obj.dparents = dparents;
            
        end
        function par = parents(obj)
            par = obj.dparents;
        end
    end

end