classdef continuous_cpd < cpd
%CONTINUOUS_CPD abstract continuous conditional probability distribution class 
%
%   Each continuous cpd is specified by its continuous parents, discrete parents,
%   and its parameters. 
%
%   SEE ALSO:
%       gaussian_cpd
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: continuous_cpd.m,v $
%
    
    methods
        function obj = continuous_cpd(child,cparents,dparents)
            % constructor
            
            obj.child = child;
            obj.cparents = cparents;
            obj.dparents = dparents;
            
        end
     
        function n = states(obj)
            % number of states of a continuous cpd is 1 by default
            n = 1;
        end
    end

end