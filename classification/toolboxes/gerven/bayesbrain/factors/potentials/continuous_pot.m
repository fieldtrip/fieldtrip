classdef continuous_pot < pot
% CONTINUOUS_POT continuous potential class 
%
%   SEE ALSO:
%       canonical_pot
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: continuous_pot.m,v $
%

    methods
        function obj = continuous_pot(cdomain,ddomain,dsize)
            % constructor

            obj = obj@pot(cdomain,ddomain,dsize);
        end
        function n = states(obj)
            % number of states of a continuous pot is 1 by default
            n = 1;
        end
    end
end