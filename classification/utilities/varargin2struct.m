function cfg = varargin2struct(va)
% VARARGIN2STRUCT transform varargin to a structure
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: varargin2struct.m,v $
%

    cfg = [];

    for i=1:2:length(va)
        cfg.(va{i}) = va{i+1};
    end
    
end