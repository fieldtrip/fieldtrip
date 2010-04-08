classdef param < handle
%PARAM parameter class
%
%  Derived from handle and contains a cell array; used to represent parameters and
%  allows sharing of parameters (i.e., equivalence classes).
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: param.m,v $
%

   properties
       value     % can be anything 
   end

   methods
       function obj = param(value)

           obj.value = value;
   
       end
   end
end 
