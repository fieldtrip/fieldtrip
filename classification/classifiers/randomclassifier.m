classdef randomclassifier < classifier
%RANDOMCLASSIFIER returns a random classification
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: randomclassifier.m,v $
%

    methods
        function obj = randomclassifier(varargin)                  
           obj = obj@classifier(varargin{:});                      
        end
        
        function obj = train(obj,data,design)
            % does nothing
            
            [~,design] = obj.check_input(data,design);
            
            if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
            
        end
        
        function post = test(obj,data)
           
          data = obj.check_input(data);
          
          % random classification
          post = rand(size(data,1),obj.nclasses);
          post = post ./ repmat(sum(post,2),[1 size(post,2)]);
          
        end

    end
end 
