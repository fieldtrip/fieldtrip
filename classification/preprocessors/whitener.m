classdef whitener < preprocessor
%WHITENER whitens the data
%
%   data is assumed to be standardized
%
%   Copyright (c) 2009, Marcel van Gerven
%

    properties
      
      rdim; % number of components
      
      wmat;
           
    end

    methods
    
        function obj = whitener(varargin)
           
          obj = obj@preprocessor(varargin{:});
          
        end
        
        function obj = train(obj,data,design)

          if isempty(obj.rdim)
            obj.rdim = data.nfeatures;
          end
          
          obj.wmat = dataset.whitening_transform(data.X,obj.rdim);
                    
        end
        
        function data = test(obj,data)

          data = dataset(data.X*obj.wmat');
           
        end
       
    end
end
