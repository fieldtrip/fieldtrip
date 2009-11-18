classdef ptsner < preprocessor
%PTSNER applies parametric t-SNE dimensionality reduction algorithm
%
%   NOTE: 
%     - This preprocessor currently only works for discrete input data!
%
%   Copyright (c) 2009, Marcel van Gerven, Laurens van der Maaten
%
%   $Log: ptsner.m,v $
%

    properties
      
      layers = [100 100 200 2]; %[500 500 2000 2]; % layer structure of underlying RBM

      network; % network structure
      
    end

    methods
    
        function obj = ptsner(varargin)
           
            obj = obj@preprocessor(varargin{:});                        
            
        end
        
        function obj = train(obj,data,design)
                  
          % perform tsne on the training data to get an initial solution
          if iscell(data)

           for c=1:length(data)
             obj.network{c} = train_par_tsne(data{c},design{c},data{c}(1,:),design{c}(1,1),obj.layers, 'CD1');
           end
           
          else
    
              % Train the parametric t-SNE network
              obj.network = train_par_tsne(data, design(:,1), data(1,:), design(1,1), obj.layers, 'CD1');
            
          end

        end
        
        function data = test(obj,data)

          % Construct embedding
          if iscell(data)
            
            for c=1:length(data)
              data{c} = run_data_through_network(obj.network{c}, data{c});
            end
            
          else            
            data = run_data_through_network(obj.network, data);            
          end
          
        end
        
    end
end 
