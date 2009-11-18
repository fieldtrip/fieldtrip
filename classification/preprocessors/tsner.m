classdef tsner < preprocessor
%TSNER applies t-SNE dimensionality reduction algorithm
%
%   NOTE: 
%     - mapping between spaces is learnt post hoc using backpropagation     
%
%   Copyright (c) 2009, Marcel van Gerven, Laurens van der Maaten
%
%   $Log: tsner.m,v $
%

    properties
      
      dims = 3; % final dimensions
      initial_dims = 30; % initial dimensions
      perplexity = 10; % perplexity during training
                  
      initial_solution % initial solution for train data
      traindesign % training labels
      
      net % trained backprop network
      nhidden = 10; % number of hidden units
            
      fast = false; % try fast tsne
      landmarks = []; % landmarks to use
      path = '~/code'; % path to classification module
      
    end

    methods
    
        function obj = tsner(varargin)
           
            obj = obj@preprocessor(varargin{:});                        
            
            if obj.fast
              path = input('for fast tsne, please specify the path to the classification module: ');
              if ~isempty(path)
                obj.path = path;
              end
            end
            
        end
        
        function obj = train(obj,data,design)

          obj.traindesign = design;
                    
          % perform tsne on the training data to get an initial solution
          if iscell(data)

            obj.initial_solution = cell(1,length(data));
            for c=1:length(data)

              if obj.verbose
                labels = design{c}(:,1);
              else
                labels = [];
              end

              if obj.fast
                error('not implemented yet');
              else
                obj.initial_solution{c} = tsne(data{c},labels,obj.dims,obj.initial_dims,obj.perplexity);
              end
              
              % train neural network
              obj.net{c} = newff(data{c}',obj.initial_solution{c}',obj.nhidden);
              obj.net{c} = train(obj.net{c},data{c}',obj.initial_solution{c}');
              nntraintool('close');
              
            end

          else

            if obj.verbose
              labels = design(:,1);
            else
              labels = [];
            end
            
            if obj.fast 
              
              % we must jump to the directory and then jump back
              p = pwd;
              cd([obj.path '/classification/toolboxes/maaten/tsne']);
              obj.initial_solution = fast_tsne(data, obj.dims, obj.initial_dims, obj.landmarks, obj.perplexity);
              cd(p);
              
            else
              obj.initial_solution = tsne(data,labels,obj.dims,obj.initial_dims,obj.perplexity);
            end
  
            % train neural network
            obj.net = newff(data',obj.initial_solution',obj.nhidden);
            obj.net = train(obj.net,data',obj.initial_solution');
            nntraintool('close');
          
          end
          
        end
        
        function data = test(obj,data)
          
          if iscell(data)

            for c=1:length(data)

              data{c} = transpose(sim(obj.net{c},data{c}'));
  
              if obj.verbose
                
                X = [obj.initial_solution{c}; data{c}];
                labels = [obj.traindesign{c}; zeros(size(data{c},1),1)];
                
                if obj.dims == 1
                  scatter(X, X, 40, labels, 'filled');
                elseif obj.dims == 2
                  scatter(X(:,1), X(:,2), 40, labels, 'filled');
                elseif obj.dims == 3
                  scatter3(X(:,1), X(:,2), X(:,3), 40, labels, 'filled');
                end
                if obj.dims < 4
                  axis tight
                  axis off
                  drawnow
                end
              end
              
            end
            
          else

            data = transpose(sim(obj.net,data'));
  
            if obj.verbose
            
              X = [obj.initial_solution; data];
              labels = [obj.traindesign; zeros(size(data,1),1)];
              
              if obj.dims == 1
                scatter(X, X, 40, labels, 'filled');
              elseif obj.dims == 2
                scatter(X(:,1), X(:,2), 40, labels, 'filled');
              elseif obj.dims == 3
                scatter3(X(:,1), X(:,2), X(:,3), 40, labels, 'filled');
              end
              if obj.dims < 4
                axis tight
                axis off
                drawnow
              end
            end
            
          end
                    
        end
        
    end
end 
