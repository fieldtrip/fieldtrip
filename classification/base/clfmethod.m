classdef clfmethod
%clfmethod base class for classification toolbox methods
%   
%   This base class contains common properties
%   which may be called by all classification methods
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: clfmethod.m,v $
%

    properties
      verbose = false;
    end
    
    methods            
      
      function m = getmodel(obj,label,dims)
        % default behaviour when we ask for a model (override in subclass)
        
        if obj.verbose
          fprintf('don`t know how to return model for object of type %s; returning empty model\n',class(obj));
        end
        
        m = [];
        
      end
      
      function plot(obj,data,design)
        % plot object, possibly wrt data and/or a design
        
        if nargin < 2, data = []; end
        if nargin < 3, design = []; end
        
        figure;
        obj.plot_data_distribution(data,design);
        
      end            
     
    end
    
    methods(Abstract)
        obj = train(obj,data,design);
        post = test(obj,data);
    end
    
    methods(Access = protected,Static)
    
      function plot_data_distribution(data,design)
        
        if ~isempty(data)
          
          if isempty(design), design = ones(size(data,1),1); end
          
          if size(data,2) == 1
          
            scatter(data(:,1),ones(size(data,1),1),50,design,'filled')            
            xlabel('feature 1');
            ylim([0.8 1.2]);
            set(gca,'YTick',[]);
            title('data distribution');
            
          elseif size(data,2) == 2
            
            scatter(data(:,1),data(:,2),50,design,'filled')
            xlabel('feature 1');
            ylabel('feature 2');
            title('data distribution');
          
          elseif size(data,2) == 3
            
            scatter3(data(:,1),data(:,2),data(:,3),50,design,'filled')
            xlabel('feature 1');
            ylabel('feature 2');
            zlabel('feature 3');            
            title('data distribution');
            
          else
            fprintf('too many dimensions to plot\n');
          end
                    
        end
            
      end
    
    end
      
end