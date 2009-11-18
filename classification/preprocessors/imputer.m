classdef imputer < preprocessor
 %IMPUTER handles missing data in various ways
%
%   all data columns containing just nans will be removed
%
%   Options:
%   'imputation' = imputation method []:
%                   []            : no imputation
%                   scalar value  : impute that value
%                   'rand'        : impute a random value between 0 and 1
%                   'randn'       : impute a standard normal random value 
%
%   TODO:
%   other imputation methods
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: imputer.m,v $
%

    properties
    
      nancols % the columns containing only nans as a bitvector
      
      imputation % how to impute values
      
    end    

    methods
        function obj = imputer(varargin)
            
            obj = obj@preprocessor(varargin{:});                        
            
        end
        function obj = train(obj,data,design)
        
          % determine columns with just nan values
          if iscell(data)
            obj.nancols = cell(1,length(data));
            for c=1:length(data)
              obj.nancols{c} = isnan(sum(data{c},1));
            end

            if obj.verbose && ~isempty(obj.nancols{1})
              fprintf('removing features with fully missing data\n');
            end

          else
              obj.nancols = isnan(sum(data,1));  
              
              if obj.verbose && ~isempty(obj.nancols)
                fprintf('removing features with fully missing data\n');
              end
                        
          end
                    
        end
        
        function data = test(obj,data)                        
            
          if iscell(data)
            
            cnancols = obj.nancols;
            for c=1:length(data)
              obj.nancols = cnancols{c};
              data{c} = obj.test(data{c});
            end
            obj.nancols = cnancols;
          
          else  
            
            data = data(:,~obj.nancols);

            if ~isempty(obj.imputation)
              
                nanidx = isnan(data(:)); 
           
                if ~isempty(nanidx)
                  nnan = length(nanidx);
                  
                  if isscalar(obj.imputation)
                    data(nanidx) = obj.imputation;
                  else
                    switch obj.imputation
                      case 'rand'
                        data(nanidx) = rand(nnan,1);
                      case 'randn'
                        data(nanidx) = randn(nnan,1);
                    end
                  end
                end
            end
          end
          
        end
        
    end
end 
