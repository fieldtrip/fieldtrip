classdef ft_mv_meta < ft_mv_method
%FT_MV_META abstract class 
%
% This class takes an MVA as property and performs some form of
% meta-analysis on it.
%
% Copyright (c) 2011, Marcel van Gerven
    
 properties

   mva % used multivariate analysis
    
  end

  methods

    function obj = ft_mv_meta(varargin)

     obj@ft_mv_method(varargin{:});

     if iscell(obj.mva)
       for c=1:length(obj.mva)
         if ~isa(obj.mva{c},'ft_mv_analysis')
           obj.mva{c} = ft_mv_analysis(obj.mva{c});
         end
       end
     else
       if ~isa(obj.mva,'ft_mv_analysis')
         obj.mva = ft_mv_analysis(obj.mva);
       end
     end
     
    end
    
  end

end