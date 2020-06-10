function cmap = ft_colormap(varargin)

% FT_COLORMAP is a wrapper around the MATLAB COLORMAP function. It has the same
% usage as COLORMAP, but also knows about the colormaps returned by BREWERMAP.
% The latter can be specified as a string, e.g. 'RdBu', or as a 2-element cell,
% e.g. {'RdBu', 15} or {'*RdBu', 15}, where the second element specifies the number of colors. See BREWERMAP for more information.

switch nargin
  case 0
    cmap = colormap;
  case 1
    if ishandle(varargin{1})
      cmap = colormap(varargin{1});
    else
      if iscell(varargin{1})
        assert(numel(varargin{1})==2&&ischar(varargin{1}{1})&&isscalar(varargin{1}{2}));
         
        % this is assumed to be a pair of arguments, where the first one specifies the brewermap scheme,
        % and the second scalar the number of colors
        ft_hastoolbox('brewermap', 1);
        cmap = brewermap(varargin{1}{2}, varargin{1}{1});
      elseif ischar(varargin{1})
        % this can be a situation, where the string argument is either a permitted string for MATLAB's colormap, 
        % or for the brewermap, explicitly check the brewermap ones, assume a permitted string for colormap otherwise
        ft_hastoolbox('brewermap', 1);
        list = repmat(brewermap('list'), [2 1]);
        for k = 1:numel(list)/2, list{k} = sprintf('*%s',list{k}); end
                
        if any(strcmp(list, varargin{1}))
          cmap = brewermap(64, varargin{1});
        else    
          cmap = colormap(varargin{1});
        end
      else
        cmap = colormap(varargin{1});
      end
    end
  case 2
    if ishandle(varargin{1}) 
      cmap = ft_colormap(varargin{2});
    else
      ft_error('unexpected input to ft_colormap');
    end
  otherwise
    ft_error('wrong number of input arguments for ft_colormap');
end

if nargin==2
  colormap(varargin{1}, cmap);
else
  colormap(cmap);
end

