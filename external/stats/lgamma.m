function varargout = lgamma(varargin)
varargout{1:nargout} = gammaln(varargin{:});
end