function varargout = ft_inside_vol(varargin)

% This function is a backward compatibility wrapper for existing MATLAB scripts
% that call a function that is not part of the FieldTrip toolbox any more.
%
% Please update your code to make it future-proof.

oldname = mfilename;
newname = 'ft_inside_headmodel';

ft_warning('%s is only a backward compatibility wrapper, which will soon be removed. Please call %s instead.', upper(oldname), upper(newname));

[varargout{1:nargout}] = feval(newname, varargin{:});
