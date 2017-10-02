function [gradnew] = ft_convert_grad(grad, amplitude, distance, scaling)

% FT_CONVERT_GRAD is deprecated, please use ft_datatype_sens

% DEPRECATED by roboos on 1 July 2013
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=963 for more details
% this function can be completely removed medio 2014
warning('this function is deprecated, please use ft_datatype_sens');

% at the time of writing this code, this is only supported in the "upcoming" version, not the "latest/2011v2"
gradnew = ft_datatype_sens(grad, 'version', 'upcoming', 'amplitude', amplitude, 'distance', distance, 'scaling', scaling);
