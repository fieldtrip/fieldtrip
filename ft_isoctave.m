function tf=ft_isoctave()
% return true if and only if the platform is GNU/Octave
%
% tf=moxunit_util_platform_is_octave()
%
% Output:
%   tf      true if the platform on which this function is run is
%           GNU/Octave, false otherwise (typically, if the platform is
%           Matlab)
%
% NNO Jan 2014

    persistent cached_tf;

    if islogical(cached_tf)
        tf=cached_tf;
        return;
    end

    tf=logical(exist('OCTAVE_VERSION', 'builtin'));
    cached_tf=tf;
