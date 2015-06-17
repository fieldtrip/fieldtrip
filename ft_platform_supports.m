function tf=ft_platform_supports(what,varargin)
% return whether the current platform supports a specific capability
%
% Usage:
%   tf=ft_platform_supports(what)
%   tf=ft_platform_supports('matlabversion',min_version,max_version)
%
% The following values are allowed for the 'what' parameter:
%   value                           means that the following is supported:
%
%   'which-all'                     which(...,'all')
%   'onCleanup'                     onCleanup(...)
%   'int32_logical_operations'      bitand(a,b) with a, b of type int32
%   'graphics_objects'              graphics sysem is object-oriented
%   'libmx_c_interface'             libmx is supported through mex in the
%                                   C-language (recent Matlab versions only
%                                   support C++)
%   'program_invocation_name'       program_invocation_name() (GNU Octave)
%   'singleCompThread'              start Matlab with -singleCompThread
%   'nosplash'                                        -nosplash
%   'nodisplay'                                       -nodisplay
%   'nojvm'                                           -nojvm
%   'no-gui'                        start GNU Octave with --no-gui
%   'RandStream.setGlobalStream'    RandStream.setGlobalStream(...)
%   'RandStream.setDefaultStream'   RandStream.setDefaultStream(...)
%   'rng'                           rng(...)
%   'rand-state'                    rand('state')
%
%
if ~ischar(what)
    error('first argument must be a string');
end

switch what
    case 'matlabversion'
        tf=is_matlab() && matlabversion(varargin{:});

    case 'which-all'
        tf=is_matlab();

    case 'onCleanup'
        tf=is_octave() || matlabversion(7.8, Inf);

    case 'int32_logical_operations'
        % earlier version of Matlab don't support bitand (and similar)
        % operations on int32
        tf=is_octave() || ~matlabversion(-inf, '2012a');

    case 'graphics_objects'
        % introduced in Matlab 2014b, graphics is handled through objects;
        % previous versions use numeric handles
        tf=is_matlab() && matlabversion('2014b', Inf);

    case 'libmx_c_interface'
        % removed after 2013b
        %
        tf=matlabversion(-Inf, '2013b');

    case 'program_invocation_name'
        % Octave supports program_invocation_name, which returns the path
        % of the binary that was run to start Octave
        tf=is_octave();

    case 'singleCompThread'
        tf=is_matlab() && matlabversion(7.8, inf);

    case {'nosplash','nodisplay','nojvm'}
        % Only on Matlab
        tf=is_matlab();

    case 'no-gui'
        % Only on Octave
        tf=is_octave();

    case 'RandStream.setGlobalStream'
        tf=is_matlab() && matlabversion('2008b', '2011b');

    case 'RandStream.setDefaultStream'
        tf=is_matlab() && matlabversion('2012a', inf);

    case 'randomized_PRNG_on_startup'
        tf=is_octave() || ~matlabversion(-Inf,'7.3');

    case 'rng'
        % recent Matlab versions
        tf=is_matlab() && matlabversion('7.12',Inf);

    case 'rand-state'
        % GNU Octave
        tf=is_octave();

    otherwise
        error('unsupported value for first argument: %s', what);

end % switch

end % function

function tf=is_matlab()
    tf=~is_octave();
end

function tf=is_octave()
    persistent cached_tf;

    if isempty(cached_tf)
        cached_tf=logical(exist('OCTAVE_VERSION', 'builtin'));
    end

    tf=cached_tf;
end

