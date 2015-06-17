function tf=ft_platform_supports(what,varargin)
% return whether the current platform supports a specific capability

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


    otherwise
        error('unsupported value for first argument: %s', what);

end % switch

end % function

function tf=environment_is_matlab()
    tf=~environment_is_octave();
end

function tf=environment_is_octave()
    persistent cached_tf;

    if isempty(cached_tf)
        cached_tf=logical(exist('OCTAVE_VERSION', 'builtin'));
    end

    tf=cached_tf;
end

