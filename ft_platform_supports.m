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

