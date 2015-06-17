function tf=ft_platform_supports(what)
% return whether the current platform supports a specific capability

if ~ischar(what)
    error('first argument must be a string');
end

switch what
   
    otherwise
        error('unsupported value for first argument: %s', what);
    
end % switch


end % function

