function a = tcpread(sock, siz, type)

% TCPREAD reads the available data from an open TCP socket,
% concatenates it to an existing local buffer and subsequently formats
% a small section of that buffer into a 'int16' or whatever data type.
% This allows you to read-access a TCP port in small chuncks similar
% as "fread" without loosing any bytes due to the TCP buffer itself
% overrunning.
%
% Use as
%   dat = tcpread(sock, size, format)
% where the socket is previously opened using pnet and format is
% something like 'uint8', 'char', or 'double'.
%
% To read null-terminated strings of unknown length, you can use
%   dat = tcpread(sock, char(0), 'char')
% 
% To clear the persistent buffer that is maintained inside this
% function for subsequent calls, you should "clear tcpread".
% 
% See also PNET

% Copyright (C) 2008-2009, Robert Oostenveld 

persistent offset buf bufnull

if isempty(offset)
    offset  = 0;
end

if isempty(buf)
    buf = [];
end

if isempty(bufnull)
    bufnull = [];
end

% this for reading a null-terminated string of unknown length
parsenull = isequal(siz, char(0));

% read all the available new data from the TCP stream
new = pnet(sock, 'read', 2^16, 'uint8', 'native', 'noblock');
buf = cat(2, buf, new);

%  fprintf('offset = %d, numel(buf) = %d\n', offset, numel(buf));

if ~parsenull
    siz = double(siz);
    if (length(siz)==1)
        siz = [1 siz];
    end
end

if ~parsenull
    % determine how many bytes should be converted
    switch type
        case {'char'}
            n = prod(siz) * 1;
        case {'int8' 'uint8'}
            n = prod(siz) * 1;
        case {'int16' 'uint16'}
            n = prod(siz) * 2;
        case {'int32' 'uint32' 'float32' 'single'}
            n = prod(siz) * 4;
        case {'int64' 'uint64' 'float64' 'double'}
            n = prod(siz) * 8;
        otherwise
            ft_error('unsupported type');
    end
    % read additional data from the TCP stream until there is enough
    while numel(buf)<(offset+n)
        new = pnet(sock, 'read', 2^16, 'uint8', 'native', 'noblock');
        buf = cat(2, buf, new);
    end
else
    if numel(bufnull)~=numel(buf)
        % this is to speed up subsequent parsing of a string based on null-characters
        % it is likely that a whole sequence of channel labels has to be parsed
        bufnull = (~buf);
    end
    % the purpose is to return a single null-terminated string
    n = find(bufnull((offset+1):end), 1, 'first') - 1;
end

if parsenull
    b = buf(offset + (1:n));
    a = cast(b, 'char');
    % on the next call also skip the null-character
    offset = offset + n + 1;
elseif strcmp(type, 'char')
    b = buf(offset + (1:n));
    a = reshape(cast(b, type), siz);
    % on the next call skip the bytes that have been processed sofar
    offset = offset + n;
else
    b = buf(offset + (1:n));
    a = reshape(typecast(b, type), siz);
    % on the next call skip the bytes that have been processed sofar
    offset = offset + n;
end

if offset>2^14
    % remove the oldest part from the buffer
    buf    = buf((offset+1):end);
    offset = 0;
    % fprintf('offset = %d, numel(buf) = %d\n', offset, numel(buf));
end

