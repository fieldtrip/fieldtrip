function y=afni_swapbytes(x)
% Perform byteswap for input data of type int, single or double
%
% (Octave 3.8 does not support a byteswap for data of class 'single')

    c=class(x);

    switch c
        case {'int8','uint8'}
            % single byte, no swap needed
            y=x;
            return;

        case {'int16','uint16'}
            nbytes=2;

        case {'int32','uint32','single'}
            nbytes=4;

        case {'int64','uint64','double'}
            nbytes=8;

        otherwise
            error('Cannot swap bytes for input of type ''%s''', c);
    end

    % convert x to bytes, with the k-th column the byte representation of
    % the k-th element of x
    x_bytes=reshape(typecast(x(:),'uint8'),nbytes,[]);

    % swap each columna and convert back to original format
    y=reshape(typecast(reshape(x_bytes(end:-1:1,:),[],1),c),size(x));

