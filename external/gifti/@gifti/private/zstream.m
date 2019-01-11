function Z = zstream(action,D)
% Compress/decompress stream of bytes using Deflate/Inflate
% FORMAT Z = zstream('C',D)
% D        - data stream to compress (converted to uint8 if needed)
% Z        - compressed data stream (uint8 vector)
% FORMAT D = zstream('D',Z)
% Z        - data stream to decompress (uint8 vector)
% D        - decompressed data stream (uint8 vector)
%__________________________________________________________________________
%
% This C-MEX file relies on:
% * miniz, by Rich Geldreich
%   http://code.google.com/p/miniz/
% Fallback Java implementation is adapted from:
% * dzip/dunzip, by Michael Kleder
%   http://www.mathworks.com/matlabcentral/fileexchange/8899
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: zstream.m 6417 2015-04-21 16:03:44Z guillaume $


if exist('OCTAVE_VERSION','builtin')
    error('zstream.c not compiled - see Makefile');
end

switch upper(action)
    case 'C'
        D = typecast(D(:),'uint8');
        f = java.io.ByteArrayOutputStream();
        g = java.util.zip.DeflaterOutputStream(f);
        g.write(D);
        g.close;
        Z = typecast(f.toByteArray,'uint8');
        f.close;
        
    case 'D'
        import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
        a   = java.io.ByteArrayInputStream(D);
        b   = java.util.zip.InflaterInputStream(a);
        isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
        c   = java.io.ByteArrayOutputStream;
        isc.copyStream(b,c);
        Z   = c.toByteArray;
        
    otherwise
        error('Unknown action "%s".',action);
end
