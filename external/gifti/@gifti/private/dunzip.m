function M = dunzip(Z)
% DUNZIP - decompress gzipped stream of bytes
% FORMAT M = dzip(Z)
% Z  -  compressed variable to decompress (uint8 vector)
% M  -  decompressed output
% 
% See also DZIP

% Carefully tested, but no warranty; use at your own risk.
% Michael Kleder, Nov 2005
% Modified by Guillaume Flandin, May 2008

import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
a   = java.io.ByteArrayInputStream(Z);
b   = java.util.zip.InflaterInputStream(a);
isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
c   = java.io.ByteArrayOutputStream;
isc.copyStream(b,c);
M   = c.toByteArray;
