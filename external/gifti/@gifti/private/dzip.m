function Z = dzip(M)
% DZIP - losslessly compress data using gzip
% FORMAT Z = dzip(M)
% M  -  byte stream to compress
% Z  -  compressed output
%
% See also DUNZIP

% This function uses the public domain ZLIB Deflater algorithm.
% Carefully tested, but no warranty; use at your own risk.
% Michael Kleder, Nov 2005
% Modified by Guillaume Flandin, May 2008

M = typecast(M(:),'uint8');
f = java.io.ByteArrayOutputStream();
%f = javaObject('java.io.ByteArrayOutputStream');
g = java.util.zip.DeflaterOutputStream(f);
%g = javaObject('java.util.zip.DeflaterOutputStream',f);
g.write(M);
g.close;
Z = typecast(f.toByteArray,'uint8');
f.close;
