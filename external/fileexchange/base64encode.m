function output = base64encode(input)
%BASE64ENCODE Encode a byte array using Base64 codec.
%
%    output = base64encode(input)
%
% The function takes a char, int8, or uint8 array INPUT and returns Base64
% encoded string OUTPUT. JAVA must be running to use this function. Note
% that encoding doesn't preserve input dimensions.
%
% See also base64decode

error(nargchk(1, 1, nargin));
error(javachk('jvm'));
if ischar(input), input = uint8(input); end

output = char(org.apache.commons.codec.binary.Base64.encodeBase64Chunked(input))';

end
