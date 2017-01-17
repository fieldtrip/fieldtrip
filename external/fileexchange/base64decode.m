function output = base64decode(input)
%BASE64DECODE Decode Base64 string to a byte array.
%
%    output = base64decode(input)
%
% The function takes a Base64 string INPUT and returns a uint8 array
% OUTPUT. JAVA must be running to use this function. The result is always
% given as a 1-by-N array, and doesn't retrieve the original dimensions.
%
% See also base64encode

error(nargchk(1, 1, nargin));
error(javachk('jvm'));
if ischar(input), input = uint8(input); end

output = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(input), 'uint8')';

end

