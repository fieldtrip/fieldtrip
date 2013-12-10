function y = base64decode(x)
%BASE64DECODE Perform base64 decoding on a string.
%
%   BASE64DECODE(STR) decodes the given base64 string STR.
%
%   Any character not part of the 65-character base64 subset set is silently
%   ignored.
%
%   This function is used to decode strings from the Base64 encoding specified
%   in RFC 2045 - MIME (Multipurpose Internet Mail Extensions).  The Base64
%   encoding is designed to represent arbitrary sequences of octets in a form
%   that need not be humanly readable.  A 65-character subset ([A-Za-z0-9+/=])
%   of US-ASCII is used, enabling 6 bits to be represented per printable
%   character.
%
%   See also BASE64ENCODE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-09-20 08:20:50 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

%   Modified by Guillaume Flandin, May 2008


% Perform the following mapping
%--------------------------------------------------------------------------
%   A-Z  ->  0  - 25         a-z  ->  26 - 51         0-9  ->  52 - 61
%   +    ->  62              /    ->  63              =    ->  64
%   anything else -> NaN

base64chars = NaN(1,256);
base64chars('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=') = 0:64;
x = base64chars(x);

% Remove/ignore any characters not in the base64 characters list or '='
%--------------------------------------------------------------------------

x = x(~isnan(x));

% Replace any incoming padding ('=' -> 64) with a zero pad
%--------------------------------------------------------------------------

if     x(end-1) == 64, p = 2; x(end-1:end) = 0;
elseif x(end)   == 64, p = 1; x(end) = 0;
else                   p = 0;
end

% Allocate decoded data array
%--------------------------------------------------------------------------

n = length(x) / 4;                               % number of groups
x = reshape(uint8(x), 4, n);                     % input data
y = zeros(3, n, 'uint8');                        % decoded data

% Rearrange every 4 bytes into 3 bytes
%--------------------------------------------------------------------------
%    00aaaaaa 00bbbbbb 00cccccc 00dddddd
%
% to form
%
%    aaaaaabb bbbbcccc ccdddddd

y(1,:) = bitshift(x(1,:), 2);                    % 6 highest bits of y(1,:)
y(1,:) = bitor(y(1,:), bitshift(x(2,:), -4));    % 2 lowest bits of y(1,:)

y(2,:) = bitshift(x(2,:), 4);                    % 4 highest bits of y(2,:)
y(2,:) = bitor(y(2,:), bitshift(x(3,:), -2));    % 4 lowest bits of y(2,:)

y(3,:) = bitshift(x(3,:), 6);                    % 2 highest bits of y(3,:)
y(3,:) = bitor(y(3,:), x(4,:));                  % 6 lowest bits of y(3,:)

% Remove any zero pad that was added to make this a multiple of 24 bits
%--------------------------------------------------------------------------

if p, y(end-p+1:end) = []; end

% Reshape to a row vector
%--------------------------------------------------------------------------

y = reshape(y, 1, []);
