function y = base64encode(x, eol)
%BASE64ENCODE Perform base64 encoding on a string.
%
%   BASE64ENCODE(STR, EOL) encode the given string STR.  EOL is the line ending
%   sequence to use; it is optional and defaults to '\n' (ASCII decimal 10).
%   The returned encoded string is broken into lines of no more than 76
%   characters each, and each line will end with EOL unless it is empty.  Let
%   EOL be empty if you do not want the encoded string broken into lines.
%
%   STR and EOL don't have to be strings (i.e., char arrays).  The only
%   requirement is that they are vectors containing values in the range 0-255.
%
%   This function may be used to encode strings into the Base64 encoding
%   specified in RFC 2045 - MIME (Multipurpose Internet Mail Extensions).  The
%   Base64 encoding is designed to represent arbitrary sequences of octets in a
%   form that need not be humanly readable.  A 65-character subset
%   ([A-Za-z0-9+/=]) of US-ASCII is used, enabling 6 bits to be represented per
%   printable character.
%
%   Examples
%   --------
%
%   If you want to encode a large file, you should encode it in chunks that are
%   a multiple of 57 bytes.  This ensures that the base64 lines line up and
%   that you do not end up with padding in the middle.  57 bytes of data fills
%   one complete base64 line (76 == 57*4/3):
%
%   If ifid and ofid are two file identifiers opened for reading and writing,
%   respectively, then you can base64 encode the data with
%
%      while ~feof(ifid)
%         fwrite(ofid, base64encode(fread(ifid, 60*57)));
%      end
%
%   or, if you have enough memory,
%
%      fwrite(ofid, base64encode(fread(ifid)));
%
%   See also BASE64DECODE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-02-03 21:36:56 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam


   % make sure we have the EOL value
   if nargin < 2
      eol = ''; %sprintf('\n');
   else
      if sum(size(eol) > 1) > 1
         error('EOL must be a vector.');
      end
      if any(eol(:) > 255)
         error('EOL can not contain values larger than 255.');
      end
   end

   if sum(size(x) > 1) > 1
      error('STR must be a vector.');
   end

   x   = uint8(x);
   eol = uint8(eol);

   ndbytes = length(x);                 % number of decoded bytes
   nchunks = ceil(ndbytes / 3);         % number of chunks/groups
   nebytes = 4 * nchunks;               % number of encoded bytes

   % add padding if necessary, to make the length of x a multiple of 3
   if rem(ndbytes, 3)
      x(end+1 : 3*nchunks) = 0;
   end

   x = reshape(x, [3, nchunks]);        % reshape the data
   y = repmat(uint8(0), 4, nchunks);    % for the encoded data

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Split up every 3 bytes into 4 pieces
   %
   %    aaaaaabb bbbbcccc ccdddddd
   %
   % to form
   %
   %    00aaaaaa 00bbbbbb 00cccccc 00dddddd
   %
   y(1,:) = bitshift(x(1,:), -2);                  % 6 highest bits of x(1,:)

   y(2,:) = bitshift(bitand(x(1,:), 3), 4);        % 2 lowest bits of x(1,:)
   y(2,:) = bitor(y(2,:), bitshift(x(2,:), -4));   % 4 highest bits of x(2,:)

   y(3,:) = bitshift(bitand(x(2,:), 15), 2);       % 4 lowest bits of x(2,:)
   y(3,:) = bitor(y(3,:), bitshift(x(3,:), -6));   % 2 highest bits of x(3,:)

   y(4,:) = bitand(x(3,:), 63);                    % 6 lowest bits of x(3,:)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Now perform the following mapping
   %
   %   0  - 25  ->  A-Z
   %   26 - 51  ->  a-z
   %   52 - 61  ->  0-9
   %   62       ->  +
   %   63       ->  /
   %
   % We could use a mapping vector like
   %
   %   ['A':'Z', 'a':'z', '0':'9', '+/']
   %
   % but that would require an index vector of class double.
   %
   z = repmat(uint8(0), size(y));
   i =           y <= 25;  z(i) = 'A'      + double(y(i));
   i = 26 <= y & y <= 51;  z(i) = 'a' - 26 + double(y(i));
   i = 52 <= y & y <= 61;  z(i) = '0' - 52 + double(y(i));
   i =           y == 62;  z(i) = '+';
   i =           y == 63;  z(i) = '/';
   y = z;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Add padding if necessary.
   %
   npbytes = 3 * nchunks - ndbytes;     % number of padding bytes
   if npbytes
      y(end-npbytes+1 : end) = '=';     % '=' is used for padding
   end

   if isempty(eol)

      % reshape to a row vector
      y = reshape(y, [1, nebytes]);

   else

      nlines = ceil(nebytes / 76);      % number of lines
      neolbytes = length(eol);          % number of bytes in eol string

      % pad data so it becomes a multiple of 76 elements
      y(nebytes + 1 : 76 * nlines) = 0;
      y = reshape(y, 76, nlines);

      % insert eol strings
      eol = eol(:);
      y(end + 1 : end + neolbytes, :) = eol(:, ones(1, nlines));

      % remove padding, but keep the last eol string
      m = nebytes + neolbytes * (nlines - 1);
      n = (76+neolbytes)*nlines - neolbytes;
      y(m+1 : n) = '';

      % extract and reshape to row vector
      y = reshape(y, 1, m+neolbytes);

   end

   % output is a character array
   y = char(y);
