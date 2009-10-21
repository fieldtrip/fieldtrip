function varargout = filetype_check_uri(filename, ftyp)

% FILETYPE_CHECK_URI
%
% Supported URIs are
%   shm://<filename>
%   fifo://<filename>
%   buffer://<host>:<port>
%   tcp://<host>:<port>
%   udp://<host>:<port>
%   mysql://<user>:<password>@<host>:<port>
%   rfb://<password>@<host>:<port>
%   serial:<port>?key1=value1&key2=value2&...
%
% The URI schemes supproted by these function are not the official schemes.
% See the documentation included inside this function for more details.
% RFC4395 defines an IANA-maintained registry of URI Schemes. See also
% http://www.iana.org/assignments/uri-schemes.html and
% http://en.wikipedia.org/wiki/URI_scheme#Generic_syntax.

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: filetype_check_uri.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.7  2008/12/19 14:40:23  marvger
% added support for udp, tcp and fifo
%
% Revision 1.6  2008/02/19 10:07:34  roboos
% replaced tcpsocket with buffer
%
% Revision 1.5  2007/11/07 10:45:56  roboos
% made port optional for mysql
%
% Revision 1.4  2007/10/15 16:02:46  roboos
% added rfb://<password>@<host>:<port>
%
% Revision 1.3  2007/07/27 12:18:50  roboos
% added ctf_shm
%
% Revision 1.2  2007/06/19 11:13:59  chrhes
% small change in how the serial port case is processed to allow for the
% possibility of no options being passed
%
% Revision 1.1  2007/06/06 07:10:36  roboos
% initial implementation
%

sep = find(filename==':');
if ~isempty(sep)
  scheme = filename(1:(sep(1)-1));
else
  scheme = [];
end

if nargin>1
  % only return true or false
  varargout{1} = strcmp(scheme, ftyp);
else
  % return the full details of this URI scheme
  switch scheme
    case 'shm'
      % shm://<filename>
      % the filename is optional, usually it can and will be read from shared memory 
      if length(filename)>6
        varargout{1} = filename(7:end);
      else
        varargout{1} = [];
      end

    case 'fifo'
      % fifo://<filename>
      varargout{1} = filename(8:end);

    case 'buffer'
      % buffer://<host>:<port>
      tok = tokenize(filename(10:end), ':');
      varargout{1} = tok{1};
      varargout{2} = str2num(tok{2});

    case 'tcp'
      % tcp://<host>:<port>
      tok = tokenize(filename(7:end), ':');
      varargout{1} = tok{1};
      varargout{2} = str2num(tok{2});

    case 'udp'
      % udp://<host>:<port>
      tok = tokenize(filename(7:end), ':');
      varargout{1} = tok{1};
      varargout{2} = str2num(tok{2});

    case 'rfb'
      % rfb://<password>@<host>:<port>
      tok0 = tokenize(filename(7:end), '@');
      tok1 = tokenize(tok0{2}, ':');
      varargout{1} = tok0{1};
      varargout{2} = tok1{1};
      varargout{3} = str2num(tok1{2});

    case 'mysql'
      % mysql://<user>:<password>@<host>:<port>
      tok0 = tokenize(filename(9:end), '@');
      tok1 = tokenize(tok0{1}, ':');
      tok2 = tokenize(tok0{2}, ':');
      varargout{1} = tok1{1};
      varargout{2} = tok1{2};
      varargout{3} = tok2{1};
      if length(tok2)>1
        varargout{4} = str2num(tok2{2});
      else
        varargout{4} = [];
      end

    case 'serial'
      % serial:<Port>?key1=value1&key2=value2&... 
      % the supported optional arguments are
      %   BaudRate
      %   DataBits
      %   DataTerminalReady
      %   FlowControl
      %   Parity
      %   Port
      %   ReadAsyncMode
      %   RequestToSend
      %   StopBits
      %   Terminator
      tok0 = tokenize(filename(8:end), '?');
      varargout{1} = tok0{1};
      varargout{2} = [];
      if length(tok0)>1
        tok1 = tokenize(tok0{2}, '&');
        for i=1:length(tok1)
          tok2 = tokenize(tok1{i}, '=');
          opt{(i-1)*2+1} = tok2{1};
          opt{(i-1)*2+2} = tok2{2};
        end
        varargout{2} = opt;
      end
      

    otherwise
      error('unsupported scheme in URI')
  end
end

% RFC4395 defines an IANA-maintained registry of URI Schemes.
% see also http://www.iana.org/assignments/uri-schemes.html
% and http://en.wikipedia.org/wiki/URI_scheme#Generic_syntax
%
% Scheme Description Reference
% -------------------------------------------------------------
% aaa	Diameter Protocol	[RFC3588]
% aaas	Diameter Protocol with Secure Transport	[RFC3588]
% acap	application configuration access protocol	[RFC2244]
% cap	Calendar Access Protocol	[RFC4324]
% cid	content identifier	[RFC2392]
% crid	TV-Anytime Content Reference Identifier	[RFC4078]
% data	data	[RFC2397]
% dav	dav	[RFC-ietf-webdav-rfc2518bis-18.txt]
% dict	dictionary service protocol	[RFC2229]
% dns	Domain Name System	[RFC4501]
% fax	fax	[RFC3966]
% file	Host-specific file names	[RFC1738]
% ftp	File Transfer Protocol	[RFC1738]
% go	go	[RFC3368]
% gopher	The Gopher Protocol	[RFC4266]
% h323	H.323	[RFC3508]
% http	Hypertext Transfer Protocol	[RFC2616]
% https	Hypertext Transfer Protocol Secure	[RFC2818]
% icap	Internet Content Adaptation Protocol	[RFC3507]
% im	Instant Messaging	[RFC3860]
% imap	internet message access protocol	[RFC2192]
% info	Information Assets with Identifiers in Public Namespaces	[RFC4452]
% ipp	Internet Printing Protocol	[RFC3510]
% iris	Internet Registry Information Service	[RFC3981]
% iris.beep	iris.beep	[RFC3983]
% iris.xpc	iris.xpc	[RFC-ietf-crisp-iris-xpc-06.txt]
% iris.xpcs	iris.xpcs	[RFC-ietf-crisp-iris-xpc-06.txt]
% iris.lwz	iris.lwz	[RFC-ietf-crisp-iris-lwz-08.txt]
% ldap	Lightweight Directory Access Protocol	[RFC4516]
% mailto	Electronic mail address	[RFC2368]
% mid	message identifier	[RFC2392]
% modem	modem	[RFC3966]
% msrp	Message Session Relay Protocol	[RFC-ietf-simple-message-sessions-19.txt]
% msrps	Message Session Relay Protocol Secure	[RFC-ietf-simple-message-sessions-19.txt]
% mtqp	Message Tracking Query Protocol	[RFC3887]
% mupdate	Mailbox Update (MUPDATE) Protocol	[RFC3656]
% news	USENET news	[RFC1738]
% nfs	network file system protocol	[RFC2224]
% nntp	USENET news using NNTP access	[RFC1738]
% opaquelocktoken	opaquelocktokent	[RFC-ietf-webdav-rfc2518bis-18.txt]
% pop	Post Office Protocol v3	[RFC2384]
% pres	Presence	[RFC3859]
% rtsp	real time streaming protocol	[RFC2326]
% service	service location	[RFC2609]
% shttp	Secure Hypertext Transfer Protocol	[RFC2660]
% sip	session initiation protocol	[RFC3261]
% sips	secure session initiation protocol	[RFC3261]
% snmp	Simple Network Management Protocol	[RFC4088]
% soap.beep	soap.beep	[RFC3288]
% soap.beeps	soap.beeps	[RFC3288]
% tag	tag	[RFC4151]
% tel	telephone	[RFC3966]
% telnet	Reference to interactive sessions	[RFC4248]
% tftp	Trivial File Transfer Protocol	[RFC3617]
% thismessage	multipart/related relative reference resolution	[RFC2557]
% tip	Transaction Internet Protocol	[RFC2371]
% tv	TV Broadcasts	[RFC2838]
% urn	Uniform Resource Names (click for registry)	[RFC2141]
% vemmi	versatile multimedia interface	[RFC2122]
% xmlrpc.beep	xmlrpc.beep	[RFC3529]
% xmlrpc.beeps	xmlrpc.beeps	[RFC3529]
% xmpp	Extensible Messaging and Presence Protocol	[RFC4622]
% z39.50r	Z39.50 Retrieval	[RFC2056]
% z39.50s	Z39.50 Session	[RFC2056]
