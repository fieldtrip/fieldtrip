function response = json2couch(jsonfile, couchdburl, dbname, docname, options)
%
%    json2couch(jsonfile, servername, dbname, docname, options)
%
%    uploading JSON-encoded data to a CouchDB database (a NoSQL database)
%    as a document
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        jsonfile: the path to the .json file
%        couchdburl: the URL of the CouchDB server, usually it is
%                  http://servername:5984 where servername is your own
%                  server's domain name
%        dbname: the database for whcih the file is uploaded to, must be
%                  created first
%        docname: the document name for the uploaded JSON data
%        options: a options structure created by weboptions(), defining
%                 Username, Password, ContentType when such information is
%                 desired to access the server; by default,
%                 options.ContentType is set to 'json' and
%                 options.RequestMethod is set to 'POST'
%
%                 if options is a string, it defines a template for a curl
%                 command, for example 'curl -X PUT -d @$f $u/$d/$s'.
%                 Variables (start with $) are expanded as
%                   $f -> path of the json file (first input)
%                   $u -> couchdb full URL (second input)
%                   $d -> database name (third input)
%                   $s -> document name (forth input)
%
%    examples:
%        values = inputdlg({'Username:', 'Password:'});
%        options = weboptions('ContentType', 'json', 'RequestMethod', 'POST', 'Username',values{1},'Password',values{2});
%        json2couch('ds001.json', 'https://example.com:5984', 'bids-samples', 'ds001', options)
%        json2couch('ds001.json', sprintf('https://%s:%s@example.com:5984', values{1}, values{2}), 'bids-samples', 'ds001', 'curl -X PUT -d @$f $u/$d/$s')
%
%    license:
%        BSD license, see LICENSE_BSD.txt files for details
%
% -- this function is part of JBIDS toolbox (https://neurojson.org/#software)
%

if (nargin < 5)
    options = weboptions('');
end

if (~ischar(options) && ~isa(options, 'string'))
    options.ContentType = 'json';
    options.RequestMethod = 'POST';
    response = webwrite([couchdburl '/' dbname '/' docname], fileread(jsonfile), options);
else
    options = regexprep(options, '\$f', ['''' jsonfile '''']);
    options = regexprep(options, '\$u', couchdburl);
    options = regexprep(options, '\$d', dbname);
    options = regexprep(options, '\$s', docname);
    [status, response] = system(options);
    if (status ~= 0)
        error('command failed:\n%s\n', response);
    end
end
