function tree = xml_parser(xmlstr)
% XML (eXtensible Markup Language) Processor
% FORMAT tree = xml_parser(xmlstr)
%
% xmlstr  - XML string to parse
% tree    - tree structure corresponding to the XML file
%__________________________________________________________________________
%
% xml_parser.m is an XML 1.0 (http://www.w3.org/TR/REC-xml) parser.
% It aims to be fully conforming. It is currently not a validating 
% XML processor.
%
% A description of the tree structure provided in output is detailed in 
% the header of this m-file.
%__________________________________________________________________________
% Copyright (C) 2002-2015  http://www.artefact.tk/

% Guillaume Flandin
% $Id$

% XML Processor for GNU Octave and MATLAB (The Mathworks, Inc.)
% Copyright (C) 2002-2015 Guillaume Flandin <Guillaume@artefact.tk>
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation Inc, 59 Temple Pl. - Suite 330, Boston, MA 02111-1307, USA.
%--------------------------------------------------------------------------

% Suggestions for improvement and fixes are always welcome, although no
% guarantee is made whether and when they will be implemented.
% Send requests to <Guillaume@artefact.tk>
% Check also the latest developments on the following webpage:
%           <http://www.artefact.tk/software/matlab/xml/>
%--------------------------------------------------------------------------

% The implementation of this XML parser is much inspired from a 
% Javascript parser that used to be available at <http://www.jeremie.com/>

% A C-MEX file xml_findstr.c is also required, to encompass some
% limitations of the built-in FINDSTR function.
% Compile it on your architecture using 'mex -O xml_findstr.c' command
% if the compiled version for your system is not provided.
% If this function does not behave as expected, comment the line
% '#define __HACK_MXCHAR__' in xml_findstr.c and compile it again.
%--------------------------------------------------------------------------

% Structure of the output tree:
% There are 5 types of nodes in an XML file: element, chardata, cdata,
% pi and comment.
% Each of them contains an UID (Unique Identifier): an integer between
% 1 and the number of nodes of the XML file.
%
%    element (a tag <name key="value"> [contents] </name>
%       |_ type:       'element'
%       |_ name:       string
%       |_ attributes: cell array of struct 'key' and 'value' or []
%       |_ contents:   double array of uid's or [] if empty
%       |_ parent:     uid of the parent ([] if root)
%       |_ uid:        double
%
%    chardata (a character array)
%       |_ type:   'chardata'
%       |_ value:  string
%       |_ parent: uid of the parent
%       |_ uid:    double
%
%    cdata (a litteral string <![CDATA[value]]>)
%       |_ type:   'cdata'
%       |_ value:  string
%       |_ parent: uid of the parent
%       |_ uid:    double
%
%      pi (a processing instruction <?target value ?>)
%       |_ type:   'pi' 
%       |_ target: string (may be empty)
%       |_ value:  string
%       |_ parent: uid of the parent
%       |_ uid:    double
%
%    comment (a comment <!-- value -->)
%       |_ type:   'comment'
%       |_ value:  string
%       |_ parent: uid of the parent
%       |_ uid:    double
%
%--------------------------------------------------------------------------

% TODO/BUG/FEATURES:
%  - [compile] only a warning if TagStart is empty ?
%  - [attribution] should look for " and ' rather than only "
%  - [main] with normalize as a preprocessing, CDATA are modified
%  - [prolog] look for a DOCTYPE in the whole string even if it occurs
%    only in a far CDATA tag, bug even if the doctype is inside a comment
%  - [tag_element] erode should replace normalize here
%  - remove globals? uppercase globals  rather persistent (clear mfile)?
%  - xml_findstr is indeed xml_strfind according to Mathworks vocabulary
%  - problem with entities: do we need to convert them here? (&eacute;)
%--------------------------------------------------------------------------

%- XML string to parse and number of tags read
global xmlstring Xparse_count xtree;

%- Check input arguments
%error(nargchk(1,1,nargin));
if isempty(xmlstr)
    error('[XML] Not enough parameters.')
elseif ~ischar(xmlstr) || sum(size(xmlstr)>1)>1
    error('[XML] Input must be a string.')
end

%- Initialize number of tags (<=> uid)
Xparse_count = 0;

%- Remove prolog and white space characters from the XML string
xmlstring = normalize(prolog(xmlstr));

%- Initialize the XML tree
xtree = {};
tree = fragment;
tree.str = 1;
tree.parent = 0;

%- Parse the XML string
tree = compile(tree);

%- Return the XML tree
tree = xtree;

%- Remove global variables from the workspace
clear global xmlstring Xparse_count xtree;

%==========================================================================
% SUBFUNCTIONS

%--------------------------------------------------------------------------
function frag = compile(frag)
    global xmlstring xtree Xparse_count;
    
    while 1,
        if length(xmlstring)<=frag.str || ...
           (frag.str == length(xmlstring)-1 && strcmp(xmlstring(frag.str:end),' '))
            return
        end
        TagStart = xml_findstr(xmlstring,'<',frag.str,1);
        if isempty(TagStart)
            %- Character data
            error('[XML] Unknown data at the end of the XML file.');
            Xparse_count = Xparse_count + 1;
            xtree{Xparse_count} = chardata;
            xtree{Xparse_count}.value = erode(entity(xmlstring(frag.str:end)));
            xtree{Xparse_count}.parent = frag.parent;
            xtree{frag.parent}.contents = [xtree{frag.parent}.contents Xparse_count];
            frag.str = '';
        elseif TagStart > frag.str
            if strcmp(xmlstring(frag.str:TagStart-1),' ')
                %- A single white space before a tag (ignore)
                frag.str = TagStart;
            else
                %- Character data
                Xparse_count = Xparse_count + 1;
                xtree{Xparse_count} = chardata;
                xtree{Xparse_count}.value = erode(entity(xmlstring(frag.str:TagStart-1)));
                xtree{Xparse_count}.parent = frag.parent;
                xtree{frag.parent}.contents = [xtree{frag.parent}.contents Xparse_count];
                frag.str = TagStart;
            end
        else 
            if strcmp(xmlstring(frag.str+1),'?')
                %- Processing instruction
                frag = tag_pi(frag);
            else
                if length(xmlstring)-frag.str>4 && strcmp(xmlstring(frag.str+1:frag.str+3),'!--')
                    %- Comment
                    frag = tag_comment(frag);
                else
                    if length(xmlstring)-frag.str>9 && strcmp(xmlstring(frag.str+1:frag.str+8),'![CDATA[')
                        %- Litteral data
                        frag = tag_cdata(frag);
                    else
                        %- A tag element (empty (<.../>) or not)
                        if ~isempty(frag.end)
                            endmk = ['/' frag.end '>'];
                        else 
                            endmk = '/>';
                        end
                        if strcmp(xmlstring(frag.str+1:frag.str+length(frag.end)+2),endmk) || ...
                            strcmp(strip(xmlstring(frag.str+1:frag.str+length(frag.end)+2)),endmk)
                            frag.str = frag.str + length(frag.end)+3;
                            return
                        else
                            frag = tag_element(frag);
                        end
                    end
                end
            end
        end
    end

%--------------------------------------------------------------------------
function frag = tag_element(frag)
    global xmlstring xtree Xparse_count;
    close =  xml_findstr(xmlstring,'>',frag.str,1);
    if isempty(close)
        error('[XML] Tag < opened but not closed.');
    else
        empty = strcmp(xmlstring(close-1:close),'/>');
        if empty
            close = close - 1;
        end
        starttag = normalize(xmlstring(frag.str+1:close-1));
        nextspace = xml_findstr(starttag,' ',1,1);
        attribs = '';
        if isempty(nextspace)
            name = starttag;
        else
            name = starttag(1:nextspace-1);
            attribs = starttag(nextspace+1:end);
        end
        Xparse_count = Xparse_count + 1;
        xtree{Xparse_count} = element;
        xtree{Xparse_count}.name = strip(name);
        if frag.parent
            xtree{Xparse_count}.parent = frag.parent;
            xtree{frag.parent}.contents = [xtree{frag.parent}.contents Xparse_count];
        end
        if ~isempty(attribs)
            xtree{Xparse_count}.attributes = attribution(attribs);
        end
        if ~empty
            contents = fragment;
            contents.str = close+1;
            contents.end = name;
            contents.parent = Xparse_count;
            contents = compile(contents);
            frag.str = contents.str;
        else
            frag.str = close+2;
        end
    end

%--------------------------------------------------------------------------
function frag = tag_pi(frag)
    global xmlstring xtree Xparse_count;
    close = xml_findstr(xmlstring,'?>',frag.str,1);
    if isempty(close)
        warning('[XML] Tag <? opened but not closed.')
    else
        nextspace = xml_findstr(xmlstring,' ',frag.str,1);
        Xparse_count = Xparse_count + 1;
        xtree{Xparse_count} = pri;
        if nextspace > close || nextspace == frag.str+2
            xtree{Xparse_count}.value = erode(xmlstring(frag.str+2:close-1));
        else
            xtree{Xparse_count}.value = erode(xmlstring(nextspace+1:close-1));
            xtree{Xparse_count}.target = erode(xmlstring(frag.str+2:nextspace));
        end
        if frag.parent
            xtree{frag.parent}.contents = [xtree{frag.parent}.contents Xparse_count];
            xtree{Xparse_count}.parent = frag.parent;
        end
        frag.str = close+2;
    end

%--------------------------------------------------------------------------
function frag = tag_comment(frag)
    global xmlstring xtree Xparse_count;
    close = xml_findstr(xmlstring,'-->',frag.str,1);
    if isempty(close)
        warning('[XML] Tag <!-- opened but not closed.')
    else
        Xparse_count = Xparse_count + 1;
        xtree{Xparse_count} = comment;
        xtree{Xparse_count}.value = erode(xmlstring(frag.str+4:close-1));
        if frag.parent
            xtree{frag.parent}.contents = [xtree{frag.parent}.contents Xparse_count];
            xtree{Xparse_count}.parent = frag.parent;
        end
        frag.str = close+3;
    end

%--------------------------------------------------------------------------
function frag = tag_cdata(frag)
    global xmlstring xtree Xparse_count;
    close = xml_findstr(xmlstring,']]>',frag.str,1);
    if isempty(close)
        warning('[XML] Tag <![CDATA[ opened but not closed.')
    else
        Xparse_count = Xparse_count + 1;
        xtree{Xparse_count} = cdata;
        xtree{Xparse_count}.value = xmlstring(frag.str+9:close-1);
        if frag.parent
            xtree{frag.parent}.contents = [xtree{frag.parent}.contents Xparse_count];
            xtree{Xparse_count}.parent = frag.parent;
        end
        frag.str = close+3;
    end

%--------------------------------------------------------------------------
function all = attribution(str)
    %- Initialize attributs
    nbattr = 0;
    all = cell(nbattr);
    %- Look for 'key="value"' substrings
    while 1,
        eq = xml_findstr(str,'=',1,1);
        if isempty(str) || isempty(eq), return; end
        id = sort([xml_findstr(str,'"',1,1),xml_findstr(str,'''',1,1)]); id=id(1);
        nextid = sort([xml_findstr(str,'"',id+1,1),xml_findstr(str,'''',id+1,1)]);nextid=nextid(1);
        nbattr = nbattr + 1;
        all{nbattr}.key = strip(str(1:(eq-1)));
        all{nbattr}.val = entity(str((id+1):(nextid-1)));
        str = str((nextid+1):end);
    end

%--------------------------------------------------------------------------
function elm = element
    global Xparse_count;
    elm = struct('type','element','name','','attributes',[],'contents',[],'parent',[],'uid',Xparse_count);
   
%--------------------------------------------------------------------------
function cdat = chardata
    global Xparse_count;
    cdat = struct('type','chardata','value','','parent',[],'uid',Xparse_count);
   
%--------------------------------------------------------------------------
function cdat = cdata
    global Xparse_count;
    cdat = struct('type','cdata','value','','parent',[],'uid',Xparse_count);
   
%--------------------------------------------------------------------------
function proce = pri
    global Xparse_count;
    proce = struct('type','pi','value','','target','','parent',[],'uid',Xparse_count);

%--------------------------------------------------------------------------
function commt = comment
    global Xparse_count;
    commt = struct('type','comment','value','','parent',[],'uid',Xparse_count);

%--------------------------------------------------------------------------
function frg = fragment
    frg = struct('str','','parent','','end','');

%--------------------------------------------------------------------------
function str = prolog(str)
    %- Initialize beginning index of elements tree
    b = 1;
    %- Initial tag
    start = xml_findstr(str,'<',1,1);
    if isempty(start) 
        error('[XML] No tag found.')
    end
    %- Header (<?xml version="1.0" ... ?>)
    if strcmpi(str(start:start+2),'<?x')
        close = xml_findstr(str,'?>',1,1);
        if ~isempty(close) 
            b = close + 2;
        else 
            warning('[XML] Header tag incomplete.')
        end
    end
    %- Doctype (<!DOCTYPE type ... [ declarations ]>)
    start = xml_findstr(str,'<!DOCTYPE',b,1);  % length('<!DOCTYPE') = 9
    if ~isempty(start) 
        close = xml_findstr(str,'>',start+9,1);
        if ~isempty(close)
            b = close + 1;
            dp = xml_findstr(str,'[',start+9,1);
            if (~isempty(dp) && dp < b)
                k = xml_findstr(str,']>',start+9,1);
                if ~isempty(k)
                    b = k + 2;
                else
                    warning('[XML] Tag [ in DOCTYPE opened but not closed.')
                end
            end
        else
            warning('[XML] Tag DOCTYPE opened but not closed.')
        end
    end
    %- Skip prolog from the xml string
    str = str(b:end);

%--------------------------------------------------------------------------
function str = strip(str)
    str(isspace(str)) = '';

%--------------------------------------------------------------------------
function str = normalize(str)
    % Find white characters (space, newline, carriage return, tabs, ...)
    i = isspace(str);
    i = find(i == 1);
    str(i) = ' ';
    % replace several white characters by only one
    if ~isempty(i)
        j = i - [i(2:end) i(end)];
        str(i(j == -1)) = [];
    end

%--------------------------------------------------------------------------
function str = entity(str)
    str = strrep(str,'&lt;','<');
    str = strrep(str,'&gt;','>');
    str = strrep(str,'&quot;','"');
    str = strrep(str,'&apos;','''');
    str = strrep(str,'&amp;','&');
   
%--------------------------------------------------------------------------
function str = erode(str)
    if ~isempty(str) && str(1)==' ', str(1)=''; end;
    if ~isempty(str) && str(end)==' ', str(end)=''; end;
