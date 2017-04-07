function tree = xmltree(varargin)
% XMLTREE/XMLTREE Constructor of the XMLTree class
% FORMAT tree = xmltree(varargin)
% 
% varargin - XML filename or XML string
% tree     - XMLTree Object
%
%     tree = xmltree;             % creates a minimal XML tree: '<tag/>'
%     tree = xmltree('foo.xml');  % creates a tree from XML file 'foo.xml'
%     tree = xmltree('<tag>content</tag>') % creates a tree from string
%__________________________________________________________________________
%
% This is the constructor of the XMLTree class. 
% It creates a tree of an XML 1.0 file (after parsing) that is stored 
% using a Document Object Model (DOM) representation.
% See http://www.w3.org/TR/REC-xml for details about XML 1.0.
% See http://www.w3.org/DOM/ for details about DOM platform.
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: xmltree.m 4460 2011-09-05 14:52:16Z guillaume $

switch(nargin)
    case 0
        tree.tree{1} = struct('type','element',...
                              'name','tag',...
                              'attributes',[],...
                              'contents',[],...
                              'parent',[],...
                              'uid',1);
        tree.filename = '';
        tree = class(tree,'xmltree');
    case 1
        if isa(varargin{1},'xmltree')
            tree = varargin{1};
        elseif ischar(varargin{1})
            % Input argument is an XML string
            if (~exist(varargin{1},'file') && ...
                ~isempty(xml_findstr(varargin{1},'<',1,1)))
                tree.tree = xml_parser(varargin{1});
                tree.filename = '';
            % Input argument is an XML filename
            else
                fid = fopen(varargin{1},'rt');
                if (fid == -1) 
                    error(['[XMLTree] Cannot open ' varargin{1}]);
                end
                xmlstr = fread(fid,'*char')';
                %xmlstr = fscanf(fid,'%c');
                fclose(fid);
                tree.tree = xml_parser(xmlstr);
                tree.filename = varargin{1};
            end
            tree = class(tree,'xmltree');
        else 
            error('[XMLTree] Bad input argument');
        end
    otherwise
        error('[XMLTree] Too many input arguments');
end
