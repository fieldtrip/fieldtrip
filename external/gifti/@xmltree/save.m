function varargout = save(tree, filename)
% XMLTREE/SAVE Save an XML tree in an XML file
% FORMAT varargout = save(tree,filename)
%
% tree      - XMLTree
% filename  - XML output filename
% varargout - XML string
%__________________________________________________________________________
%
% Convert an XML tree into a well-formed XML string and write it into
% a file or return it as a string if no filename is provided.
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: save.m 4460 2011-09-05 14:52:16Z guillaume $


%error(nargchk(1,2,nargin));

prolog = '<?xml version="1.0" ?>\n';

%- Return the XML tree as a string
if nargin == 1
    varargout{1} = [sprintf(prolog) ...
        print_subtree(tree,'',root(tree))];
%- Output specified
else
    %- Filename provided
    if ischar(filename)
        [fid, msg] = fopen(filename,'w');
        if fid==-1, error(msg); end
        if isempty(tree.filename), tree.filename = filename; end
    %- File identifier provided
    elseif isnumeric(filename) && numel(filename) == 1
        fid = filename;
        prolog = ''; %- With this option, do not write any prolog
    else
        error('[XMLTree] Invalid argument.');
    end
    fprintf(fid,prolog);
    save_subtree(tree,fid,root(tree));
    if ischar(filename), fclose(fid); end
    if nargout == 1
        varargout{1} = print_subtree(tree,'',root(tree));
    end
end

%==========================================================================
function xmlstr = print_subtree(tree,xmlstr,uid,order)
    if nargin < 4, order = 0; end
    xmlstr = [xmlstr blanks(3*order)];
    switch tree.tree{uid}.type
        case 'element'
            xmlstr = sprintf('%s<%s',xmlstr,tree.tree{uid}.name);
            for i=1:length(tree.tree{uid}.attributes)
                xmlstr = sprintf('%s %s="%s"', xmlstr, ...
                    tree.tree{uid}.attributes{i}.key,...
                    tree.tree{uid}.attributes{i}.val);
            end
            if isempty(tree.tree{uid}.contents)
                xmlstr = sprintf('%s/>\n',xmlstr);
            else
                xmlstr = sprintf('%s>\n',xmlstr);
                for i=1:length(tree.tree{uid}.contents)
                    xmlstr = print_subtree(tree,xmlstr, ...
                        tree.tree{uid}.contents(i),order+1);
                end
                xmlstr = [xmlstr blanks(3*order)];
                xmlstr = sprintf('%s</%s>\n',xmlstr,...
                    tree.tree{uid}.name);
            end
        case 'chardata'
            xmlstr = sprintf('%s%s\n',xmlstr, ...
                entity(tree.tree{uid}.value));
        case 'cdata'
            xmlstr = sprintf('%s<![CDATA[%s]]>\n',xmlstr, ...
                tree.tree{uid}.value);
        case 'pi'
            xmlstr = sprintf('%s<?%s %s?>\n',xmlstr, ...
                tree.tree{uid}.target, tree.tree{uid}.value);
        case 'comment'
            xmlstr = sprintf('%s<!-- %s -->\n',xmlstr,...
                tree.tree{uid}.value);
        otherwise
            warning(sprintf('Type %s unknown: not saved', ...
                tree.tree{uid}.type));
    end

%==========================================================================
function save_subtree(tree,fid,uid,order)
    if nargin < 4, order = 0; end
    fprintf(fid,blanks(3*order));
    switch tree.tree{uid}.type
        case 'element'
            fprintf(fid,'<%s',tree.tree{uid}.name);
            for i=1:length(tree.tree{uid}.attributes)
                fprintf(fid,' %s="%s"',...
                tree.tree{uid}.attributes{i}.key, ...
                tree.tree{uid}.attributes{i}.val);
            end
            if isempty(tree.tree{uid}.contents)
                fprintf(fid,'/>\n');
            else
                fprintf(fid,'>\n');
                for i=1:length(tree.tree{uid}.contents)
                    save_subtree(tree,fid,...
                        tree.tree{uid}.contents(i),order+1)
                end
                fprintf(fid,blanks(3*order));
                fprintf(fid,'</%s>\n',tree.tree{uid}.name);
            end
        case 'chardata'
            fprintf(fid,'%s\n',entity(tree.tree{uid}.value));
        case 'cdata'
                fprintf(fid,'<![CDATA[%s]]>\n',tree.tree{uid}.value);
        case 'pi'
            fprintf(fid,'<?%s %s?>\n',tree.tree{uid}.target, ...
                tree.tree{uid}.value);
        case 'comment'
            fprintf(fid,'<!-- %s -->\n',tree.tree{uid}.value);
        otherwise
            warning(sprintf('[XMLTree] Type %s unknown: not saved', ...
                tree.tree{uid}.type));
    end


%==========================================================================
function str = entity(str)
    % This has the side effect of strtrim'ming the char array.
    str = char(strrep(cellstr(str), '&',  '&amp;' ));
    str = char(strrep(cellstr(str), '<',  '&lt;'  ));
    str = char(strrep(cellstr(str), '>',  '&gt;'  ));
    str = char(strrep(cellstr(str), '"',  '&quot;'));
    str = char(strrep(cellstr(str), '''', '&apos;'));
