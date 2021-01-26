function pretty_print_struct(st, indent, option, logger)
spaces = '';

if ~exist('st','var') || isempty(st)
    return;
end
if ~exist('indent','var') || isempty(indent)
    indent = 0;
end
if ~exist('option','var') || isempty(option)
    option = 1;
end
if ~exist('logger','var')
    logger = [];
end
logger = InitLogger(logger);

if iswholenum(indent)
    spaces = blanks(indent);
elseif ischar(indent)
    spaces = blanks(length(indent));
end

if isstruct(st) || isobject(st)
    s = evalc('disp(st)');
    c = str2cell_fast(s, char(10));
    for ii=1:length(c)
        if option==1
            logger.Write(sprintf('%s%s\n', spaces, strtrim_improve(c{ii})));
        elseif option==2
            logger.Write(sprintf('%s%s\n', spaces, c{ii}));
        end
    end
else
    str = '';
    for jj=1:ndims(st)
        if jj==1
            str = num2str(size(st,jj));
        else
            str = sprintf('%sx%s', str, num2str(size(st,jj)));
        end
    end    
    logger.Write(sprintf('        %s: [%s %s]\n', inputname(1), str, class(st)));
end

