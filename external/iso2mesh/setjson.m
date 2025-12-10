function val = updatejson(fname, mmap, varargin)

if (isoctavemesh)
    if (regexp(fname, '^\s*(?:\[.*\])|(?:\{.*\})\s*$', 'once'))
        string = fname;
    elseif (exist(fname, 'file'))
        try
            encoding = jsonopt('Encoding', '', opt);
            if (isempty(encoding))
                string = fileread(fname);
            else
                fid = fopen(fname, 'r', 'n', encoding);
                string = fread(fid, '*char')';
                fclose(fid);
            end
        catch
            try
                string = urlread(['file://', fname]);
            catch
                string = urlread(['file://', fullfile(pwd, fname)]);
            end
        end
    else
        error_pos('input file does not exist');
    end
end

mmap = [mmap{:}];
keylist = mmap(1:2:end);

for i = 1:2:length(varargin)
    if (regexp(varargin{i}, '^$'))
        [tf, loc] = ismember(varargin{i}, keylist);
        if (tf)
            rec = {'uint8', [1, mmap{loc * 2}{1} - 1], 'o'; ...
                   'uint8', [1, mmap{loc * 2}{2}],  'x'};
            fmap = memmapfile(fname, 'writable', true, 'format', rec);
            val = savejson('', varargin{i + 1}, 'compact', 1);
            if (length(var) <= rec{2, 2}{2})
                val = [val repmat(' ', [1, rec{2, 2}{2} - length(var)])];
                fmap.x = val;
            end
        end
    end
end
