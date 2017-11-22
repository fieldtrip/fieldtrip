function om_save_sym(data,filename,format)

% OM_SAVE_SYM   Save symmetric Matrix
%
%   Save symmetric Matrix
%
%   SYNTAX
%       OM_SAVE_SYM(DATA,FILENAME,FORMAT)
%
%       FORMAT : can be 'ascii' or 'binary' or 'matlab' (default)
%

% Copyright (C) 2010-2017, OpenMEEG developers

me = 'OM_SAVE_SYM';

if nargin == 0
    eval(['help ',lower(me)])
    return
end

if nargin < 3
    format = 'matlab';
end

dims = size(data);
assert(dims(1) == dims(2),'Matrix non square')
assert(isempty(find(data ~= data')),'Matrix non symmetric')

switch format
    case 'matlab'
        file = fopen(filename,'w');
        dim=length(data);
        data_raw=struct('symmatrix',struct('size',dim,'data',data(triu(ones(dim))>0)));
        save(filename,'-mat','-struct','data_raw','-v6')
        fclose(file);
        clear data_raw;
    case 'binary'
        disp(['Saving file ',filename])
        file = fopen(filename,'w');
        dim = dims(1);
        fwrite(file,dim,'uint32','ieee-le');
        fwrite(file,data(triu(ones(dim,dim)) > 0),'double','ieee-le');
        fclose(file);
    case 'ascii'
        for i=1:dims(1)
            if i == 1
                dlmwrite(filename, data(i,i:end), 'delimiter', '\t','precision','%.18e');
            else
                dlmwrite(filename, data(i,i:end), 'delimiter', '\t','-append','precision','%.18e');
            end
        end
    otherwise
        error([me,' : Unknown file format'])
end