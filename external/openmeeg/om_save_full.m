function om_save_full(data,filename,format)

% OM_SAVE_FULL   Save full Matrix
%
%   Save full Matrix
%
%   SYNTAX
%       OM_SAVE_FULL(DATA,FILENAME,FORMAT)
%
%       FORMAT : can be 'ascii' or 'binary' or 'matlab' (default)
%

% Copyright (C) 2010-2017, OpenMEEG developers

me = 'OM_SAVE_FULL';

if nargin == 0
    eval(['help ',lower(me)])
    return
end

if nargin < 3
    format = 'matlab';
end

switch format
    case 'matlab'
        file = fopen(filename,'w');
        data_raw=struct('linop',data);
        save(filename,'-mat','-struct','data_raw','-v7')
        fclose(file);
        clear data_raw;
    case 'binary'
        disp(['Saving file ',filename])
        file = fopen(filename,'w');
        dims = size(data);
        fwrite(file,dims,'uint32','ieee-le');
        fwrite(file,data(:),'double','ieee-le');
        fclose(file);
    case 'ascii'
        data = double(data);
        save(filename,'data','-ascii','-double')
    otherwise
        error([me,' : Unknown file format'])
end

