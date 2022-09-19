function niml=afni_niml_read(fn)
% simple function to read niml files
%
% P=AFNI_NIML_READ(FN) reads a file FN and returns a NIML structure P.
%
% NNO Dec 2009 <n.oosterhof@bangor.ac.uk>

    bytes=read_bytes(fn);
    niml=afni_niml_parse(bytes);


function s=read_bytes(fn)
    fid=fopen(fn);

    if fid==-1
        if ~exist(fn,'file')
            error('File ''%s'' does not seem to exist', fn);
        else
            error('Error reading from file ''%s''', fn);
        end
    end

    file_closer=onCleanup(@()fclose(fid));

    s=fread(fid,inf,'*uint8')';
