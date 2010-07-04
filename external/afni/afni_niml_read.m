function [p,s]=afni_niml_read(fn)
% simple function to read niml files
%
% P=AFNI_NIML_READ(FN) reads a file FN and returns a NIML structure P.
% [P,S]=AFNI_NIML_READ(FN) also returns the file contents S
%
% If FN refers to a file in another format than ASCII NIML, an attempt is
% made to convert the file to this format with AFNI's ConvertDset
%
% NNO Dec 2009 <n.oosterhof@bangor.ac.uk>

s=simpleread(fn);

[b1, ext]=is_non_niml_ascii(fn); % if the extension suggests that it is not ASCII NIML
if b1 || is_non_niml_ascii(s);   % ... or the file contents suggest this
                                 % then we try to convert the file   
    warning('Input file seems to be not a NIML ASCII file, will try to convert it to ascii');
    s=binaryToASCII(s,ext);
end
p=afni_niml_parse(s);


function s=simpleread(fn)
fprintf('Reading %s\n', fn);
fid=fopen(fn);
if fid==-1
    error('Error reading from file %s\n', fn);
end
s=fread(fid,inf,'char=>char');
fclose(fid);
s=s';


function sa=binaryToASCII(sb,ext)
% converts a non-ASCII NIML file to an ASCII NIML file
% ext is the extension of the input file

% keep looking for a temporary file name that does not exist yet
rand('twister',sum(100*clock)); % initialize random number generator
while true
    idx=floor(rand()*1e6);
    prefix=sprintf('__TMP_%d',idx);
    tmpfn1=sprintf('%s_in%s',prefix,ext);          % input 
    tmpfn2=sprintf('%s_out%s',prefix,'.niml.dset'); % output
    if ~exist(tmpfn1,'file') && ~exist(tmpfn2,'file')
        break;
    end
end

fid=fopen(tmpfn1,'w');
if fid<0 
    error('Could not open temporary file for writing');
end

fwrite(fid,sb);
fclose(fid);

cmd=sprintf('ConvertDset -input %s -o_niml_asc -prefix %s', tmpfn1, tmpfn2);
fprintf('Running command: %s\n', cmd);
[a,w]=unix(cmd);

if a ~= 0
    fprintf(w)
    error('Error when trying to convert binary to ascii');
end

fprintf('Conversion was successful\n');

% read the output file
sa=simpleread(tmpfn2);

% clean up
delete(tmpfn1);
delete(tmpfn2);


function [b,ext]=is_non_niml_ascii(s)
% tries to determine whether the file is niml ascii or not
% - input argument can either be a filename or the contents of the file
% - ouput contains whether it's a NIML ASCII file, and the extension
[ffff,fn,ext1]=fileparts(s);
if exist(s,'file')
    % s is a filename
    exts={'.gii','1D'};
    b=~isempty(strmatch(lower(ext1),exts));
    [ffff,fffff,ext2]=fileparts(fn);
    ext=[ext2 ext1]; % in case there is a double extension such as .niml.dset
else
    % s contains the file contents
    expr='<(?<h>\w+).*ni_form\s*=\s*"(?<rhs>[^"]+)".*</\1>';
    hh=regexp(s,expr,'names');
    b=false;
    for k=1:numel(hh)
        b=~isempty(findstr(hh(k).rhs,'binary.'));
        if b
            return
        end
    end

    ext=NaN; %unknown
end



    
