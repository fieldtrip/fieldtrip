function afni_niml_write(p,fn,format)
% simple function to write niml files in ASCII format
%
% AFNI_NIML_WRITE(P,FN) writes a NIML struct to file FN. P should be in the
% form as returned by AFNI_NIML_READ.
%
% The order of the arguments may be reversed. If FN is numeric, then it is
% assumed to be a file identifier. If P is omitted, output is written to
% stdout.
%
% This function requires a tree-like structure similar to the NIML format.
% To write a 'simple' struct, use AFNI_NIML_WRITESIMPLE.
%
% NNO Dec 2009 <n.oosterhof@bangor.ac.uk>

if nargin<3
    format='binary';
end

if nargin<2
    fn=1;
end

if isstruct(fn) || iscell(fn) % allow for wrong order of arguments; swap
    tmp=p;
    p=fn;
    fn=tmp;
end

% conversion to string
s=afni_niml_print(p, format);

% see if output is a file identifier
if isnumeric(fn) && round(fn)==fn
    fid=fn;
    fn=sprintf('FID %d',fid);
else
    fid=fopen(fn,'w');
end

if fid==0
    error('Could not write to %s\n', fn);
end

fwrite(fid,s);

if fid>2 % we don't close standard input, output, or error
    fclose(fid);
end
