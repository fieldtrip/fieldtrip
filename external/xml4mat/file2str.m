function y=file2str(x)

%FILE2STR reads textfile intoa single long string
%
% Syntax: y=file2str(x)
%
% Description
%   x is a filename
%   y is the long string with all contents
%
% Jonas Almeida, almeidaj@musc.edu, 30 June 2003, MAT4NAT Tbox


fid=fopen(x,'r');
i=1;
while ~feof(fid)
   y{i}=fgetl(fid);
   i=i+1;
end
fclose(fid);
y=strcat(y{:});