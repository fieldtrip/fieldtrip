function [flg] = filexist (filename)
%
%	[flg] = filexist (filename)
%
%  if filename exists then flg = 1, else flg = 0
%

flg=0;

if exist(filename,'file')
    flg=1;
end
