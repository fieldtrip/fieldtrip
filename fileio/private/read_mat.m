function [matrix, extra] = read_mat(fn)

% READ_MAT reads a matrix from an ascii or binary MBF format file
%
%  Usage: m         = loadmat('file');
%     or  [m,extra] = loadmat('file');
%  
%  LOADMAT('file') returns the matrix stored in 'file' and
%  the extra information stored at the bottom of that file.
%  LOADMAT works for binary as well as asci matrix files.
%  
%  See also WRITE_MAT

% Copyright (C) 1998, Thom Oostendorp 
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

f=fopen(fn);
if (f==-1)
  fprintf('\nCannot open %s\n\n', fn);
  result=0;
  extra='';
  return;
end

[N,nr]=fscanf(f,'%d',2);
if (nr~=2)
  fclose(f);
  f=fopen(fn);
  [magic ,nr]=fread(f,8,'char');
  if (char(magic')==';;mbfmat')
    fread(f,1,'char');
    hs=fread(f,1,'long');
    fread(f,1,'char');
    fread(f,1,'char');
    fread(f,1,'char');
    N=fread(f,2,'long');
    M=fread(f,[N(2),N(1)],'double');
  else
    fclose(f);
    f=fopen(fn);
    N=fread(f,2,'long');
    M=fread(f,[N(2),N(1)],'float');
  end
else
  M=fscanf(f,'%f',[N(2) N(1)]);
end

[extra,nextra]=fread(f,1000,'char');
fclose(f);
S=sprintf('\n%s contains %d rows and %d columns\n', fn, N(1), N(2));
disp(S);
if (nextra~=0)
  S=sprintf('%s contains the following extra information:\n', fn);
  disp(S);
  disp(char(extra'));
end

matrix=M';
extra=char(extra');

