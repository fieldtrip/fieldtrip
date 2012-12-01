function I=processvarargin(I,D)
%function I2=processvarargin(I,D)
%
%PURPOSE
%
%Subfunction of various functions in Icasso.
%
%INPUT
%
% I (cell array) input {'id1' value1 'id2' value2 ...}
% D (cell array) defaults {'id1 value1 'id2' value2 'id3' value3 ...}
%
%OUTPUT
%
% I2 (cell array) combination of inputs I and defaults D.
%
%DETAILS
%
%Checks that the format of input I is correct: 
%Error if the input array I is not of form 
% {string whatever string whatever ...}. 
%That is, there are 2N elements in the array and elements with odd
%index are strings (identifiers). 
%
%Adds defaults, i.e., adds those id, value pairs in 
%cell array D whose id do _not_ appear in input I.

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.2 johan 100105

I=I(:); 
D=D(:);
Ninput=length(I);

%%% Check number of pairs
if fix(Ninput/2)~=Ninput/2,
  error('Optional ''identifier'',value pairs mismatch!');
end

%%% Take default ids and values from the cell array

for i=1:2:Ninput,
  if ~ischar(I{i})
    error('Optional ''identifier'',value pairs mismatch!');
  end
  %% Check if some default is redefined, and drop the default away
  
  index2default=find(strcmp(lower(D(1:2:end)),lower(I{i})));
  if any(index2default);
    index2default=index2default(:)*2-1;
    index2default=[index2default;index2default+1];
    D(index2default)=[];
  end
end

I=[I;D]';
  
  
  