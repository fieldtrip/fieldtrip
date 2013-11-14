function a = softmax2(n)
%SOFTMAX2     Softmax transfer function
%
%	Description
%	A = SOFTMAX2(N) takes a matrix N of network outputs and 
%       transfers it through a sotmax function to get A.
%       
%       Code:
%       temp = exp(n);
%       a = temp./(sum(temp, 2)*ones(1,size(n,2)));
%
%	See also
%	MLP2, MLP2PAK, MLP2UNPAK, MLP2ERR, MLP2BKP, MLP2GRAD
%

% Copyright (c) 1999 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

temp = exp(n);
a = temp./(sum(temp, 2)*ones(1,size(n,2)));
