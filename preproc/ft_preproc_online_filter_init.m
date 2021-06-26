function state = ft_preproc_online_filter_init(B, A, x)

% FT_PREPROC_ONLINE_FILTER_INIT initialize an IIR filter model with coefficients B
% and A, as used in filter and butter etc. The most recent sample of the signal must
% be given as a column vector.
%
% Use as
%   state = ft_preproc_online_filter_init(B, A, dat)
%
% This function will calculate the filter delay states such that the initial response
% will be as if the filter already have been applied forever.
%
% See also PREPROC, FT_PREPROC_ONLINE_FILTER_APPLY

% Copyright (C) 2010, Stefan Klanke
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

% Normalize filter coefficients if not already done so
A = A(:); % use column vector
B = B(:); % use column vector

if A(1)~=1
  B = B/A(1);
  A = A/A(1);
end

La = length(A);
Lb = length(B);

state = [];
state.d = size(x,1);

if La<Lb
  % pad A with zeros
  A = [A; zeros(Lb-La,1)];
  state.N = Lb-1;	% filter order
elseif La>Lb
  % pad B with zeros
  B = [B; zeros(La-Lb,1)];
  state.N = La-1;
else
  state.N = La-1;
end

state.A  = A;
state.B  = B;
state.A2 = A(2:end);
state.B1 = B(1);
state.B2 = B(2:end);

% this would be for direct form II, but MATLAB filter uses direct form II transpose
% state.z = x*(ones(1,state.N)/sum(B));

% there might be a faster way to compute this, but I can't think of any right now
% M is the matrix that describes the evolution of the filter in a 2+N-dim space
% composed of [output; delay states; input].
% We want to find the delay states corresponding to constant input (=1).
M = [[0; -A(2:end); 0],[eye(state.N);zeros(2,state.N)],[B;1]];
n = null(M-eye(2+state.N));  % = eigenvector of M corresponding to eigenvalue=1, that is n=M*n
n = n(:,1);               % ensure that it is only the first eigenvector
z = n(2:end-1);           % delay state part of it
z = z*(1-B(1))/z(1);      % scale appropiately
state.z = x*z';
