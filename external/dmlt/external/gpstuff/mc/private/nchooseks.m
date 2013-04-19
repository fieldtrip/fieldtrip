function c = nchooseks(N,K)
%NCHOOSEKS    Multinomial coefficient
%
%   C = NCHOOSEKS(N,K) where N is a non-negative integer and K is a vector
%   of k non-negative integers returns N!/(K(1)!...K(k))!)
%   This is the number of combinations of N things grouped in k groups with
%   sizes in vector K. SUM(K) must be be N.
%   When a coefficient is greater than 10^15, a warning will be produced
%   indicating possible inexact results. In such cases, the result is good 
%   to 15 digits.
%
%   See also NCHOOSEK

% Copyright (c) 2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if any(K < 0) || any(K ~= round(K)) || min(size(K)) ~=1
  error('The second input has to be a vector of non-negative integers'); 
end

[m, n] = size(N);
if m ~= 1 && n ~= 1
   error('The first argument has to be a non-negative integer');
end

if sum(K) ~= N, error('SUM(K) must be equal to N'); end 

nums=1:N;
ii=1;
for i=1:numel(K)
  for j=1:K(i)
    nums(ii)=ii/j;
    ii=ii+1;
  end
end
c=round(prod(nums));  

if c > 1e+015
  warning('MATLAB:nchoosek:LargeCoefficient',['Result may not be ' ...
                      'exact. Coefficient is greater than 10^15,\n         and ' ...
                      'is only good to 15 digits.']);
end
