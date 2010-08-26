function FM = online_filter_init(B, A, x)
% function FM = online_filter_init(B, A, x)
%
% Initialize an IIR filter model with coefficients B and A, as used in filter and butter etc.
% One sample x of the signal must be given as a column vector.
%
% This function will calculate the filter delay states such that the initial response 
% is as if 'x' would have been applied since forever.

% 2010 S. Klanke

% Normalize filter coefficients if not already done so
A = A(:); % use column vector
B = B(:); % use column vector

if A(1)~=1
	B = B/A(1);
	A = A/A(1);
end

La = length(A);
Lb = length(B);

FM = [];
FM.d = size(x,1);

if La<Lb
   % pad A with zeros
   A = [A; zeros(Lb-La,1)];
   FM.N = Lb-1;	% filter order
elseif La>Lb
   % pad B with zeros
   B = [B; zeros(La-Lb,1)];
   FM.N = La-1;
else
   FM.N = La-1;
end

FM.A  = A;
FM.B  = B;
FM.A2 = A(2:end);
FM.B1 = B(1);
FM.B2 = B(2:end);
FM.z = x*(ones(1,FM.N)/sum(B));