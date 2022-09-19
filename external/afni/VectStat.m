function [Stat] = VectStat (M)
%
%   [Stat] = VectStat (M)
%
%Purpose:
%   calculates very basic stats of vectors
%
%
%Input Parameters:
%   M : MxN matrix
%      NAN or Inf values are not considered and the user is warned of their existence
%
%Output Parameters:
%
%   Stat is an Mx1 vector of structures with the following fields
%     Each Stat(i) has the following structures
%      .N : the number of elements in M(i,:)
%      .M : the mean of vector M(i,:)
%      .Md: the median of vector M(i,:)
%      .S : the (unbiased) standard deviation of M(i,:)
%      .V : the variance (S^2) of M(i,:)
%      .mx: The maximum value
%      .mn: The minimum value
%
%     if the input vector is empty , NaN values are returned in Stat
%More Info :
%
%
%
%
%     Author : Ziad Saad
%     Date : Tue Apr 20 17:13:23 CDT 1999


%Define the function name for easy referencing
FuncName = 'VectStat';

%Debug Flag
DBG = 1;
Stat.N = NaN;
Stat.M = NaN;
Stat.Md = NaN;
Stat.S = NaN;
Stat.V = NaN;



if (is_row(M) ~= -1),
	sZ = 1;
	sz2 = length(M);
	M = M(:)'; %turn M into a row vector
else
	sZ = size(M,1);
	sz2 = size(M,2);
end

if (sZ == 0 | sz2 ==0),
	ErrEval(FuncName,'Wrn_Empty M returning NaN structure for result');
	return;
end

%initialize, once
	Stat(sZ).N = 0;Stat(sZ).M = 0;Stat(sZ).Md = 0;Stat(sZ).S = 0;Stat(sZ).V = 0; Stat(sZ).mn = 0; Stat(sZ).mx = 0;

for (i=1:1:sZ),
	igood = find(isfinite(M(i,:)));
	N_igood = length(igood);
	if (N_igood < sz2),
		stmp = sprintf ('%s Warning: Inf or Nan ignored in row %g\n', FuncName, i);
	end
	
	Stat(i).N = N_igood;
	Stat(i).M = mean(M(i,igood));
	Stat(i).Md = median(M(i,igood));
	Stat(i).S = std(M(i,igood),0); %unbiased estimator
	Stat(i).V = Stat(i).S .^ 2;
	Stat(sZ).mn = min(M(i,igood));
	Stat(sZ).mx = max(M(i,igood));
end


return;

