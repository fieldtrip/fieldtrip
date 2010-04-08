function [xapp,xtest] = normalizemeanstd(xapp,xtest)

% NORMALIZEMEANSTD normalizes data xapp to mean=0 and standard deviation=1
% (z-scores) and applies the resulting transformation to xtest
% 
% Use as
%           [xapp,xtest] = normalizemeanstd(xapp,xtest) OR
%           [xapp] = normalizemeanstd(xapp) 
%
% INPUT/OUTPUT
%           xapp  - training (primary) feature data (observations in rows and variables in columns)
%           xtest - test (secondary) feature data (observations in rows and variables in columns)
% 

meanxapp=mean(xapp);
stdxapp=std(xapp);
nbxapp=size(xapp,1);

nbvar=size(xapp,2);

xapp= (xapp - ones(nbxapp,1)*meanxapp)./ (ones(nbxapp,1)*stdxapp) ;

if nargin >1
    nbxtest=size(xtest,1);
    xtest= (xtest - ones(nbxtest,1)*meanxapp)./ (ones(nbxtest,1)*stdxapp );

end;