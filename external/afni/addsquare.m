function [err, hndl] = addsquare(Pll, Pur, C)
%
% [err, hndl] = addsquare(Pll, Pur, C)
%
%
%Purpose:
%  adds a rectangle on current plot
%   
%   
%   
%Input Parameters:
%  Pll : [x y] of lower left point
%  Pur : [x y] of upper right point
%  C  : rgb triplet, between 0 and 1  rectangle color
%   
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   hndl : the handle to the rectangle object
%   
%      
%Key Terms:
%   
%More Info :
%   
%   add the square first, before making the plots,
%   otherwise, the patch will obscure the plots
%   
%
%     Author : Ziad Saad
%     Date : Fri Oct 06 15:47:14 EDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'addsquare';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

X = [Pll(1) Pur(1) Pur(1) Pll(1)]';
Y = [Pll(2) Pll(2) Pur(2) Pur(2)]'; 

hold on;
hndl = patch(X,Y, C);
set(hndl, 'EdgeColor', C);




err = 0;
return;

