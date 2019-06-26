function [ans, Prefix] = PrefixStatus (Prefix)
%
%   [Status, Prefix] = PrefixStatus (Name)
%
%Purpose:
%   Determines if the name chosen is a valid for an AFNI Brick
%   If you pass a name containing prefix only and no view then
%   View is assumed to be +orig 
%   
%   
%Input Parameters:
%   Name: a string containing an AFNI Prefix, 
%          or Prefix and view moon+orig or 
%           Prefix, view and extension moon+acpc.HEAD
%          
%   
%   
%Output Parameters:
%   Status : -1 error
%          : 0 Read Only
%          : 1 Read Write 
%   Prefix: The Prefix
%   
%      
%Key Terms:
%   
%More Info :
%   
%   The function checks for the existence of a brick with a similar
%   prefix 
%   
%
%     Author : Ziad Saad
%     Date : Wed Sep 19 10:30:05 EDT 2001
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'CheckPrefix';

%Debug Flag
DBG = 1;

%initailize return variables
ans = 0;

%make sure you have no HEAD or BRIK as an extension
[Prefix, Ext] = RemoveExtension(Prefix, '.HEAD|.BRIK');

%make sure you have no view
[Prefix, ViewName] = RemoveExtension(Prefix, '+orig|+acpc|+tlrc');

%check for conflict
if (~isempty(Ext) & isempty(ViewName))
	ans = -1; 
	ErrEval(FuncName,'Err_Expected a correct view (+orig, +acpc or +tlrc) but found none.');
	return;
end 

if ~isempty(ViewName), View = ViewName; else View = '+orig'; end

%look for various files
if (exist(sprintf('%s%s.HEAD', Prefix, View)) == 2), return; end
if (exist(sprintf('%s%s.BRIK', Prefix, View)) == 2), return; end




ans = 1;
return;

