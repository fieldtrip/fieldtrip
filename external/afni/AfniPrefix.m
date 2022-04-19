function [Prefix, ViewName, Ext] = AfniPrefix (Prefix)
%
%   [Prefix, View, Ext] = AfniPrefix (Name)
%
%Purpose:
%   returns the afni prefix, view and extension from a brick name
%   No checking is done to make sure that View or Ext are acceptable
%   see PrefixStatus instead
%
%Input Parameters:
%   Name : 'elvis' or 'elvis+orig' or 'elvis+orig.HEAD'
%
%
%Output Parameters:
%   Prefix: 'elvis'
%   View: '+orig'
%   Ext: .HEAD
%
%Key Terms:
%
%More Info :
%   PrefixStatus
%
%
%
%     Author : Ziad Saad
%     Date : Tue Dec 4 15:24:09 EST 2001
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'AfniPrefix';

[Prefix, Ext] = RemoveExtension(Prefix, '.HEAD|.BRIK');
[Prefix, ViewName] = RemoveExtension(Prefix, '+orig|+acpc|+tlrc');






err = 0;
return;

