function [Prefix, Ext] = RemoveNIMLExtension (fname)
%
%   [Prefix, Ext] = RemoveNIMLExtension (fname)
%
%Purpose:
%   removes known .NIML filename extensions.
%   Look at code more more details.
%
%Input Parameters:
%   fname: filename
%
%
%Output Parameters:
%   Prefix: filename without extension
%   Ext: recognized extension of NIML filename, if any
%        if Ext is empty then no know extensions were found
%
%
%Key Terms:
%
%More Info :
%   RemoveExtension
%
%
%
%     Author : Ziad Saad
%     Date : Tue Nov 9 09:59:58 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'RemoveNIMLExtension';

%Debug Flag
DBG = 1;

Prefix = ''; Ext = '';

if (isempty(fname)) return; end
[Prefix, Ext] = RemoveExtension(fname,...
            '.niml|.niml.dset');

return;

