function [err, Bstrct] = WaverBasisOption (BasisFunc, N_Basis, BasisOpt)
%
%   [err,] = WaverBasisOption ()
%
%Purpose:
%
%
%
%Input Parameters:
%
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%
%
%
%Key Terms:
%
%More Info :
%
%
%
%
%     Author : Ziad Saad
%     Date : Fri Jul 18 14:35:24 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'WaverBasisOption';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
Bstrct = [];

if (strcmp (lower(BasisFunc),'tent')),
   if (~isfield(BasisOpt,'tSpan') | isempty(BasisOpt.tSpan)),
      fprintf (2,'Error %s: tSpan must be set with tent option.\n',...
         FuncName);
      return;
   end

   for (i=1:1:N_Basis),
      Bstrct(i).BareOpt = sprintf ('tent((t-%g)/%g)',...
         (i-1).*BasisOpt.tSpan./N_Basis, BasisOpt.tSpan./N_Basis);
      Bstrct(i).opt = sprintf ('-EXPR ''%s''',...
         Bstrct(i).BareOpt);
      Bstrct(i).minlag = 0;
      Bstrct(i).maxlag = 0;
   end
   err = 0;
   return;
elseif (strcmp (lower(BasisFunc),'gam') | strcmp (lower(BasisFunc),'mgh')),
   if (N_Basis ~= 1),
     fprintf (2,'Error %s: Only 1 basis for GAM or MGH\n', FuncName);
     return;
   end
   Bstrct(1).BareOpt = '';
   Bstrct(1).opt = '-GAM';
	Bstrct(1).power = BasisOpt.gamb;
	Bstrct(1).scale = BasisOpt.gamc;
	Bstrct(1).delay = BasisOpt.gamd;
   if (~isfield(BasisOpt,'minlag') | isempty(BasisOpt.minlag)) BasisOpt.minlag = 0; end
   if (~isfield(BasisOpt,'maxlag') | isempty(BasisOpt.maxlag)) BasisOpt.maxlag = 0; end
   Bstrct(1).minlag = BasisOpt.minlag;
   Bstrct(1).maxlag = BasisOpt.maxlag;
   err = 0;
   return;
else
   if (strcmp (lower(BasisFunc),'spm1') | strcmp (lower(BasisFunc),'spm2')),
	   for (i=1:1:N_Basis),
         Bstrct(i).BareOpt = '';
         Bstrct(i).minlag = 0;
         Bstrct(i).maxlag = 0;
		end
		err = 0;
		return;
   else
   fprintf (2,'Error %s: Bad Basis %s\n', FuncName, BasisFunc);
   return;
	end
end

err = 0;
return;
