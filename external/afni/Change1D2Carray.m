function Change1D2Carray(f1d, decplaces, fc)
% take a 1D file and write it out as a C arrary
%  f1d : 1D file name
%  decplaces: -1 use all precision
%             0  ints
%             > 0 number of decimal places
[err,v] = Read_1D(f1d);

if (decplaces >= 0),
   v = round(v .* 10.^decplaces)/10.^decplaces;
end

fcid = fopen(fc, 'a');

varname = zdeblank(f1d);
ip = find (varname == '.');
if (~isempty(ip)), varname(ip) = '_'; end

fprintf(fcid, '/* Translation from %s using function Change1D2Carray.m */\n', f1d);
fprintf(fcid, 'int d1_%s = %d;\n', varname , size(v,1));
fprintf(fcid, 'int d2_%s = %d;\n', varname , size(v,2));
if (decplaces ~= 0),  fprintf(fcid, 'float %s[%d][%d] = {\n', varname , size(v,1), size(v,2));
else fprintf(fcid, 'int %s[%d][%d] = {\n', varname, size(v,1), size(v,2));
end
for(i=1:1:size(v,1)),
   fprintf(fcid, '   { ');
   for(j=1:1:size(v,2)),
      if (j<size(v,2)) fprintf(fcid, '%g, ', v(i,j));
      else fprintf(fcid, '%g ', v(i,j));
      end
   end
   if (i<size(v,1)) fprintf(fcid, '},\n');
   else fprintf(fcid, '}\n');
   end
end
fprintf(fcid, '};\n');
fclose(fcid);
