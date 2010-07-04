function cod = afni_ni_rowtype_name_to_code(t)
   global nidef;
   
   if (isempty(nidef)) afni_ni_defs(); end
   nt = size(t,1);
   cod = -ones(1,nt);
   for (i=1:1:nt),
      mm = strmatch(t(i,:), nidef.type_alias, 'exact');
      if (isempty(mm)),
         mm = strmatch(t(i,:), nidef.type_name, 'exact');
      end
      if (~isempty(mm))   
         cod(i) = mm;
      end
   end
   return;

