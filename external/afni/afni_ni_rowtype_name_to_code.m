function cod = afni_ni_rowtype_name_to_code(t)
   ni_def=afni_ni_defs();

   nt = size(t,1);
   cod = -ones(1,nt);
   for (i=1:1:nt),
      mm = strmatch(t(i,:), ni_def.type_alias, 'exact');
      if (isempty(mm)),
         mm = strmatch(t(i,:), ni_def.type_name, 'exact');
      end
      if (~isempty(mm))
         cod(i) = mm;
      end
   end
   return;

