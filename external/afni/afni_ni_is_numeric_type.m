function isnum = afni_ni_is_numeric_type(t)
   ni_def=afni_ni_defs();

   if (isnumeric(t)),
      isnum = ( (t) >= 0 & (t) <= ni_def.NI_RGBA );
   elseif (isnumeric(afni_ni_rowtype_name_to_code(t))),
      isnum = ( (t) >= 0 & (t) <= ni_def.NI_RGBA );
   end

   isnum = prod(single(isnum));

   return;


