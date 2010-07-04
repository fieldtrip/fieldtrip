function isnum = afni_ni_is_numeric_type(t) 
   global nidef;
   
   if (isempty(nidef)) afni_ni_defs(); end
   
   if (isnumeric(t)),
      isnum = ( (t) >= 0 & (t) <= nidef.NI_RGBA );
   elseif (isnumeric(afni_ni_rowtype_name_to_code(t))),
      isnum = ( (t) >= 0 & (t) <= nidef.NI_RGBA );
   end

   isnum = prod(single(isnum));
   
   return;


