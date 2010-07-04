function [nidef] = afni_ni_defs()
   global nidef;
   
   %work out types /* definitions from niml.h 
   nidef.NI_BYTE     =   0        ;       %/* == MRI_byte    */
   nidef.NI_SHORT    =   1       ;        %/* == MRI_short   */
   nidef.NI_INT      =   2       ;        %/* == MRI_int     */
   nidef.NI_FLOAT32  =   3       ;        %/* == MRI_float   */
   nidef.NI_FLOAT    =   nidef.NI_FLOAT32;
   nidef.NI_FLOAT64  =   4      ;         %/* == MRI_double  */
   nidef.NI_DOUBLE   =   nidef.NI_FLOAT64;
   nidef.NI_COMPLEX64=   5   ;            %/* == MRI_complex */
   nidef.NI_COMPLEX  =   nidef.NI_COMPLEX64;
   nidef.NI_RGB      =   6   ;            %/* == MRI_rgb     */
   nidef.NI_RGBA     =   7  ;             %/* == MRI_rgba    */
   nidef.NI_STRING   =   8  ;             %/* after "basic" types */

%/*! Number of types of fixed size ("basic" types).
%    Note that if this changes,
%    the NI_rowtype stuff must be altered accordingly. */
   nidef.NI_NUM_BASIC_TYPES = 8;
   
%/*! One more than the last NI_ data type code defined above. */
   nidef.NI_NUM_TYPES       = 9;

   nidef.type_name = strvcat({
      'byte'  ; 'short'  ; 'int'     ;
      'float' ; 'double' ; 'complex' ;
      'rgb'   ; 'rgba'   ; 'String' });

   nidef.type_alias = strvcat({   %/* aliases */
     'uint8'   ; 'int16'   ; 'int32'     ;
     'float32' ; 'float64' ; 'complex64' ;
     'rgb8'    ; 'rgba8'   ; 'CString' });

   nidef.init = 1;
   return


