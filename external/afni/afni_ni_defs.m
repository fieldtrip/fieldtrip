function ni_def = afni_ni_defs()
    persistent cached_ni_def

    % for efficency reasons generate the definitions only once;
    % upon each subsequent call the cached data is used
    if ~isstruct(cached_ni_def)
        cached_ni_def=generate_ni_def();
    end

    ni_def=cached_ni_def;

function ni_def=generate_ni_def()
    ni_def=struct();
    %work out types /* definitions from niml.h
    ni_def.NI_BYTE     =   0        ;       %/* == MRI_byte    */
    ni_def.NI_SHORT    =   1       ;        %/* == MRI_short   */
    ni_def.NI_INT      =   2       ;        %/* == MRI_int     */
    ni_def.NI_FLOAT32  =   3       ;        %/* == MRI_float   */
    ni_def.NI_FLOAT    =   ni_def.NI_FLOAT32;
    ni_def.NI_FLOAT64  =   4      ;         %/* == MRI_double  */
    ni_def.NI_DOUBLE   =   ni_def.NI_FLOAT64;
    ni_def.NI_COMPLEX64=   5   ;            %/* == MRI_complex */
    ni_def.NI_COMPLEX  =   ni_def.NI_COMPLEX64;
    ni_def.NI_RGB      =   6   ;            %/* == MRI_rgb     */
    ni_def.NI_RGBA     =   7  ;             %/* == MRI_rgba    */
    ni_def.NI_STRING   =   8  ;             %/* after "basic" types */

    %/*! Number of types of fixed size ("basic" types).
    %    Note that if this changes,
    %    the NI_rowtype stuff must be altered accordingly. */
    ni_def.NI_NUM_BASIC_TYPES = 8;

    %/*! One more than the last NI_ data type code defined above. */
    ni_def.NI_NUM_TYPES       = 9;

    ni_def.type_name = strvcat({
      'byte'  ; 'short'  ; 'int'     ;
      'float' ; 'double' ; 'complex' ;
      'rgb'   ; 'rgba'   ; 'String' });

    ni_def.type_alias = strvcat({   %/* aliases */
     'uint8'   ; 'int16'   ; 'int32'     ;
     'float32' ; 'float64' ; 'complex64' ;
     'rgb8'    ; 'rgba8'   ; 'CString' });

    % TODO: see how to deal with other types
    undef=@()error('Unsupported type');

    ni_def.type = {...
       @uint8  ; @int16    ; @int32      ; ...
       @single ; @double   ; undef       ; ...
       undef   ; undef     ; @char       };

    ni_def.size_bytes=[...
       1       ; 2          ; 4          ; ...
       4       ; 8          ; NaN        ; ...
       NaN     ; NaN        ; 1          ];

    ni_def.init = 1;

    ni_def.escape={'&lt;'  ,'<' ;...
                   '&gt;'  ,'>' ;...
                   '&quot;','"' ;...
                   '&apos;','''';...
                   '&amp;' ,'&' };



