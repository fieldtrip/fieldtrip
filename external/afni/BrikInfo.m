function [err,Info, BRIKinfo] = BrikInfo (BrickName)
%
%   [err,Info] = BrikInfo (BrickName)
%
%Purpose:
%   
%   returns some field values in the .HEAD files
%   
%Input Parameters:
%   Brick filename 
%  (also works for 1D files but Info has very limited info.
%   You should really avoid using BrikInfo on 1D file
%   because the 1D file is read in its entirety (and then cleared)
%   before Info can be returned. Better used BrikLoad directly or
%   Read_1D)
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   Info is a structure with the following fields
%           Fieldnames in UPPER CASE correspond to the fields in AFNI
%           Field names in lower case correspond to some interpretation of the uppercase ones
%      .RootName : the brik name without .HEAD or .BRIK. 
%      .TypeName: The interpretation of BRICK_TYPES (short, int, floats, etc...)
%      .TypeBytes : The corresponding byte size of BRICK_TYPES
%		 .ByteOrder: interpretation of BYTEORDER_STRING
%      .Orientation : interpretation of ORIENT_SPECIFIC, Orientation(:,1)' forms the three letter orientation code 
%        The ith row of Orientation describes the orientation along that dimension.
%   see AFNI's README.attributes for a complete description of these fields
%   MANDATORY FIELDS AS SPECIFIED IN: ~cox/README.attributes (added April 6 2001)
%      .DATASET_RANK : ASK BOB. I think the first is always 3 (3D) and the other is the number of sub-bricks
%      .DATASET_DIMENSIONS : Number of voxels in each of the three dimensions
%      .TYPESTRING: Determines if the data set is of type Anat (underlay) or Func (overlay)
%      .SCENE_DATA : The three integer codes describing the dataset type                
%		 .ORIENT_SPECIFIC : orintation code 
%      .ORIGIN : The xyz coordinates of the center of the (0, 0, 0) voxel in the data set
%      .DELTA : the increment (in mm) to go from one voxel to the next (could be +ve or -ve depending on slices)
%      .TAXIS_NUMS: see readme file
%      .TAXIS_FLOATS
%      .TAXIS_OFFSETS
%   ALMOST MANDATORY FIELDS   
%      .IDCODE_STRING
%      .IDCODE_DATE
%      .BYTEORDER_STRING : Byte order string
%      .BRICK_STATS : range of values in brick (min to max) 
%         (Do not apply the scaling factor to these values)
%      .BRICK_TYPES : types of values in .BRIK 
%      .BRICK_FLOAT_FACS : float factors to apply to bricks for recovering original values
%      .BRICK_LABS : The Sub-brick labels (~delimited)
%      .BRICK_STATAUX : Auxilliary Statistical Information
%      .STAT_AUX
%      .HISTORY_NOTE
%      .NOTES_COUNT
%      .NOTE_NUMBER_001
%      .TAGALIGN_MATVEC
%      .VOLREG_CENTER_OLD
%      .VOLREG_CENTER_BASE
%      .VOLREG_ROTPARENT_IDCODE
%      .VOLREG_ROTPARENT_NAME
%      .VOLREG_GRIDPARENT_IDCODE
%      .VOLREG_GRIDPARENT_NAME
%      .VOLREG_INPUT_IDCODE
%      .VOLREG_INPUT_NAME
%      .VOLREG_BASE_IDCODE
%      .VOLREG_BASE_NAME
%      .VOLREG_ROTCOM_NUM
%      .VOLREG_MATVEC_xxxxxx
%      .VOLREG_ROTCOM_xxxxxx      
%      .IDCODE_ANAT_PARENT
%      .TO3D_ZPAD
%      .IDCODE_WARP_PARENT
%      .WARP_TYPE
%      .WARP_DATA
%      .MARKS_XYZ
%      .MARKS_LAB
%      .MARKS_HELP 
%      .MARKS_FLAGS
%      .TAGSET_NUM  
%      .TAGSET_FLOATS
%      .TAGSET_LABELS 
%      .LABEL_1
%      .LABEL_2 
%      .DATASET_NAME
%      .DATASET_KEYWORDS
%      .BRICK_KEYWORDS 
%      .HISTORY_NOTE
%      .NOTES_COUNT
%      .NOTE_NUMBER_xxx
%      
%  To implement in the future,
%      VOLREG_MATVEC_xxxxxx , VOLREG_ROTCOM_xxxxxx
%      
%  The following fields were added to support 1D file format. 
%      .FileFormat: 'BRIK' ('1D' is allowed for 1D files but 1D files do not use BrikInfo)
%      .Extension_1D: The extension of the 1D filename
%  BrikLoad and WriteBrik now read and write 1D files
%
%Key Terms:
%   The 1st, second and third dimensions refer to the dimensions the slices were entered into to3d 
%More Info :
%    afni's README.attributes
%
%     Author : Ziad Saad
%     Date : Sun Oct 17 00:13:49 CDT 1999 , Major modification, April 06 2001


%Define the function name for easy referencing
FuncName = 'BrikInfo';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
Info = [];

%check if this is a 1D file ..
is1D = 0;
[St, xtr] = Remove1DExtension(BrickName);
if (~isempty(xtr)),
   is1D = 1;
end

if (is1D), % 1D land
   [err, V, Info] = Read_1D(BrickName, 1);
   if (err), 
      ErrMessage = sprintf ('%s: Failed to read %s file', FuncName, BrickName);
      err = ErrEval(FuncName,'Err_1D file could not be read');
      return;
   end
   clear V; %biggest waste since 2002
   return;
end

	
	vtmp = findstr(BrickName,'.BRIK');
	if (~isempty(vtmp)), %remove .BRIK
		BrickName = BrickName(1:vtmp(1)-1);
	end
	
	vtmp = findstr(BrickName,'.HEAD');
	if (isempty(vtmp)), %add .HEAD
		BrickName = sprintf('%s.HEAD',BrickName);
	end
	
	
	if (exist(BrickName) ~= 2),
		err = ErrEval(FuncName,'Err_Could not find data set');
		return;
	end
	
	%store the name without the extension
		vtmp = findstr(BrickName,'.HEAD');
		Info.RootName = BrickName(1:vtmp-1);

	fidin=fopen(BrickName,'r');
	if (fidin < 0),	err = ErrEval(FuncName,'Err_Could not read .HEAD file'); return;	end
		BRIKinfo = fscanf(fidin,'%c');
	fclose(fidin);
	N_BRIKinfo = length(BRIKinfo);

%get subBrik info
	%Nx, Ny, Nz
		[err,Info.DATASET_DIMENSIONS] = BrikInfo_SectionValue (BRIKinfo, 'DATASET_DIMENSIONS');
	
	%DATASET_RANK
		[err, Info.DATASET_RANK] = BrikInfo_SectionValue(BRIKinfo, 'DATASET_RANK');

	%BRICK_TYPES
		[err, Info.BRICK_TYPES] = BrikInfo_SectionValue(BRIKinfo, 'BRICK_TYPES');
	
	%BRICK_STATS
		[err, Info.BRICK_STATS] = BrikInfo_SectionValue(BRIKinfo, 'BRICK_STATS');

	%BRICK_FLOAT_FACS
		[err, Info.BRICK_FLOAT_FACS] = BrikInfo_SectionValue(BRIKinfo, 'BRICK_FLOAT_FACS');
	
	%BYTEORDER_STRING
		[err, Info.BYTEORDER_STRING] = BrikInfo_SectionValue(BRIKinfo, 'BYTEORDER_STRING');
			
	%ORIENT_SPECIFIC
		[err, Info.ORIENT_SPECIFIC] = BrikInfo_SectionValue(BRIKinfo, 'ORIENT_SPECIFIC');
		
	%ORIGIN
		[err, Info.ORIGIN] = BrikInfo_SectionValue(BRIKinfo, 'ORIGIN');
	
	%DELTA
		[err, Info.DELTA] = BrikInfo_SectionValue(BRIKinfo, 'DELTA');	
	
	%BRICK_LABS
		[err, Info.BRICK_LABS] = BrikInfo_SectionValue(BRIKinfo, 'BRICK_LABS');	
	
	%BRICK_KEYWORDS
		[err, Info.BRICK_KEYWORDS] = BrikInfo_SectionValue(BRIKinfo, 'BRICK_KEYWORDS');	
		
	%SCENE_DATA
		[err, Info.SCENE_DATA] = BrikInfo_SectionValue(BRIKinfo, 'SCENE_DATA');	

	%TYPESTRING
		[err, Info.TYPESTRING] = BrikInfo_SectionValue(BRIKinfo, 'TYPESTRING');	

	%TAXIS_NUMS
		[err, Info.TAXIS_NUMS] = BrikInfo_SectionValue(BRIKinfo, 'TAXIS_NUMS');	
	
	%TAXIS_FLOATS
		[err, Info.TAXIS_FLOATS] = BrikInfo_SectionValue(BRIKinfo, 'TAXIS_FLOATS');	
		
	%TAXIS_OFFSETS
		[err, Info.TAXIS_OFFSETS] = BrikInfo_SectionValue(BRIKinfo, 'TAXIS_OFFSETS');	
	
	%IDCODE_STRING
		[err, Info.IDCODE_STRING] = BrikInfo_SectionValue(BRIKinfo, 'IDCODE_STRING');	
	
	%IDCODE_DATE
		[err, Info.IDCODE_DATE] = BrikInfo_SectionValue(BRIKinfo, 'IDCODE_DATE');	
	
	%BRICK_STATAUX
		[err, Info.BRICK_STATAUX] = BrikInfo_SectionValue(BRIKinfo, 'BRICK_STATAUX');	
	
	%STAT_AUX
		[err, Info.STAT_AUX] = BrikInfo_SectionValue(BRIKinfo, 'STAT_AUX');	
	
	%HISTORY_NOTE
		[err, Info.HISTORY_NOTE] = BrikInfo_SectionValue(BRIKinfo, 'HISTORY_NOTE');	
	
	%NOTES_COUNT
		[err, Info.NOTES_COUNT] = BrikInfo_SectionValue(BRIKinfo, 'NOTES_COUNT');	
	
	%NOTE_NUMBER_001
		[err, Info.NOTE_NUMBER_001] = BrikInfo_SectionValue(BRIKinfo, 'NOTE_NUMBER_001');	
	
	%TAGALIGN_MATVEC
		[err, Info.TAGALIGN_MATVEC] = BrikInfo_SectionValue(BRIKinfo, 'TAGALIGN_MATVEC');	
	
	%VOLREG_CENTER_OLD
		[err, Info.VOLREG_CENTER_OLD] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_CENTER_OLD');	
	
	%VOLREG_CENTER_BASE
		[err, Info.VOLREG_CENTER_BASE] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_CENTER_BASE');	
	
	%VOLREG_ROTPARENT_IDCODE
		[err, Info.VOLREG_ROTPARENT_IDCODE] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_ROTPARENT_IDCODE');	
	
	%VOLREG_ROTPARENT_NAME
		[err, Info.VOLREG_ROTPARENT_NAME] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_ROTPARENT_NAME');	
	
	%VOLREG_GRIDPARENT_IDCODE
		[err, Info.VOLREG_GRIDPARENT_IDCODE] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_GRIDPARENT_IDCODE');	
	
	%VOLREG_GRIDPARENT_NAME
		[err, Info.VOLREG_GRIDPARENT_NAME] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_GRIDPARENT_NAME');	
	
	%VOLREG_INPUT_IDCODE
		[err, Info.VOLREG_INPUT_IDCODE] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_INPUT_IDCODE');	
	
	%VOLREG_INPUT_NAME
		[err, Info.VOLREG_INPUT_NAME] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_INPUT_NAME');	
	
	%VOLREG_BASE_IDCODE
		[err, Info.VOLREG_BASE_IDCODE] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_BASE_IDCODE');	
	
	%VOLREG_BASE_NAME
		[err, Info.VOLREG_BASE_NAME] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_BASE_NAME');	
	
	%VOLREG_ROTCOM_NUM
		[err, Info.VOLREG_ROTCOM_NUM] = BrikInfo_SectionValue(BRIKinfo, 'VOLREG_ROTCOM_NUM');	
	
	if (~err),
		for (i=1:1:Info.VOLREG_ROTCOM_NUM),
			sp = pad_strn(num2str(i-1), '0', 6, 1);
			spad = sprintf('[err, Info.VOLREG_MATVEC_%s] = BrikInfo_SectionValue(BRIKinfo, ''VOLREG_MATVEC_%s'');' , sp, sp);
			eval([ spad ]);
			spad = sprintf('[err, Info.VOLREG_ROTCOM_%s] = BrikInfo_SectionValue(BRIKinfo, ''VOLREG_ROTCOM_%s'');' , sp, sp);
			eval([ spad ]);
		end
	end	
   
		    
	%IDCODE_ANAT_PARENT
		[err, Info.IDCODE_ANAT_PARENT] = BrikInfo_SectionValue(BRIKinfo, 'IDCODE_ANAT_PARENT');	
	
	%TO3D_ZPAD
		[err, Info.TO3D_ZPAD] = BrikInfo_SectionValue(BRIKinfo, 'TO3D_ZPAD');	
	
	%IDCODE_WARP_PARENT
		[err, Info.IDCODE_WARP_PARENT] = BrikInfo_SectionValue(BRIKinfo, 'IDCODE_WARP_PARENT');	
	
	%WARP_TYPE
		[err, Info.WARP_TYPE] = BrikInfo_SectionValue(BRIKinfo, 'WARP_TYPE');	
	
	%WARP_DATA
		[err, Info.WARP_DATA] = BrikInfo_SectionValue(BRIKinfo, 'WARP_DATA');	
	
	%MARKS_XYZ
		[err, Info.MARKS_XYZ] = BrikInfo_SectionValue(BRIKinfo, 'MARKS_XYZ');	
	
	%MARKS_LAB
		[err, Info.MARKS_LAB] = BrikInfo_SectionValue(BRIKinfo, 'MARKS_LAB');	
	
	%MARKS_HELP
		[err, Info.MARKS_HELP] = BrikInfo_SectionValue(BRIKinfo, 'MARKS_HELP');	
	
	%MARKS_FLAGS
		[err, Info.MARKS_FLAGS] = BrikInfo_SectionValue(BRIKinfo, 'MARKS_FLAGS');	
	
	%TAGSET_NUM
		[err, Info.TAGSET_NUM] = BrikInfo_SectionValue(BRIKinfo, 'TAGSET_NUM');	
	
	%TAGSET_FLOATS
		[err, Info.TAGSET_FLOATS] = BrikInfo_SectionValue(BRIKinfo, 'TAGSET_FLOATS');	
	
	%TAGSET_LABELS
		[err, Info.TAGSET_LABELS] = BrikInfo_SectionValue(BRIKinfo, 'TAGSET_LABELS');	
	
	%LABEL_1
		[err, Info.LABEL_1] = BrikInfo_SectionValue(BRIKinfo, 'LABEL_1');	
	
	%LABEL_2
		[err, Info.LABEL_2] = BrikInfo_SectionValue(BRIKinfo, 'LABEL_2');	
	
	%DATASET_NAME
		[err, Info.DATASET_NAME] = BrikInfo_SectionValue(BRIKinfo, 'DATASET_NAME');	
	
	%DATASET_KEYWORDS
		[err, Info.DATASET_KEYWORDS] = BrikInfo_SectionValue(BRIKinfo, 'DATASET_KEYWORDS');	
	
	%BRICK_KEYWORDS
		[err, Info.BRICK_KEYWORDS] = BrikInfo_SectionValue(BRIKinfo, 'BRICK_KEYWORDS');	
	
	%HISTORY_NOTE
		%[err, Info.HISTORY_NOTE] = BrikInfo_SectionValue(BRIKinfo, 'HISTORY_NOTE');	
	
	%NOTES_COUNT
		%[err, Info.NOTES_COUNT] = BrikInfo_SectionValue(BRIKinfo, 'NOTES_COUNT');	
	
	if (~err),
		for (i=1:1:Info.NOTES_COUNT),
			sp = pad_strn(num2str(i), '0', 3, 1);
			spad = sprintf('[err, Info.NOTE_NUMBER_%s] = BrikInfo_SectionValue(BRIKinfo, ''NOTE_NUMBER_%s'');' , sp, sp);
			eval([ spad ]);
		end
	end	
	             
	%SOMETHING
		%[err, Info.SOMETHING] = BrikInfo_SectionValue(BRIKinfo, 'SOMETHING');	
      
	%DOF
		[err, Info.WORSLEY_DF] = BrikInfo_SectionValue(BRIKinfo, 'WORSLEY_DF');	
      
	%NONJ
		[err, Info.WORSLEY_NCONJ] = BrikInfo_SectionValue(BRIKinfo, 'WORSLEY_NCONJ');	
	
	%FWHM
		[err, Info.WORSLEY_FWHM] = BrikInfo_SectionValue(BRIKinfo, 'WORSLEY_FWHM');	
	
	itype = unique(Info.BRICK_TYPES);
 
 	if (length(itype) > 1),
		Info.TypeName = 'Mutliple Types';
	else
		switch itype,
			case 0
				Info.TypeName = 'byte';
				Info.TypeBytes = 1;
			case 1
				Info.TypeName = 'short';
				Info.TypeBytes = 2; %Platform dependent
			case 2
				Info.TypeName = 'int';
				Info.TypeBytes = 4; %Platform dependent
			case 3
				Info.TypeName = 'float';
				Info.TypeBytes = 4; %Platform dependent
			otherwise
				Info.TypeName = 'Dunno';
				Info.TypeBytes = 0;
		end
	end

	 if (isempty(Info.BYTEORDER_STRING)), %field not found go native
	 	Info.ByteOrder = 'unspecified';
	 else
		 if (~isempty(strmatch('MSB_FIRST', Info.BYTEORDER_STRING))),
				Info.ByteOrder = 'ieee-be'; %Big Endian
		else 
			if (~isempty(strmatch('LSB_FIRST', Info.BYTEORDER_STRING))),
					Info.ByteOrder = 'ieee-le'; %Little Endian
			else
				err = ErrEval(FuncName,'Err_Could not understand BYTEORDER_STRING'); 
				return;	
			end
		end
	 end
	
	for (i=1:1:3),
		switch Info.ORIENT_SPECIFIC(i)
			case 0
				Info.Orientation(i,:) = 'RL';	%right to left 
			case 1
				Info.Orientation(i,:) = 'LR';
			case 2
				Info.Orientation(i,:) = 'PA';
			case 3
				Info.Orientation(i,:) = 'AP';
			case 4
				Info.Orientation(i,:) = 'IS';
			case 5
				Info.Orientation(i,:) = 'SI';
			otherwise,
				err = ErrEval(FuncName,'Err_Cannot understand Orientation code');
				return;
		end
	end
	
   Info.FileFormat = 'BRIK';
   Info.Extension_1D = '';	
err = 0;
return;

