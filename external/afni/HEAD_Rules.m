function [err, ErrMessage, Rules] = HEAD_Rules (FieldNames)
%
%   [err, ErrMessage, Rules] = HEAD_Rules (FieldNames)
%
%Purpose:
%   Retrieve the rules for the different fields in AFNI
%   
%   
%Input Parameters:
%   FieldNames is a ~ delirited list of N field names
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   Rules is a vector of N structures with the following fields
%     .Name is the field's name
%     .isNum if it is supposed to be a string or a number (0 (string)/1 (int)/2(float), -1 for fields not in the database)
%     .Length length of number vector or string (leave empty if not applicable)
%     .minLength minimum length of vector or string (leave empty if not applicable)
%   
%      
%Key Terms:
%   
%More Info :
%   AFNI's README.attributes
%   
%   
%
%     Author : Ziad Saad
%     Date : Mon Apr 9 15:58:24 PDT 2001
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
	FuncName = 'HEAD_Rules';
%Debug Flag
	DBG = 1;

%initailize return variables
	err = 1;
	ErrMessage = [];
	Rules = [];


if (nargin == 0),
	FieldNames = '';
end

if (isempty(FieldNames)),
	%Need to retrieve rules for all fields 
	%This is a ~ delimited list of all the fields present in this function
		tmp1 = 'BRICK_FLOAT_FACS~BRICK_LABS~BRICK_KEYWORDS~BRICK_STATAUX~BRICK_STATS~BRICK_TYPES~BYTEORDER_STRING~';
		tmp2 = 'DATASET_DIMENSIONS~DATASET_KEYWORDS~DATASET_NAME~DATASET_RANK~DELTA~HISTORY_NOTE~IDCODE_ANAT_PARENT~';
		tmp3 = 'IDCODE_DATE~IDCODE_STRING~IDCODE_WARP_PARENT~LABEL_1~LABEL_2~MARKS_FLAGS~MARKS_HELP~MARKS_LAB~MARKS_XYZ~';			
		tmp4 = 'NOTES_COUNT~ORIENT_SPECIFIC~ORIGIN~SCENE_DATA~STAT_AUX~TAGALIGN_MATVEC~TAGSET_FLOATS~TAGSET_LABELS~TAGSET_NUM~';
		tmp5 = 'TAXIS_FLOATS~TAXIS_NUMS~TAXIS_OFFSETS~TO3D_ZPAD~TYPESTRING~VOLREG_BASE_IDCODE~VOLREG_BASE_NAME~VOLREG_CENTER_BASE~';
		tmp6 = 'VOLREG_CENTER_OLD~VOLREG_GRIDPARENT_IDCODE~VOLREG_GRIDPARENT_NAME~VOLREG_INPUT_IDCODE~VOLREG_INPUT_NAME~';
		tmp7 = 'VOLREG_ROTCOM_NUM~VOLREG_ROTPARENT_IDCODE~VOLREG_ROTPARENT_NAME~WARP_DATA~WARP_TYPE~';
		AllFields = zdeblank(sprintf('%s%s%s%s%s%s%s', tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7));
		FieldNames = AllFields;
else
	FieldNames = zdeblank (FieldNames);
end

%remove last ~ delimiter
	Ntmp = length(FieldNames);
	if (FieldNames(Ntmp) == '~'), FieldNames = FieldNames(1:Ntmp-1); end

Nfields = WordCount(FieldNames, '~');
%initialize structure fields
	Rules(Nfields).Name = '';
	Rules(Nfields).isNum = [];
	Rules(Nfields).Length = [];
	Rules(Nfields).minLength = [];

for (ir=1:1:Nfields),
	[err,CurName] = GetWord(FieldNames, ir, '~');
	isNum = [];
	Length = [];
	minLength = [];
	if (~isempty(CurName)),
		switch CurName,
			case 'BRICK_FLOAT_FACS',
				isNum = 2;
			case 'BRICK_LABS',
				isNum = 0;
			case 'BRICK_KEYWORDS',
				isNum = 0;
			case 'BRICK_STATAUX',
				isNum = 2;
			case 'BRICK_STATS',
				isNum = 2;
			case 'BRICK_TYPES',
				isNum = 1;
			case 'BYTEORDER_STRING',
				isNum = 0;
			case 'DATASET_DIMENSIONS',
				isNum = 1;
				minLength = 3;
			case 'DATASET_KEYWORDS',
				isNum = 0;
			case 'DATASET_NAME',
				isNum = 0;
			case 'DATASET_RANK',
				isNum = 1;
				minLength = 2;
			case 'DELTA',
				isNum = 2;
				Length = 3;
			case 'HISTORY_NOTE',
				isNum = 0;
			case 'IDCODE_ANAT_PARENT',
				isNum = 0;
			case 'IDCODE_DATE',
				isNum = 0;
			case 'IDCODE_STRING',
				isNum = 0;
			case 'IDCODE_WARP_PARENT',
				isNum = 0;
			case 'LABEL_1',
				isNum = 0;
			case 'LABEL_2',
				isNum = 0;
			case 'MARKS_FLAGS',
				isNum = 1;
				minLength = 2;
			case 'MARKS_HELP',
				isNum = 0;
				Length = 2559;
			case 'MARKS_LAB',
				isNum = 0;
				Length = 199;
			case 'MARKS_XYZ',
				isNum = 2;
				Length = 30;
			case 'NOTES_COUNT',
				isNum = 1;
				Length = 1;
			case 'ORIENT_SPECIFIC',
				isNum = 1;
				Length = 3;
			case 'ORIGIN',
				isNum = 2;
				Length = 3;
			case 'SCENE_DATA',
				isNum = 1;
				minLength = 3;
			case 'STAT_AUX',
				isNum = 2;
			case 'TAGALIGN_MATVEC',
				isNum = 2;
				Length = 12;
			case 'TAGSET_FLOATS',
				isNum = 2;
			case 'TAGSET_LABELS',
				isNum = 0;
			case 'TAGSET_NUM',
				isNum = 1;
				Length = 2;
			case 'TAXIS_FLOATS',
				isNum = 2;
				minLength = 5;
			case 'TAXIS_NUMS',
				isNum = 1;
				minLength = 3;
			case 'TAXIS_OFFSETS',
				isNum = 2;
			case 'TO3D_ZPAD',
				isNum = 1;
				Length = 3;
			case 'TYPESTRING',
				isNum = 0;
			case 'VOLREG_BASE_IDCODE',
				isNum = 0;
				Length = 26;
			case 'VOLREG_BASE_NAME',
				isNum = 0;
			case 'VOLREG_CENTER_BASE',
				isNum = 2;
				Length = 3;
			case 'VOLREG_CENTER_OLD',
				isNum = 2;
				Length = 3;
			case 'VOLREG_GRIDPARENT_IDCODE',
				isNum = 0;
			case 'VOLREG_GRIDPARENT_NAME',
				isNum = 0;
			case 'VOLREG_INPUT_IDCODE',
				isNum = 0;
				Length = 26;
			case 'VOLREG_INPUT_NAME',
				isNum = 0;
			case 'VOLREG_ROTCOM_NUM',
				isNum = 1;
				Length = 1;
			case 'VOLREG_ROTPARENT_IDCODE',
				isNum = 0;
				Length = 26;
			case 'VOLREG_ROTPARENT_NAME',
				isNum = 0;
			case 'WARP_DATA',
				isNum = 2;
			case 'WARP_TYPE',
				isNum = 1;
				minLength = 2;
			otherwise,
				isNum = -1;
		end
		Rules(ir).Name = CurName;
		Rules(ir).isNum = isNum;
		Rules(ir).Length = Length;
		Rules(ir).minLength = minLength;
	else
		err = 1; ErrMessage = sprintf('Error %s: Empty String in Field Name List', FuncName); errordlg(ErrMessage); return;	
	end %isempty(CurName)
end %for cnt



	 
 


 


err = 0;
ErrMessage = '';

return;

