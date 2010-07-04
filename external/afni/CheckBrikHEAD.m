function [err, ErrMessage, Info] = CheckBrikHEAD (Info)
%
%   [err, ErrMessage, InfoOut] = CheckBrikHEAD (Info)
%
%Purpose:
%   Checks to determine if the fields in Info are appropriate AFNI style
%   
%   
%Input Parameters:
%   Info is a structure containing all the Header fields, 
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   ErrMessage
%   Vector Indexing in error messages is done a la matlab. 
%   so DATASET_RANK(1) in the error message refers to the first value in DATASET_RANK
%   which would be DATASET_RANK[0] in AFNI
%   InfoOut
%   same as Info    , but with the empty optional fields trimmed out
%
%Key Terms:
%   
%More Info :
%   BrikInfo for list of Mandatory fields and AFNI's README.attributes
%   
%   
%
%     Author : Ziad Saad
%     Date : Fri Apr 6 09:48:49 PDT 2001
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'CheckBrikHEAD';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

ErrMessage = '';

%Mandatory Fields Specs 

MandatoryFieldNames = 'DATASET_RANK~DATASET_DIMENSIONS~TYPESTRING~SCENE_DATA~ORIENT_SPECIFIC~ORIGIN~DELTA';
N_Mandatory = WordCount(MandatoryFieldNames, '~');

%check if the mandatory fields are present
for (im = 1:1:N_Mandatory),
	%Check that all the Mandatory Fields have been specified.
	[err, CurName] = GetWord(MandatoryFieldNames, im, '~');
	if(~isfield(Info, CurName)),
		err = 1; ErrMessage = sprintf('Error %s: Field %s must be specified for the Header to be proper.', FuncName, CurName);
		warndlg(ErrMessage);
		return;
	end
end

%Now Get all the rules
[err, ErrMessage, Rules] = HEAD_Rules; 
if (err), ErrMessage = sprintf('Error %s: Error in HEAD_RULES.\n%s', FuncName, ErrMessage); errordlg(ErrMessage); return; end
N_Rules = length(Rules);

%auto checks and cleanup of empty fields
for (ir = 1:1:N_Rules),
			%fprintf(1,'Checking %s ...', Rules(ir).Name);
	%Check that all the Remaining Fields are OK.
	if (isfield(Info, Rules(ir).Name)), 
		if (isempty(getfield(Info, Rules(ir).Name))),
			%fprintf(1,' Empty .\n');
			%remove this field if it is empty, this reduces the clutter of the ascii .HEAD file
			Info = rmfield(Info, Rules(ir).Name);
		else 
				%fprintf(1,'Not Empty .\n');
			%check for type coherence
			if (ischar(getfield(Info, Rules(ir).Name))),
				if (Rules(ir).isNum),
					err = 1; ErrMessage = sprintf('Error %s: Field %s type (string or numerical) is wrong.', FuncName, Rules(ir).Name);warndlg(ErrMessage);return;
				end
			else %a number, verify type
				if (Rules(ir).isNum == 1 & ~isint(getfield(Info, Rules(ir).Name))),
					err = 1; ErrMessage = sprintf('Error %s: Field %s type must be an integer.', FuncName, Rules(ir).Name);
				end
			end
			%check for length specs
			if (~isempty(Rules(ir).Length)),
				if (length(getfield(Info,Rules(ir).Name)) ~= Rules(ir).Length),
					err = 1; ErrMessage = sprintf('Error %s: Field %s length must be %d.', FuncName, Rules(ir).Name, Rules(ir).Length);warndlg(ErrMessage);return;
				end
			end
			%check for minimum length specs
			if (~isempty(Rules(ir).minLength)),
				if (length(getfield(Info,Rules(ir).Name)) < Rules(ir).minLength),
					err = 1; ErrMessage = sprintf('Error %s: Field %s length must be at least %d.', FuncName, Rules(ir).Name, Rules(ir).minLength);warndlg(ErrMessage);return;
				end
			end
		end
	end
end



%Further Field specfic checks
	%For Mandatory Fields
		if (Info.DATASET_RANK(1) ~= 3),
			err = 1; ErrMessage = sprintf('Error %s: DATASET_RANK(1) must be 3', FuncName); errordlg(ErrMessage); return;
		end
		if (length(Info.DATASET_RANK) < 8),
			v = zeros(1,8); v(1:length(Info.DATASET_RANK)) = Info.DATASET_RANK;
			Info.DATASET_RANK = v;
		end
		if (min(Info.DATASET_DIMENSIONS(1:3)) < 1),
			err = 1; ErrMessage = sprintf('Error %s: All values in DATASET_DIMENSIONS must be >= 1', FuncName); errordlg(ErrMessage); return;
		end
		if (length(Info.DATASET_DIMENSIONS) < 5),
			v = zeros(1,5); v(1:length(Info.DATASET_DIMENSIONS)) = Info.DATASET_DIMENSIONS;
			Info.DATASET_DIMENSIONS = v;
		end
		
		
		List = {'3DIM_HEAD_ANAT', '3DIM_HEAD_FUNC', '3DIM_GEN_ANAT', '3DIM_GEN_FUNC'}; 
		s = 0;
		for (il = 1:1:length(List)),
			s = s + strcmp(Info.TYPESTRING, char(List(il)));
		end
		if (s ~= 1),
			err = 1; ErrMessage = sprintf('Error %s: TYPESTRING must be one of \n3DIM_HEAD_ANAT, 3DIM_HEAD_FUNC, 3DIM_GEN_ANAT or 3DIM_GEN_FUNC', FuncName); errordlg(ErrMessage); return;	
		end

		if (Info.SCENE_DATA(1) < 0 |  Info.SCENE_DATA(1) > 2),
			err = 1; ErrMessage = sprintf('Error %s: SCENE_DATA(1) must be between 0 and 2', FuncName); errordlg(ErrMessage); return;
		end
		if (Info.SCENE_DATA(2) < 0 |  Info.SCENE_DATA(2) > 11),
			err = 1; ErrMessage = sprintf('Error %s: SCENE_DATA(2) must be between 0 and 11', FuncName); errordlg(ErrMessage); return;
		end
		if (Info.SCENE_DATA(3) < 0 |  Info.SCENE_DATA(3) > 3),
			err = 1; ErrMessage = sprintf('Error %s: SCENE_DATA(3) must be between 0 and 3', FuncName); errordlg(ErrMessage); return;
		end

		if (Info.ORIENT_SPECIFIC(1) < 0 |  Info.ORIENT_SPECIFIC(1) > 5),
			err = 1; ErrMessage = sprintf('Error %s: ORIENT_SPECIFIC(1) must be between 0 and 5', FuncName); errordlg(ErrMessage); return;
		end


		%For remaining fields

		if (isfield(Info, 'BRICK_STATS')),
			if (length(Info.BRICK_STATS) ~= 2.*Info.DATASET_RANK(2)), 
				err = 1; ErrMessage = sprintf('Error %s: Info.BRICK_STATS should contain %d values (%d found).', FuncName, 2.*Info.DATASET_RANK(2), length(Info.BRICK_STATS)); errordlg(ErrMessage); return; 
			end
		end

		if (isfield(Info, 'BRICK_TYPES')),
			if (length(Info.BRICK_TYPES) ~= Info.DATASET_RANK(2)), 
				err = 1; ErrMessage = sprintf('Error %s: Info.BRICK_TYPES should contain %d values (%d found).', FuncName, Info.DATASET_RANK(2), length(Info.BRICK_TYPES)); errordlg(ErrMessage); return; 	
			end
		end

		if (isfield(Info, 'BRICK_FLOAT_FACS')),
			if (length(Info.BRICK_FLOAT_FACS) ~= Info.DATASET_RANK(2)), 
				err = 1; ErrMessage = sprintf('Error %s: Info.BRICK_FLOAT_FACS should contain %d values (%d found).', FuncName, Info.DATASET_RANK(2), length(Info.BRICK_FLOAT_FACS)); errordlg(ErrMessage); return; 	
			end
		end

		if (isfield(Info, 'TAXIS_NUMS')),
			if (~isfield(Info, 'TAXIS_FLOATS')),
				err = 1; ErrMessage = sprintf('Error %s: TAXIS_FLOATS must be specified along with TAXIS_NUMS', FuncName); errordlg(ErrMessage); return;
			end
			if (Info.TAXIS_NUMS(1) ~= Info.DATASET_RANK(2)),
				err = 1; ErrMessage = sprintf('Error %s: TAXIS_NUMS(1) must be equal to DATASET_RANK(2)', FuncName); errordlg(ErrMessage); return;
			end

			if (Info.TAXIS_NUMS(2) & (~isfield(Info,'TAXIS_OFFSETS') | isempty(Info.TAXIS_OFFSETS))),
				err = 1; ErrMessage = sprintf('Error %s: TAXIS_OFFSETS must have a value if TAXIS_NUMS(2) is != 0', FuncName); errordlg(ErrMessage); return;
			end
		end

		if (isfield(Info, 'IDCODE_DATE')),
			if (length(Info.IDCODE_DATE) > 46),
				err = 1; ErrMessage = sprintf('Error %s: IDCODE_DATE must be less than 46 characters.', FuncName); errordlg(ErrMessage); return;
			end
		end

		if (isfield(Info, 'BYTEORDER_STRING')),
			if (~strcmp(Info.BYTEORDER_STRING, 'LSB_FIRST') & ~strcmp(Info.BYTEORDER_STRING, 'MSB_FIRST')),
				err = 1; ErrMessage = sprintf('Error %s: BYTEORDER_STRING must be either MSB_FIRST or LSB_FIRST', FuncName); errordlg(ErrMessage); return;
			end
		end

      if (isfield(Info, 'BRICK_LABS')),
		   ntilda = length(find(Info.BRICK_LABS == '~'));
         if (ntilda ~= Info.DATASET_RANK(2)),
			   %err = 1; ErrMessage = sprintf('Error %s: There has to be exactly %d ~ separated labels in BRICK_LABS', FuncName, Info.DATASET_RANK(2)); errordlg(ErrMessage); return;	
            if (ntilda ~= Info.DATASET_RANK(2)-1), %otherwise be quiet, it is OK
               fprintf(2, 'Warning %s: You have %d sub-bricks but %d ~ separated labels.\n', FuncName, Info.DATASET_RANK(2), ntilda); 
            end
         end	
      end

		if (isfield(Info, 'NOTES_COUNT') & (Info.NOTES_COUNT < 0 | Info.NOTES_COUNT > 999)),
			err = 1; ErrMessage = sprintf('Error %s: NOTES_COUNT must be between 0 and 999.', FuncName); errordlg(ErrMessage); return;
		end 

		if (isfield(Info, 'WARP_TYPE') & (Info.WARP_TYPE(1) ~= 0 & Info.WARP_TYPE(1) ~= 1)),
			err = 1; ErrMessage = sprintf('Error %s: WARP_TYPE(1) must be 0 or 1.', FuncName); errordlg(ErrMessage); return;
		end

		if (isfield(Info, 'WARP_DATA')),
   		if (rem(length(Info.WARP_DATA), 30)),
				err = 1; ErrMessage = sprintf('Error %s: The number of elements in WARP_DATA must be a multiple of 30.\nOtherwise BLTs would taste rotten.', FuncName); errordlg(ErrMessage); return;
			end
		end

		if (isfield(Info, 'MARKS_FLAGS')),
			if (Info.MARKS_FLAGS(1) ~= 1 & Info.MARKS_FLAGS(1) ~= 2),
				err = 1; ErrMessage = sprintf('Error %s: MARKS_FLAGS(1) must be 0 or 1.', FuncName); errordlg(ErrMessage); return;	
			end
			if (Info.MARKS_FLAGS(2) ~= 1),
				err = 1; ErrMessage = sprintf('Error %s: MARKS_FLAGS(2) must be 1.', FuncName); errordlg(ErrMessage); return;	
			end
		end

		if (isfield(Info, 'TAGSET_NUM')),
			if (Info.TAGSET_NUM(1) > 100),
				err = 1; ErrMessage = sprintf('Error %s: TAGSET_NUM(1) should not be larger than 100.', FuncName); errordlg(ErrMessage); return;
			end
		end

		if (isfield(Info, 'TAGSET_NUM') & isfield(Info, 'TAGSET_FLOATS')),
			if (length(Info.TAGSET_FLOATS) ~= (Info.TAGSET_NUM(1) .* Info.TAGSET_NUM(2)) ),   % Ask Ziad about this later
				err = 1; ErrMessage = sprintf('Error %s: Length of TAGSET_FLOATS must equal TAGSET_NUM(1) * TAGSET_NUM(2).', FuncName); errordlg(ErrMessage); return;	
			end
		end

		if (isfield(Info, 'TAGSET_NUM') & isfield(Info, 'TAGSET_LABELS')),
			if (WordCount(Info.TAGSET_LABELS, '~') ~= Info.TAGSET_NUM(1)),   % Ask Ziad later
				err = 1; ErrMessage = sprintf('Error %s: TAGSET_LABELS must contain TAGSET_NUM(1) ~ delimited strings.', FuncName); errordlg(ErrMessage); return;	
			end
		end

		if (isfield(Info, 'BRICK_TYPES')),
			if (~isempty(unique(setdiff(Info.BRICK_TYPES, [0 1 3 5]))))
				err = 1; ErrMessage = sprintf('Error %s: BRICK_TYPES only types 0 1 3 5 are supported.', FuncName); errordlg(ErrMessage); return;	
			end
		end

		%checks on the rotate fields
		if (isfield(Info, 'VOLREG_ROTCOM_NUM')),
			for (i=1:1:Info.VOLREG_ROTCOM_NUM),
					sp = pad_strn(num2str(i-1), '0', 6, 1);
					spad = sprintf('~isfield(Info, ''VOLREG_MATVEC_%s'') | isempty(Info.VOLREG_MATVEC_%s) | ischar(Info.VOLREG_MATVEC_%s) | length(Info.VOLREG_MATVEC_%s) ~= 12', sp, sp, sp, sp);
					if (eval([ spad ])),
						err = 1; ErrMessage = sprintf('Error %s: Undefined VOLREG_MATVEC_%s or bad type (strings) or bad length (~=12)', FuncName, sp);errordlg(ErrMessage); return;  
					end
			end
		end

		%check on the NOT_NUMBER_... fields
		if (~isfield(Info, 'NOTES_COUNT') & isfield(Info, 'NOTE_NUMBER_001')),
			Info = rmfield(Info, 'NOTE_NUMBER_001');
		end

		%check on the HISTORY field
		if (isfield(Info, 'HISTORY_NOTE')),
			vtmp = unique(double(Info.HISTORY_NOTE));
			if (~isempty(intersect([13 10 34 9 7 11 8], vtmp))),
				err = 1; ErrMessage = sprintf('Error %s: HISTORY_NOTE has not been TomRossed (formatted) properly.', FuncName); errordlg(ErrMessage);return;  				
			end
		end
		
%check on field groups

err = 0;
return;

