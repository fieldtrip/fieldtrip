function [err, ErrMessage, Info, V] = ResliceBrick (BrickName, NewCode, Opt)
%
%   [err, ErrMessage, Info, M] = ResliceBrick (BrickName, NewCode, [Opt])
%
%Purpose:
%  Reslices an AFNI brick kinda a la 3daxialize but in any orientation desired
%  No interpolation or shananigans are done 
%   
%Input Parameters:
%   BrickName the name of the afni brick
%   NewCode The orientation code to write the new brick as 'RPS' or 'RAI' 
%     or the name of a brick to reslice the data as. The new code is extracted from
%     that brick's header.
%   Opt is an optional options structure with the following fields
%     .Prefix : (optional) prefix of the output brick. Default is ''
%               in which case, the output prefix is the input prefix_<NewCode>
%     .WriteBrick : 0/[1] optional flag for writing the brick to disk
%               you may choose not to if you are passing it back to the 
%               calling function.
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   ErrMessage: a string containing warnings or error messages
%   Info : The header of the resliced brick
%   M : the resliced brick
%      
%Key Terms:
%   
%More Info :
%   
%   example: The images making up data set AMzsspgrco+orig were acquired coronally.
%         In certain instances (especially those involving 3dvolreg), you need to have
%         the images making up the different volumes to be in aqcuired in the same plane
%         and orientation. In this example, we reslice AMzsspgrco+orig to look like its images
%         were acquired just like those making up AHzsspgrax+orig (axial, LAI) 
%   
%     [err, ErrMessage] = ResliceBrick ('AMzsspgrco+orig.BRIK', 'AHzsspgrax+orig');
%
%     Author : Ziad Saad
%     Date : Fri May 11 14:10:50 PDT 2001 / latest update Oct 11 01
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
	FuncName = 'ResliceBrick';

%Debug Flag
	DBG = 1;

%initailize return variables
	err = 1;
	ErrMessage = '';

%make sure Code is OK
	tmp = WordNumber(AllOrient , NewCode, '|');
	if (isempty(tmp)), 
		%maybe it's a Brick
		CodeBrick = RemoveExtension(NewCode, '.HEAD|.BRIK');
		CodeBrick = sprintf('%s.HEAD', CodeBrick);
		if (exist(CodeBrick) == 2),
			[err,InfoCode] = BrikInfo (CodeBrick);
			if (err), ErrMessage = sprintf('Error %s subfunction: %s has a bad header', FuncName, CodeBrick ); return; end
			NewCode = InfoCode.ORIENT_SPECIFIC;
			%change it to string
			[err, NewCode] = AFNI_OrientCode (NewCode);
		else
			ErrMessage = sprintf('Error %s subfunction: %s . %s not a valid orientation code nore a brick.' ,FuncName, ErrMessage, NewCode); return; 
		end
	end
	
%Clean up the input Brick's header
	BrickName = RemoveExtension(BrickName, '.HEAD|.BRIK');
	[BrickName_Prefix, BrickName_View] = RemoveExtension(BrickName, '+orig|+acpc|+tlrc');

%make sure Opt is OK
	if (nargin == 2), 
		Opt.Prefix = '';
		Opt.verbose = 1;
		Opt.WriteBrick = 1;
	else
		if (~isfield(Opt, 'Prefix') | isempty(Opt.Prefix)), Opt.Prefix = ''; end
		if (~isfield(Opt, 'verbose') | isempty(Opt.verbose)), Opt.verbose = 1; end
		if (~isfield(Opt, 'WriteBrick') | isempty(Opt.WriteBrick)), Opt.WriteBrick = 1; end
	end
	
	if ((~isfield(Opt,'WriteBrick') | ~Opt.WriteBrick) & nargout < 4),
		warndlg(sprintf('Warning %s: From the number of output arguments and the Opt.WriteBrick, it does not look like you are planning on doing much.', FuncName));
		drawnow;
	end
	
	if (isempty(Opt.Prefix)),
		Opt.Prefix = sprintf('%s_%s', BrickName_Prefix, NewCode);
	end
		Info.RootName  = sprintf('%s%s', Opt.Prefix, BrickName_View);
		
	
	if ( (exist(sprintf('%s.HEAD', Info.RootName)) == 2 ) | (exist(sprintf('%s.BRIK', Info.RootName)) == 2 ) )
		ErrMessage = sprintf ('Error %s: Brick %s exists. Will not overwrite', FuncName, Info.RootName);
	end

%load the brick
	if (Opt.verbose),
		fprintf(1,'Loading brick %s from disk ...', BrickName);
	end
	[err, V, Info , ErrMessage] = BrikLoad(sprintf('%s', BrickName));
	if (err), ErrMessage = sprintf('Error %s subfunction: %s ', FuncName, ErrMessage); return; end
	if (Opt.verbose),
		fprintf(1,'\tDone.\n');
	end

%determine the required transformation
	if (Opt.verbose),
		[err, OldCode] = AFNI_OrientCode (Info.ORIENT_SPECIFIC);
		fprintf(1,'Reslicing %s to %s ...\n', OldCode, NewCode);
	end

	[err, maplocation, mapsign, Mtrans] = AFNI_CoordChange (Info.ORIENT_SPECIFIC, NewCode);

	V = permute(V, maplocation);
	Info_orig = Info;
	
	Info.DELTA = Info_orig.DELTA(maplocation).* mapsign;
	
	if (mapsign(1) < 0),
		N = size(V,1);
		V(:,:,:) = V(N:-1:1, :, :);
		Info.ORIGIN(1) = Info_orig.ORIGIN(maplocation(1)) + (Info_orig.DATASET_DIMENSIONS(maplocation(1)) -1).* Info_orig.DELTA(maplocation(1));
	else
		Info.ORIGIN(1) = Info_orig.ORIGIN(maplocation(1));
	end

	if (mapsign(2) < 0),
		N = size(V,2);
		V(:,:,:) = V(:, N:-1:1, :);
		Info.ORIGIN(2) = Info_orig.ORIGIN(maplocation(2)) + (Info_orig.DATASET_DIMENSIONS(maplocation(2)) -1).* Info_orig.DELTA(maplocation(2));
	else
		Info.ORIGIN(2) = Info_orig.ORIGIN(maplocation(2));
	end

	if (mapsign(3) < 0),
		N = size(V,3);
		V(:,:,:) = V(:, :, N:-1:1);
		Info.ORIGIN(3) = Info_orig.ORIGIN(maplocation(3)) + (Info_orig.DATASET_DIMENSIONS(maplocation(3)) -1) .* Info_orig.DELTA(maplocation(3));
	else
		Info.ORIGIN(3) = Info_orig.ORIGIN(maplocation(3));
	end

%setup the new header and write out the files
	%FIXED: Oct. 11 01 . THIS HERE NEEDS FIXING. IT ONLY WORKS WHEN THE VOLUME IS CENTERED IN ALL DIMENSIONS
	%warndlg(sprintf('%s: This function may not yet properly correct for the origin.\nCheck results in AFNI', FuncName));
	
	Info.DATASET_DIMENSIONS = Info_orig.DATASET_DIMENSIONS(maplocation);
	[err, Info.ORIENT_SPECIFIC] =  AFNI_OrientCode(NewCode);

	if (Opt.verbose),
		fprintf(1,'Writing resliced brick %s to disk ...', Info.RootName);
	end

	if (Opt.WriteBrick),
		[err, ErrMessage, Info] = WriteBrik (V, Info, Opt);
	end
	
	if (Opt.verbose),
		fprintf(1,'\tDone.\n');
	end

err = 0;
return;

