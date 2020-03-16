% This scritp allows to easily compile the MEX functions in the toolbox.
% Once the compiler has been setup properly (see help for function MEX in
% Matlab) simply run the script. Remember that, in order for the make file
% to work properly, the toolbox folders need to be included in the Matlab
% path (in the Set Path Matlab window (File>Set Path...) click the "Add 
% with Subfolders..." button, navigate to the toolbox folder, select it and
% save).

informationDir = which('information.m');
toolboxDir = informationDir(1:end-14);
% Verifying that we are compiling for the correct version
vdot_locations = strfind(toolboxDir, 'v.');
version = toolboxDir(vdot_locations(end):end);
if ~strcmp(version, 'v.1.0.5')
    error('Update the path to the correct toolbox version before proceeding with compilation.');
end

currentDir = cd;

disp('Compiling files for Information Breakdown ToolBox v.1.0.5.');
dirName = fullfile(toolboxDir, 'Extrapolations', 'Auxiliary Functions', 'partition_R');
cd(dirName);
mex -O -LargeArrayDims partition_R.c

dirName = fullfile(toolboxDir, 'Extrapolations', 'Auxiliary Functions', 'shuffle_R_across_trials');
cd(dirName);
mex -O -LargeArrayDims shuffle_R_across_trials.c

dirName = fullfile(toolboxDir, 'Methods', 'Auxiliary Functions', 'shuffle_R_across_cells');
cd(dirName);
mex -O -LargeArrayDims shuffle_R_across_cells.c

dirName = fullfile(toolboxDir, 'Methods', 'Direct Method');
cd(dirName);
mex -O -LargeArrayDims direct_method_v5b.c panzeri_and_treves_96.c

% Back to original dir:
cd(currentDir);

disp('Files compiled successfully.');