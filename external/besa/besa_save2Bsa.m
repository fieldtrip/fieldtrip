function status = besa_save2Bsa(file_path, file_name, vertex_coords)
% BESA_SAVE2BSA saves a matrix with vertices to a BESA solution file 
% (.bsa). 
%
% Parameters:
%     [file_path]
%         Full path to the folder where the file should be saved.
% 
%     [file_name]
%         The name of the file where the output should be written.
% 
%     [data_matrix]
%         A 2D data matrix with vertex coordinates in Talairach coordinate 
%         system. It should have the size [NumberOfVertices x 3]. 
% 
% Return:
%     [status] 
%         The status of the writing process. If status < 0 then the process
%         wasn't successful.
% 
% Copyright (C) 2014, BESA GmbH
%
% File name: besa_save2Bsa.m
%
% Author: Todor Jordanov
% Created: 2014-03-11

status = 1;

% Check the size of the grid
Size1 = size(vertex_coords, 1);
Size2 = size(vertex_coords, 2);

if(Size1 == 3 && Size2 ~= 3)
    
    vertex_coords = vertex_coords';
    
end

NumVertices = size(vertex_coords, 1);

FullPath = fullfile(file_path, file_name);
% Open file for writing.
fid = fopen(FullPath, 'w');

% MATLAB reserves file identifiers 0, 1, and 2 for standard input,  
% standard output (the screen), and standard error, respectively. When 
% fopen successfully opens a file, it returns a file identifier greater 
% than or equal to 3.
if(fid >= 3)

    % Write the first line of the header.
    fprintf(fid, ...
        ['BSA_1.04.19990715|TC\n'...
        'Type          x-loc       y-loc       z-loc        x-ori       y-ori       z-ori      color state size\n'...
        '=======================================================================================================\n']);

    for i = 1:NumVertices

        fprintf(fid, 'RegSrc   %g  %g  %g     0.000000    0.000000    0.000000        255     2 4.50\n', ...
            vertex_coords(i, 1), vertex_coords(i, 2), vertex_coords(i, 3));

    end

    fclose(fid);

else
    
    status = -1;
    disp('Error! Invalid file identifier.')
    
end

