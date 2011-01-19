function x = hh_create_colormap(fname,threshold,option)
% hh_create_colormap - This function create the color label for the
% segmented image
%
% Usage: x = hh_create_colormap(fname,threshold,option)
% fname:     The output color map 
% threshold: The cut off threshold
% option:    Set the tranparency of the labels 
%
% $author: Copyright (C) by Hung Dang$
% $Log: hh_create_colormap.m$
%
% Revision 1.3 Sun Sep 26 22:51:50 MDT 2010, hungptit
% Add the usage information for this function.
% 
% Revision 1.2 Wed Aug  4 12:56:23 MDT 2010, hungptit
% Set the transparency of voxels with value less than threshold *
% number of color ranges
%
% Revision 1.1 Thu Jul 22 09:41:51 MDT 2010, hungptit
% Copy from the HHSIM code.
%

SCALE = 255;
x = SCALE * [[0,0,0];colormap(hot)];   
x = int32(x);  

% Main code
N = length(x);
fid = fopen(sprintf('%s.lbl',fname),'wt');
if fid
    % Print some useful information
    fprintf(fid,'%s\n','################################################');
    fprintf(fid,'%s\n','# ITK-SnAP Label Description File');
    fprintf(fid,'%s\n','# File format:');
    fprintf(fid,'%s\n','# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL');
    fprintf(fid,'%s\n','# Fields:');
    fprintf(fid,'%s\n','#    IDX:   Zero-based index');
    fprintf(fid,'%s\n','#    -R-:   Red color component (0..255)');
    fprintf(fid,'%s\n','#    -G-:   Green color component (0..255)');
    fprintf(fid,'%s\n','#    -B-:   Blue color component (0..255)');
    fprintf(fid,'%s\n','#    -A-:   Label transparency (0.00 .. 1.00)');
    fprintf(fid,'%s\n','#    VIS:   Label visibility (0 or 1)');
    fprintf(fid,'%s\n','#    IDX:   Label mesh visibility (0 or 1)');
    fprintf(fid,'%s\n','#  LABEL:   Label description');
    fprintf(fid,'%s\n','#  Created by Hung Dang');
    fprintf(fid,'%s\n','################################################');

    % Write clear label
    fprintf(fid,['%5d \t%5d %5d %5d \t%3f %3d %3d \t "Value ' ...
                 '%2d"\n'],0,0,0,0,0,0,0,0);
    
    % Write other visible label
    switch option
      case 1
        for n = 1:N
            if (n < N * threshold)
                fprintf(fid,['%5d \t%5d %5d %5d \t%3f %3d %3d \t ' ...
                             '"region_%2d"\n'],n,x(n,1),x(n,2),x(n,3),0,0,0,n);
            else
                fprintf(fid,['%5d \t%5d %5d %5d \t%3f %3d %3d \t ' ...
                             '"region_%2d"\n'],n,x(n,1),x(n,2),x(n,3),min(n/64,1),1,1,n);
            end
        end
      otherwise
        for n = 1:N
            if (n < N * threshold)
                fprintf(fid,['%5d \t%5d %5d %5d \t%3f %3d %3d \t ' ...
                             '"region_%2d"\n'],n,x(n,1),x(n,2),x(n,3),0,0,0,n);
            else
                fprintf(fid,['%5d \t%5d %5d %5d \t%3f %3d %3d \t ' ...
                             '"region_%2d"\n'],n,x(n,1),x(n,2),x(n,3),1,1,1,n);                
            end
        end
    end
    fclose(fid);
else
    error('Could not write the color map to a lebel file');
end
