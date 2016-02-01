function [vertices, label, colortable] = Read_Brain_Annotation(filename, varargin)
%
% NAME
%
%       function [vertices, label, colortable] = ...
%                                       read_annotation(filename [, verbosity])
%
% ARGUMENTS
% INPUT
%       filename        string          name of annotation file to read
%
% OPTIONAL
%       verbosity       int             if true (>0), disp running output
%                                       + if false (==0), be quiet and do not
%                                       + display any running output
%
% OUTPUT
%       vertices        vector          vector with values running from 0 to
%                                       + size(vertices)-1
%       label           vector          lookup of annotation values for 
%                                       + corresponding vertex index.
%       colortable      struct          structure of annotation data
%                                       + see below
%       
% DESCRIPTION
%
%       This function essentially reads in a FreeSurfer annotation file
%       <filename> and returns structures and vectors that together 
%       assign each index in the surface vector to one of several 
%       structure names.
%       
% COLORTABLE STRUCTURE
% 
%       Consists of the following fields:
%       o numEntries:   number of entries
%       o orig_tab:     filename of original colortable file
%       o struct_names: cell array of structure names
%       o table:        n x 5 matrix
%                       Columns 1,2,3 are RGB values for struct color
%                       Column 4 is a flag (usually 0)
%                       Column 5 is the structure ID, calculated from
%                       R + G*2^8 + B*2^16 + flag*2^24
%                       
% LABEL VECTOR
% 
%       Each component of the <label> vector has a structureID value. To
%       match the structureID value with a structure name, lookup the row
%       index of the structureID in the 5th column of the colortable.table
%       matrix. Use this index as an offset into the struct_names field
%       to match the structureID with a string name.      
%
% PRECONDITIONS
%
%       o <filename> must be a valid FreeSurfer annotation file.
%       
% POSTCONDITIONS
%
%       o <colortable> will be an empty struct if not embedded in a
%         FreeSurfer annotation file. 
%       

%
% read_annotation.m
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision$
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

fp = fopen(filename, 'r', 'b');

verbosity = 1;
if length(varargin)
    verbosity       = varargin{1};  
end;

if(fp < 0)
   if verbosity, disp('Annotation file cannot be opened'); end;
   return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
   if verbosity, disp('No Colortable found.'); end;
   colortable = struct([]);
   fclose(fp);
   return; 
end

if(bool)
    
    %Read colortable
    numEntries = fread(fp, 1, 'int');

    if(numEntries > 0)
        
        if verbosity, disp(['Reading from Original Version']); end;
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);

        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        if verbosity
            disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        end
    else
        version = -numEntries;
        if verbosity
          if(version~=2)    
            disp(['Error! Does not handle version ' num2str(version)]);
          else
            disp(['Reading from version ' num2str(version)]);
          end
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
              if verbosity, disp(['Error! Read entry, index ' num2str(structure)]); end;
            end
            if(~isempty(colortable.struct_names{structure}))
              if verbosity, disp(['Error! Duplicate Structure ' num2str(structure)]); end;
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;       
        end
        if verbosity 
          disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        end
    end    
else
    if verbosity
        disp('Error! Should not be expecting bool = 0');    
    end;
end

fclose(fp);


