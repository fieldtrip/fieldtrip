function matlab_install(SuiteSparse_path)
%  Matlab function to compile all the c-files to mex in the GPstuff toolbox.
%
%  Some of the sparse GP functionalities in the toolbox require 
%  SuiteSparse toolbox by Tim Davis:
%    http://www.cise.ufl.edu/research/sparse/SuiteSparse/current/SuiteSparse/
%
%  This package includes the SuiteSparse version 3.4. 
% 
%  * To install without SuiteSparse run matlab_install
%  * To install with SuiteSparse run matlab_install('SuiteSparseOn')
%
%   The function matlab_install compiles the mex-files and prints on
%   the screen, which directories should be added to Matlab paths. 
    
% Copyright (c) 2008-2012 Jarno Vanhatalo
    
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


    
    if nargin < 1
        SuiteSparse_path = [];
        fprintf('\n The path to the SuiteSparse package is not provided. \n')
        fprintf('\n Installing GPstuff without compactly supported covariance functions. \n')
        fprintf(' You are not able to use the following covariance functions:  \n')
        fprintf(' gpcf_ppcs0  \n gpcf_ppcs1 \n gpcf_ppcs2 \n gpcf_ppcs3 \n\n\n')
    elseif strcmp(SuiteSparse_path, 'SuiteSparseOn')
        cdir = pwd;        
        cd SuiteSparse
        SuiteSparse_path = [pwd '/'];
        
        % Compile SuiteSparse
        fprintf('Compiling SuiteSparse. This may take a while \n \n')
        paths = SuiteSparse_install(false);
        
        cd(cdir)
        fprintf('Compiling GPstuff. This may take a while \n \n')
    else 
        error('Unknown input argument. See help matlab_install for usage.')
    end
            
    % Go to diag/ and compile the mex-functions
    fprintf('\n Compiling files in diag.\n \n')
    cd('diag')
    diag_install
    cd('..')
        
    % Go to dist/ and compile the mex-functions        
    fprintf('\n Compiling files in dist.\n \n')
    cd('dist')
    dist_install    
    cd('..')
    
    % Go to gp/ and compile the mex-functions        
    fprintf('\n Compiling files in gp.\n \n')
    cd('gp')
    gp_install(SuiteSparse_path)    
    cd('..')       
        
    % Go to mc/ and compile the mex-functions
    fprintf('\n Compiling files in mc. \n \n')
    cd('mc')
    mc_install    
    cd('..')       
    
    PP = pwd;
    S{1} = [PP '/diag']; 
    S{2} = [PP '/dist']; 
    S{3} = [PP '/gp']; 
    S{4} = [PP '/mc']; 
    S{5} = [PP '/misc']; 
    S{6} = [PP '/optim']; 
    S{7} = [PP '/xunit']; 
    
    fprintf ('\n The following paths have been added.  You may wish to add them\n') ;
    fprintf ('permanently, using the MATLAB pathtool command or copying the below\n') ;
    fprintf ('lines to your startup.m file. \n\n');
    for i = 1:length(S)
       addpath(S{i}); 
       fprintf ('addpath %s\n', S{i}) ;
    end
   
    if nargin==1
        fprintf ('\n')
        for k = 1:length (paths)
            fprintf ('addpath %s\n', paths {k}) ;
        end
    end
end
