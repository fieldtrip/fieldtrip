% CPD_MAKE compiles several CPD functions as well as FGT-related functions,
%          Withough FGT mex function, the Fast Gauss
%          Transform option will not be available (opt.fgt=0).

% Copyright (C) 2008 Andriy Myronenko (myron@csee.ogi.edu)
%
%     This file is part of the Coherent Point Drift (CPD) package.
%
%     The source code is provided under the terms of the GNU General Public License as published by
%     the Free Software Foundation version 2 of the License.
% 
%     CPD package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with CPD package; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

function cpd_make()

psave=pwd; 
p = mfilename('fullpath'); 
[pathstr, name, ext] = fileparts(p);

%%%%%%%%%%%%%%%%%%%% cpd_P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Compiling cpd_P-mex functions...');
disp('If this is the first time you use mex, it will ask you to choose a compiler.');
disp('Just choose the matlab default one (usually option #1).');
cd (pathstr); cd mex;

try
    mex cpd_P.c;
    mex cpd_Pappmex.c;
    mex cpd_Pcorrespondence.c;
    disp('Compilation of cpd_P mex functions is SUCCESSFUL.');
catch
   disp('Compilation of cpd_P mex functions failed. Try to run mex -setup to adjust your compiler.');
    
end

%%%%%%%%%%%%%%%%%%%%%% FGT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Compiling Fast Gauss Transform (FGT) related functions...');
disp('If the compilation fails, you can still perfectly run CPD only without FGT support.');
cd (pathstr); cd FGT;
disp('');

try
    mex fgt_model.c;
    mex fgt_predict.c;
    disp('Compilation of FGT mex functions is SUCCESSFUL.');
catch
   disp('Compilation of FGT mex functions failed. Try to run mex -setup to adjust your compiler.');
end


cd(psave);