function [R] = McRetroTS(varargin)
% McRetroTS calls RetroTS to generate slice based regressors
% using respiration and pulse oximetry data.
% Compilable function
% Non-compiled version example (inside Matlab) :
%  McRetroTS(Opt.Respfile=resp4444.1D,Opt.Cardfile=card2222.1D,\
%            Opt.VolTR=2, Opt.Nslices=34, Opt.PhysFS=50, \
%            Opt.SliceOrder='alt-z', \
%            Opt.option1_name=option1_value,
%            Opt.option2_name=option2_value, ... )
% command line is the same except no parentheses or commas on line
%
% old arguments: (Respfile, Cardfile, VolTR, Nslices, PhysFS, ShowGraphs)
% Version of RetroTS made for Matlab compiler
%   No variables predefined and takes simple options on command line
%
% function [Opt, R, E] = RetroTS(SN)
%    [Opt, OptR, OptE] = RetroTS(Opt)
%This function creates slice-based regressors for regressing out
% components of heart rate, respiration and respiration volume per time.
%
%  Opt is the options structure with the following fields
%     Mandatory:
%     ----------
%     Respfile: Respiration data file
%     Cardfile: Cardiac data file
%     PhysFS: Physioliogical signal sampling frequency in Hz.
%     Nslices: Number of slices
%     VolTR: Volume TR in seconds
%     Optional:
%     ---------
%     Prefix: Prefix of output file
%     SliceOffset: Vector of slice acquisition time offsets.
%                  (default is equivalent of alt+z)
%     RVTshifts: Vector of shifts in seconds of RVT signal.
%                (default is [0:5:20])
%     RespCutoffFreq: Cut off frequency in Hz for respiratory lowpass filter
%                     (default 3 Hz)
%     CardCutoffFreq: Cut off frequency in Hz for cardiac lowpass filter
%                     (default 3 Hz)
%     ResamKernel: Resampling kernel.
%                 (default is 'linear', see help interp1 for more options)
%     FIROrder: Order of FIR filter. (default is 40)
%     Quiet: 1/0  flag. (defaut is 1)
%     Demo: 1/0 flag. (default is 0)
%
%Example:
%
%  Opt.Respfile = 'Resp_epiRT_scan_14.dat'
%  Opt.Cardfile = 'ECG_epiRT_scan_14.dat'
%  Opt.VolTR = 2
%  Opt.Nslices = 20;
%  Opt.PhysFS = 1./0.02;   %20 msec sampling period
%  RetroTS(Opt);
%

%Implementation Notes:
%%%%%%%%%%%%%%%%%%%%%%
% The script is intended as a prototype for development in C or Python
% The important routines are:
%    hilbert: Easily implemented with fft and ifft
%    interp1: A table lookup interpolation
%    fir: Tool for designing filters (we can just take it's coefficients)
%    filter: function to apply fir filter parameters (easy)
%
% All of the above can be easily implemented in C. However, I find it
% very useful to be able to plot the various steps in the process as we
% will undoubtedly face problems in the future. So I would vote for
% Python, assuming library vintage is not an issue. It looks like the
%

if (nargin < 1),
   fprintf(2,'McRetroTS - slice based regressor computation\n');
   fprintf(2,'Need some input.\n');
   fprintf(2,'Two usage modes: \n');
   fprintf(2,'\nUsage mode 1 (Provide Options via Opt structure\n');
   fprintf(2,'  as in Matlab with RetroTS.m)\n');
   fprintf(2,'\n   McRetroTS Opt.Respfile=Resp_epiRT_scan_14.dat \\\n');
   fprintf(2,'          Opt.Cardfile=ECG_epiRT_scan_14.dat \\\n');
   fprintf(2,'          Opt.VolTR=2 Opt.Nslices=34, Opt.PhysFS=50 \\\n');
   fprintf(2,'          Opt.SliceOrder=alt-z ....\n');
   fprintf(2,'\nUsage mode 2 (Provide 6 fixed parameters)\n');
   fprintf(2,'   This is the original implementation.\n');
   fprintf(2,'\nMcRetroTS Respfile Cardfile VolTR Nslices PhysFS Graph\n');
   fprintf(2,'\nExample: \n');
   fprintf(2,'  McRetroTS Resp_epiRT_scan_14.dat ECG_epiRT_scan_14.dat 2 20 50 0\n');
   fprintf(2,'The second method still allows for other options to be\n');
   fprintf(2,' passed on the command line in the method of the first\n');
   fprintf(2,' usage, but they must follow the first 6 parameters\n');
   fprintf(2,'See RetroTS.m for more details\n\n');
   if(~isdeployed)
      fprintf(2,'\nMcRetroTS help: ********************************\n');
      help('McRetroTS');
      fprintf(2,'\nRetroTS help: **********************************\n');
      help('RetroTS');
   else
      fprintf('More help is available by looking at the comments\n');
      fprintf('at the beginning of the source files: \n');
      fprintf('     McRetroTS.m and RetroTS.m\n');
      fprintf('Both available in the AFNI Matlab distribution and\n');
      fprintf(' on the AFNI website\n');
   end;
   R=1;
   return;
end

Opt.Quiet = 1;
Opt.Demo = 0;

% check if this is old school
%  get exactly 6 inputs that do not start with Opt.
if(sum(~strncmp(varargin, 'Opt.', 4)) == 6)
    num_in(1:6) = 0;
    Opt.Respfile = varargin{1};
    Opt.Cardfile = varargin{2};
    % convert the other arguments to numbers
    for i=3:1:6
         if(ischar(varargin{i}))
            num_in(i) = str2num(varargin{i});
         else
            num_in(i) = varargin{i};
         end
    end
    Opt.VolTR = num_in(3);
    Opt.Nslices = num_in(4);
    Opt.PhysFS = num_in(5);
    Opt.ShowGraphs = num_in(6);

%     Opt.Respfile = Respfile;
%     Opt.Cardfile = Cardfile;
%     if (ischar(VolTR))
%        Opt.VolTR=str2num(VolTR);
%     else
%        Opt.VolTR = VolTR
%     end
%     if (ischar(Nslices))
%        Opt.Nslices=str2num(Nslices);
%     else
%        Opt.Nslices = Nslices
%     end
%
%     if (ischar(PhysFS))
%        Opt.PhysFS = str2num(PhysFS);
%     else
%        Opt.PhysFS = PhysFS;
%     end
%
%     if (ischar(ShowGraphs))
%        Opt.ShowGraphs = str2num(ShowGraphs);
%     else
%        Opt.ShowGraphs = ShowGraphs
%     end
end
% process any other options not included (could be all options)
Opt = Args_to_opts(Opt,1,varargin);

% call the original
RetroTS(Opt);
if (Opt.ShowGraphs)
   user_input = input('Close all figures to exit\n');
   %uiwait(gcf);
end

R=0;
return;
