function [hdr] = read_ns_hdr(filename)

% READ_NS_HDR read the header from a NeuroScan 3.x or 4.x AVG/EEG/CNT File
%
% [hdr] = read_ns_hdr(filename)
%
% The output data structure hdr has the fields:
%   hdr.label       - electrode labels
%   hdr.nchan       - number of channels
%   hdr.npnt        - number of samplepoints in ERP waveform
%   hdr.rate        - sample rate (Hz)
%   hdr.xmin        - prestimulus epoch start (e.g., -100 msec)
%   hdr.xmax        - poststimulus epoch end (e.g., 900 msec)
%   hdr.nsweeps     - number of accepted trials/sweeps
%   hdr.domain      - time (0) or frequency (1) domain

% Copyright (C) 2002, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

fid = fopen_or_error(filename,'r','ieee-le');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these structures are from Neuroscan sethead.h
% the first two columns are the size in bytes and the offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERAL SETUP STRUCTURE, total 900 bytes
%    0     20  |  char   rev[20];         /* Revision string                         */
%   20      1  |  char   type;            /* File type AVG=1, EEG=0                  */
%   21     20  |  char   id[20];          /* Patient ID                              */
%   41     20  |  char   oper[20];        /* Operator ID                             */
%   61     20  |  char   doctor[20];      /* Doctor ID                               */
%   81     20  |  char   referral[20];    /* Referral ID                             */
%  101     20  |  char   hospital[20];    /* Hospital ID                             */
%  121     20  |  char   patient[20];     /* Patient name                            */
%  141      2  |  short  int age;         /* Patient Age                             */
%  143      1  |  char   sex;             /* Patient Sex Male='M', Female='F'        */
%  144      1  |  char   hand;            /* Handedness Mixed='M',Rt='R', lft='L'    */
%  145     20  |  char   med[20];         /* Medications                             */
%  165     20  |  char   class[20];       /* Classification                          */
%  185     20  |  char   state[20];       /* Patient wakefulness                     */
%  205     20  |  char   label[20];       /* Session label                           */
%  225     10  |  char   date[10];        /* Session date string                     */
%  235     12  |  char   time[12];        /* Session time strin                      */
%  247      4  |  float  mean_age;        /* Mean age (Group files only)             */
%  251      4  |  float  stdev;           /* Std dev of age (Group files only)       */
%  255      2  |  short int n;            /* Number in group file                    */
%  257     38  |  char   compfile[38];    /* Path and name of comparison file        */
%  295      4  |  float  SpectWinComp;    // Spectral window compensation factor
%  299      4  |  float  MeanAccuracy;    // Average respose accuracy
%  303      4  |  float  MeanLatency;     // Average response latency
%  307     46  |  char   sortfile[46];    /* Path and name of sort file              */
%  353      4  |  int    NumEvents;       // Number of events in eventable
%  357      1  |  char   compoper;        /* Operation used in comparison            */
%  358      1  |  char   avgmode;         /* Set during online averaging             */
%  359      1  |  char   review;          /* Set during review of EEG data           */
%  360      2  |  short unsigned nsweeps;      /* Number of expected sweeps               */
%  362      2  |  short unsigned compsweeps;   /* Number of actual sweeps                 */
%  364      2  |  short unsigned acceptcnt;    /* Number of accepted sweeps               */
%  366      2  |  short unsigned rejectcnt;    /* Number of rejected sweeps               */
%  368      2  |  short unsigned pnts;         /* Number of points per waveform           */
%  370      2  |  short unsigned nchannels;    /* Number of active channels               */
%  372      2  |  short unsigned avgupdate;    /* Frequency of average update             */
%  374      1  |  char  domain;           /* Acquisition domain TIME=0, FREQ=1       */
%  375      1  |  char  variance;         /* Variance data included flag             */
%  376      2  |  unsigned short rate;    /* D-to-A rate                             */
%  378      8  |  double scale;           /* scale factor for calibration            */
%  386      1  |  char  veogcorrect;      /* VEOG corrected flag                     */
%  387      1  |  char  heogcorrect;      /* HEOG corrected flag                     */
%  388      1  |  char  aux1correct;      /* AUX1 corrected flag                     */
%  389      1  |  char  aux2correct;      /* AUX2 corrected flag                     */
%  390      4  |  float veogtrig;         /* VEOG trigger percentage                 */
%  394      4  |  float heogtrig;         /* HEOG trigger percentage                 */
%  398      4  |  float aux1trig;         /* AUX1 trigger percentage                 */
%  402      4  |  float aux2trig;         /* AUX2 trigger percentage                 */
%  406      2  |  short int heogchnl;     /* HEOG channel number                     */
%  408      2  |  short int veogchnl;     /* VEOG channel number                     */
%  410      2  |  short int aux1chnl;     /* AUX1 channel number                     */
%  412      2  |  short int aux2chnl;     /* AUX2 channel number                     */
%  414      1  |  char  veogdir;          /* VEOG trigger direction flag             */
%  415      1  |  char  heogdir;          /* HEOG trigger direction flag             */
%  416      1  |  char  aux1dir;          /* AUX1 trigger direction flag             */
%  417      1  |  char  aux2dir;          /* AUX2 trigger direction flag             */
%  418      2  |  short int veog_n;       /* Number of points per VEOG waveform      */
%  420      2  |  short int heog_n;       /* Number of points per HEOG waveform      */
%  422      2  |  short int aux1_n;       /* Number of points per AUX1 waveform      */
%  424      2  |  short int aux2_n;       /* Number of points per AUX2 waveform      */
%  426      2  |  short int veogmaxcnt;   /* Number of observations per point - VEOG */
%  428      2  |  short int heogmaxcnt;   /* Number of observations per point - HEOG */
%  430      2  |  short int aux1maxcnt;   /* Number of observations per point - AUX1 */
%  432      2  |  short int aux2maxcnt;   /* Number of observations per point - AUX2 */
%  434      1  |  char   veogmethod;      /* Method used to correct VEOG             */
%  435      1  |  char   heogmethod;      /* Method used to correct HEOG             */
%  436      1  |  char   aux1method;      /* Method used to correct AUX1             */
%  437      1  |  char   aux2method;      /* Method used to correct AUX2             */
%  438      4  |  float  AmpSensitivity;  /* External Amplifier gain                 */
%  442      1  |  char   LowPass;         /* Toggle for Amp Low pass filter          */
%  443      1  |  char   HighPass;        /* Toggle for Amp High pass filter         */
%  444      1  |  char   Notch;           /* Toggle for Amp Notch state              */
%  445      1  |  char   AutoClipAdd;     /* AutoAdd on clip                         */
%  446      1  |  char   baseline;        /* Baseline correct flag                   */
%  447      4  |  float  offstart;        /* Start point for baseline correction     */
%  451      4  |  float  offstop;         /* Stop point for baseline correction      */
%  455      1  |  char   reject;          /* Auto reject flag                        */
%  456      4  |  float  rejstart;        /* Auto reject start point                 */
%  460      4  |  float  rejstop;         /* Auto reject stop point                  */
%  464      4  |  float  rejmin;          /* Auto reject minimum value               */
%  468      4  |  float  rejmax;          /* Auto reject maximum value               */
%  472      1  |  char   trigtype;        /* Trigger type                            */
%  473      4  |  float  trigval;         /* Trigger value                           */
%  477      1  |  char   trigchnl;        /* Trigger channel                         */
%  478      2  |  short int trigmask;     /* Wait value for LPT port                 */
%  480      4  |  float trigisi;          /* Interstimulus interval (INT trigger)    */
%  484      4  |  float trigmin;          /* Min trigger out voltage (start of pulse)*/
%  488      4  |  float trigmax;          /* Max trigger out voltage (during pulse)  */
%  492      1  |  char  trigdir;          /* Duration of trigger out pulse           */
%  493      1  |  char  Autoscale;        /* Autoscale on average                    */
%  494      2  |  short int n2;           /* Number in group 2 (MANOVA)              */
%  496      1  |  char  dir;              /* Negative display up or down             */
%  497      4  |  float dispmin;          /* Display minimum (Yaxis)                 */
%  501      4  |  float dispmax;          /* Display maximum (Yaxis)                 */
%  505      4  |  float xmin;             /* X axis minimum (epoch start in sec)     */
%  509      4  |  float xmax;             /* X axis maximum (epoch stop in sec)      */
%  513      4  |  float AutoMin;          /* Autoscale minimum                       */
%  517      4  |  float AutoMax;          /* Autoscale maximum                       */
%  521      4  |  float zmin;             /* Z axis minimum - Not currently used     */
%  525      4  |  float zmax;             /* Z axis maximum - Not currently used     */
%  529      4  |  float lowcut;           /* Archival value - low cut on external amp*/
%  533      4  |  float highcut;          /* Archival value - Hi cut on external amp */
%  537      1  |  char  common;           /* Common mode rejection flag              */
%  538      1  |  char  savemode;         /* Save mode EEG AVG or BOTH               */
%  539      1  |  char  manmode;          /* Manual rejection of incomming data      */
%  540     10  |  char  ref[10];          /* Label for reference electode            */
%  550      1  |  char  Rectify;          /* Rectification on external channel       */
%  551      4  |  float DisplayXmin;      /* Minimun for X-axis display              */
%  555      4  |  float DisplayXmax;      /* Maximum for X-axis display              */
%  559      1  |  char  phase;            /* flag for phase computation              */
%  560     16  |  char  screen[16];       /* Screen overlay path name                */
%  576      2  |  short int CalMode;      /* Calibration mode                        */
%  578      2  |  short int CalMethod;    /* Calibration method                      */
%  580      2  |  short int CalUpdate;    /* Calibration update rate                 */
%  582      2  |  short int CalBaseline;  /* Baseline correction during cal          */
%  584      2  |  short int CalSweeps;    /* Number of calibration sweeps            */
%  586      4  |  float CalAttenuator;    /* Attenuator value for calibration        */
%  590      4  |  float CalPulseVolt;     /* Voltage for calibration pulse           */
%  594      4  |  float CalPulseStart;    /* Start time for pulse                    */
%  598      4  |  float CalPulseStop;     /* Stop time for pulse                     */
%  602      4  |  float CalFreq;          /* Sweep frequency                         */
%  606     34  |  char  taskfile[34];     /* Task file name                          */
%  640     34  |  char  seqfile[34];      /* Sequence file path name                 */
%  674      1  |  char  SpectMethod;      // Spectral method
%  675      1  |  char  SpectScaling;     // Scaling employed
%  676      1  |  char  SpectWindow;      // Window employed
%  677      4  |  float SpectWinLength;   // Length of window %
%  681      1  |  char  SpectOrder;       // Order of Filter for Max Entropy method
%  682      1  |  char  NotchFilter;      // Notch Filter in or out
%  683     11  |  char  unused[11];       // Free space
%  694      2  |  short  FspStopMethod;   /* FSP - Stoping mode                      */
%  696      2  |  short  FspStopMode;     /* FSP - Stoping mode                      */
%  698      4  |  float FspFValue;        /* FSP - F value to stop terminate         */
%  702      2  |  short int FspPoint;     /* FSP - Single point location             */
%  704      2  |  short int FspBlockSize; /* FSP - block size for averaging          */
%  706      2  |  unsigned short FspP1;   /* FSP - Start of window                   */
%  708      2  |  unsigned short FspP2;   /* FSP - Stop  of window                   */
%  710      4  |  float FspAlpha;         /* FSP - Alpha value                       */
%  714      4  |  float FspNoise;         /* FSP - Signal to ratio value             */
%  718      2  |  short int FspV1;        /* FSP - degrees of freedom                */
%  720     40  |  char  montage[40];      /* Montage file path name                  */
%  760     40  |  char  EventFile[40];    /* Event file path name                    */
%  800      4  |  float fratio;           /* Correction factor for spectral array    */
%  804      1  |  char  minor_rev;        /* Current minor revision                  */
%  805      2  |  short int eegupdate;    /* How often incomming eeg is refreshed    */
%  807      1  |  char   compressed;      /* Data compression flag                   */
%  808      4  |  float  xscale;          /* X position for scale box - Not used     */
%  812      4  |  float  yscale;          /* Y position for scale box - Not used     */
%  816      4  |  float  xsize;           /* Waveform size X direction               */
%  820      4  |  float  ysize;           /* Waveform size Y direction               */
%  824      1  |  char   ACmode;          /* Set SYNAP into AC mode                  */
%  825      1  |  unsigned char   CommonChnl;      /* Channel for common waveform             */
%  826      1  |  char   Xtics;           /* Scale tool- 'tic' flag in X direction   */
%  827      1  |  char   Xrange;          /* Scale tool- range (ms,sec,Hz) flag X dir*/
%  828      1  |  char   Ytics;           /* Scale tool- 'tic' flag in Y direction   */
%  829      1  |  char   Yrange;          /* Scale tool- range (uV, V) flag Y dir    */
%  830      4  |  float  XScaleValue;     /* Scale tool- value for X dir             */
%  834      4  |  float  XScaleInterval;  /* Scale tool- interval between tics X dir */
%  838      4  |  float  YScaleValue;     /* Scale tool- value for Y dir             */
%  842      4  |  float  YScaleInterval;  /* Scale tool- interval between tics Y dir */
%  846      4  |  float  ScaleToolX1;     /* Scale tool- upper left hand screen pos  */
%  850      4  |  float  ScaleToolY1;     /* Scale tool- upper left hand screen pos  */
%  854      4  |  float  ScaleToolX2;     /* Scale tool- lower right hand screen pos */
%  858      4  |  float  ScaleToolY2;     /* Scale tool- lower right hand screen pos */
%  862      2  |  short int port;         /* Port address for external triggering    */
%  864      4  |  long  NumSamples;       /* Number of samples in continuous file    */
%  868      1  |  char  FilterFlag;       /* Indicates that file has been filtered   */
%  869      4  |  float LowCutoff;        /* Low frequency cutoff                    */
%  873      2  |  short int LowPoles;     /* Number of poles                         */
%  875      4  |  float HighCutoff;       /* High frequency cutoff                   */
%  879      2  |  short int HighPoles;    /* High cutoff number of poles             */
%  881      1  |  char  FilterType;       /* Bandpass=0 Notch=1 Highpass=2 Lowpass=3 */
%  882      1  |  char  FilterDomain;     /* Frequency=0 Time=1                      */
%  883      1  |  char  SnrFlag;          /* SNR computation flag                    */
%  884      1  |  char  CoherenceFlag;    /* Coherence has been computed             */
%  885      1  |  char  ContinuousType;   /* Method used to capture events in *.cnt  */
%  886      4  |  long  EventTablePos;    /* Position of event table                 */
%  890      4  |  float ContinousSeconds; // Number of seconds to displayed per page
%  894      4  |  long  ChannelOffset;    // Block size of one channel in SYNAMPS
%  898      1  |  char  AutoCorrectFlag;  // Autocorrect of DC values
%  899      1  |  unsigned char DCThreshold; // Auto correct of DC level

% ELECTRODE STRUCTURE, total 75 bytes
%   10     0   |  char  lab[10];          /* Electrode label - last bye contains NULL */
%    1    10   |  char  reference;        /* Reference electrode number               */
%    1    11   |  char  skip;             /* Skip electrode flag ON=1 OFF=0           */
%    1    12   |  char  reject;           /* Artifact reject flag                     */
%    1    13   |  char  display;          /* Display flag for 'STACK' display         */
%    1    14   |  char  bad;              /* Bad electrode flag                       */
%    2    15   |  unsigned short int n;   /* Number of observations                   */
%    1    17   |  char  avg_reference;    /* Average reference status                 */
%    1    18   |  char  ClipAdd;          /* Automatically add to clipboard           */
%    4    19   |  float x_coord;          /* X screen coord. for 'TOP' display        */
%    4    23   |  float y_coord;          /* Y screen coord. for 'TOP' display        */
%    4    27   |  float veog_wt;          /* VEOG correction weight                   */
%    4    31   |  float veog_std;         /* VEOG std dev. for weight                 */
%    4    35   |  float snr;              /* signal-to-noise statistic                */
%    4    39   |  float heog_wt;          /* HEOG Correction weight                   */
%    4    43   |  float heog_std;         /* HEOG Std dev. for weight                 */
%    2    47   |  short int baseline;     /* Baseline correction value in raw ad units*/
%    1    49   |  char  Filtered;         /* Toggel indicating file has be filtered   */
%    1    50   |  char  Fsp;              /* Extra data                               */
%    4    51   |  float aux1_wt;          /* AUX1 Correction weight                   */
%    4    54   |  float aux1_std;         /* AUX1 Std dev. for weight                 */
%    4    59   |  float sensitivity;      /* electrode sensitivity                    */
%    1    63   |  char  Gain;             /* Amplifier gain                           */
%    1    64   |  char  HiPass;           /* Hi Pass value                            */
%    1    65   |  char  LoPass;           /* Lo Pass value                            */
%    1    66   |  unsigned char Page;     /* Display page                             */
%    1    67   |  unsigned char Size;     /* Electrode window display size            */
%    1    68   |  unsigned char Impedance;/* Impedance test                           */
%    1    69   |  unsigned char PhysicalChnl; /* Physical channel used                    */
%    1    70   |  char  Rectify;          /* Free space                               */
%    4    71   |  float calib;            /* Calibration factor                       */

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read general parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
offset_type       =  20;
offset_nsweeps    = 364;
offset_pnts       = 368;
offset_nchans     = 370;
offset_domain     = 374;
offset_variance   = 375;
offset_rate       = 376;
offset_xmin       = 505;
offset_xmax       = 509;
fseek(fid, offset_type,     'bof');        type      = fread(fid, 1, 'ushort');
fseek(fid, offset_nsweeps,  'bof');    hdr.nsweeps   = fread(fid, 1, 'ushort');
fseek(fid, offset_pnts,     'bof');    hdr.npnt      = fread(fid, 1, 'ushort');
fseek(fid, offset_nchans,   'bof');    hdr.nchan     = fread(fid, 1, 'ushort');
fseek(fid, offset_domain,   'bof');    hdr.domain    = fread(fid, 1, 'uchar');
fseek(fid, offset_variance, 'bof');    hdr.variance  = fread(fid, 1, 'uchar');
fseek(fid, offset_rate,     'bof');    hdr.rate      = fread(fid, 1, 'ushort');
fseek(fid, offset_xmin,     'bof');    hdr.xmin      = fread(fid, 1, 'float32') * 1000;
fseek(fid, offset_xmax,     'bof');    hdr.xmax      = fread(fid, 1, 'float32') * 1000;
fseek(fid, offset_xmax,     'bof');    hdr.xmax      = fread(fid, 1, 'float32') * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read electrode configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdr.label = {};
for elec = 1:hdr.nchan,
  offset_label    = 900+(elec-1)*75 + 0;
  offset_baseline = 900+(elec-1)*75 + 47;
  offset_sens     = 900+(elec-1)*75 + 59;
  offset_calib    = 900+(elec-1)*75 + 71;
  fseek(fid, offset_label,    'bof'); lab                   = fread(fid, 10, 'uchar'); lab(find(lab==0)) = ' ';
  fseek(fid, offset_baseline, 'bof'); hdr.baseline(elec)    = fread(fid, 1, 'ushort');
  fseek(fid, offset_sens,     'bof'); hdr.sensitivity(elec) = fread(fid, 1, 'float32');
  fseek(fid, offset_calib,    'bof'); hdr.calib(elec)       = fread(fid, 1, 'float32');
  hdr.label{elec}  = fliplr(deblank(fliplr(deblank(char(lab')))));
  hdr.factor(elec) = hdr.calib(elec) * hdr.sensitivity(elec) / 204.8;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to infer whether time or frequency data is represented
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hdr.domain==0
  % data in time domain
elseif hdr.domain==1
  % data in frequency domain
else
  % probably old datafile, assume the data to be in time domain
  ft_warning('assuming the data to be in time domain (domain was %d)', hdr.domain);
  hdr.domain=0;
end

fclose(fid);
