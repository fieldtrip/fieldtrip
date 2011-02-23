function fiff_write_ch_info(fid,ch)
%
% fiff_write_ch_info(fid,ch)
%
% Writes a channel information record to a fif file
%
%     fid           An open fif file descriptor
%     ch            The channel information structure to write
%
%     The type, cal, unit, and pos members are explained in Table 9.5
%     of the MNE manual
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.6  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.5  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.4  2006/04/13 17:47:50  msh
%   Make coil coordinate transformation or EEG electrode location out of the channel  info.
%
%   Revision 1.3  2006/04/12 10:29:02  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%   Revision 1.2  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%   Revision 1.1  2005/12/05 16:01:04  msh
%   Added an initial set of fiff writing routines.
%
%

me='MNE:fiff_write_ch_info';

if nargin ~= 2
    error(me,'Incorrect number of arguments');
end

FIFF_CH_INFO=203;
FIFFT_CH_INFO_STRUCT=30;
FIFFV_NEXT_SEQ=0;

%typedef struct _fiffChPosRec {
%  fiff_int_t   coil_type;          /*!< What kind of coil. */
%  fiff_float_t r0[3];              /*!< Coil coordinate system origin */
%  fiff_float_t ex[3];              /*!< Coil coordinate system x-axis unit vector */
%  fiff_float_t ey[3];              /*!< Coil coordinate system y-axis unit vector */
%  fiff_float_t ez[3];             /*!< Coil coordinate system z-axis unit vector */
%} fiffChPosRec,*fiffChPos;        /*!< Measurement channel position and coil type */


%typedef struct _fiffChInfoRec {
%  fiff_int_t    scanNo;        /*!< Scanning order # */
%  fiff_int_t    logNo;         /*!< Logical channel # */
%  fiff_int_t    kind;          /*!< Kind of channel */
%  fiff_float_t  range;         /*!< Voltmeter range (only applies to raw data ) */
%  fiff_float_t  cal;           /*!< Calibration from volts to... */
%  fiff_ch_pos_t chpos;         /*!< Channel location */
%  fiff_int_t    unit;          /*!< Unit of measurement */
%  fiff_int_t    unit_mul;      /*!< Unit multiplier exponent */
%  fiff_char_t   ch_name[16];   /*!< Descriptive name for the channel */
%} fiffChInfoRec,*fiffChInfo;   /*!< Description of one channel */

datasize=4*13 + 4*7 + 16;
count = fwrite(fid,int32(FIFF_CH_INFO),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_CH_INFO_STRUCT),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(datasize),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFV_NEXT_SEQ),'int32');
if count ~= 1
    error(me,'write failed');
end
%
%   Start writing fiffChInfoRec
%
count = fwrite(fid,int32(ch.scanno),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(ch.logno),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(ch.kind),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,single(ch.range),'single');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,single(ch.cal),'single');
if count ~= 1
    error(me,'write failed');
end
%
%   fiffChPosRec follows
%
count = fwrite(fid,int32(ch.coil_type),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,single(ch.loc),'single');
if count ~= 12
    error(me,'write failed');
end
%
%   unit and unit multiplier
%
count = fwrite(fid,int32(ch.unit),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(ch.unit_mul),'int32');
if count ~= 1
    error(me,'write failed');
end
%
%   Finally channel name
%
len=length(ch.ch_name);
if len > 15
    ch_name = ch.ch_name(1:15);
else
    ch_name = ch.ch_name;
end
len = length(ch_name);
count = fwrite(fid,ch_name,'char');
if count ~= len
    error(me,'write failed');
end
if len < 16
    dum=zeros(1,16-len);
    count = fwrite(fid,uint8(dum),'uchar');
    if count ~= 16-len
        error(me,'write failed');
    end
end
return;
