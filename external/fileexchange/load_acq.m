%LOAD_ACQ  load BIOPAC's AcqKnowledge file format for Windows PC
%
%  Usage: acq = load_acq(filename, [force_chan_by_chan_load])
%
%  acq - AcqKnowledge file structure containing ACQ header field, 
%	and data matrix.
%
%  filename - BIOPAC's AcqKnowledge file
%
%  force_chan_by_chan_load - Optional. By default, this optional flag will
%	be set to 1, which means that when you use acq = load_acq(fname),
%	the data will be loaded one channel after the other. This can
%	avoid memory crash when you load very large ACQ data. If your
%	ACQ data is not huge, I suggest that you set this optional flag to
%	0, i.e. acq = load_acq(fname, 0). In this case, the program will
%	read data depending on the data type. If the program detects that
%	the data type in ACQ file are different from channel to channel,
%	it will still read data channel by channel. Otherwise, it will
%	read whole data in one block (a lot faster than using traditional
%	way from channel to channel with the same result).
%
%  Notes:
%
%  This program is based on Application Note #156 from BIOPAC web site:
%  http://www.biopac.com/ResearchNotes.asp?Aid=&ANid=82&Level=4
%  ( 156 - ACQKNOWLEDGE FILE FORMATS FOR PC WITH WINDOWS )
%
%  The note mentioned that: "This document describes file formatting for
%  all Windows versions of AcqKnowledge 3.9.x or below". Thanks to the 
%  open Python source code provided by Nathan Vack, this program can also
%  read AcqKnowledge 4.0 & 4.1 data (with no documentation from BIOPAC).
%  Compressed data is not supported by this program.
%
%  Created on 5-APR-2007 by Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
%  Modify on 31-JUL-2007 by Julio Cruz (julio.cruz@juliocruz.info) to
%	make the program also work on MAC OS.
%
%  Modify on 01-MAY-2008 by Antonio Molins (amolins@mit.edu) to
%	make the program work with files recorded in Biopac student data 
%	acquisition system. Changes include:
%	- added reading of markers, with "markers" field included in 
%	returned structure..
%	- proper handling of foreign data (was not getting quite there with
%	the code (as it was downloaded from MathWorks).
%	- proper handling of var_sampling_rate = 0, as produced by the
%	Biopac software in occasions.
%
function acq = load_acq(filename, chan_by_chan)

   if ~exist('chan_by_chan','var')
      chan_by_chan = 1;
   end

   if ~exist('filename','var')
      error('Usage: acq = load_acq(filename)');
   end

   if ~exist(filename,'file')
      error([filename, ': Can''t open file']);
   end

   fid = fopen(filename,'r', 'l');

   if fid < 0,
      msg = sprintf('Cannot open file %s.',filename);
      error(msg);
   end

   fread(fid, 1, 'int16')';
   file_version = fread(fid, 1, 'int32')';

   %  try different endian
   %
   if file_version < 0 | file_version > 200
      fclose(fid);
      fid = fopen(filename,'r', 'b');

      fread(fid, 1, 'int16')';
      file_version = fread(fid, 1, 'int32')';
      fseek(fid, 0, 'bof');

      if file_version < 0 | file_version > 200
         error('This ACQ file is not supported');
      end
   end

   fprintf('Loading %s ', filename);

   %  read header
   %
   acq.hdr = read_acq_hdr(fid);

   %  read data
   %
   acq.data = read_acq(fid, acq.hdr, chan_by_chan);

   %  read markers written by AM (Antonio Molins) for ACQ 3
   %
   if file_version < 68
      acq.markers = read_markers(fid, acq.hdr, size(acq.data));
   end

   fclose(fid);

   fprintf(' Done!\n');

   return;					% load_acq


%---------------------------------------------------------------------
function hdr = read_acq_hdr(fid)

   fseek(fid, 0, 'bof');
   hdr.graph = graph(fid);
   acc_chan_header_len = 0;

   for i = 1:hdr.graph.num_channels
      fseek(fid, hdr.graph.ext_item_header_len+acc_chan_header_len, 'bof');
      hdr.per_chan_data(i) = per_chan_data(fid, hdr.graph.file_version);
      acc_chan_header_len = acc_chan_header_len + hdr.per_chan_data(i).chan_header_len;
   end

   % added by AM: foreign data section was being started to read one short,
   % this fseek takes care of that
   fseek(fid, hdr.graph.ext_item_header_len+acc_chan_header_len, 'bof');
      
   hdr.foreign = foreign(fid, hdr.graph.file_version);

   if hdr.graph.file_version >= 68 & hdr.graph.file_version < 80
      unused = fread(fid, hdr.foreign.length2-12, 'uint8')';
   end

   for i = 1:hdr.graph.num_channels
      hdr.per_chan_type(i) = per_chan_type(fid);
   end

   return;					% read_acq_hdr


%---------------------------------------------------------------------
function hdr = graph(fid, file_version)

   %  Struct						% off + size
   unused = fread(fid, 1, 'int16')';			% 0 + 2
   hdr.file_version = fread(fid, 1, 'int32');		% 2 + 4
   hdr.ext_item_header_len = fread(fid, 1, 'int32');	% 6 + 4
   hdr.num_channels = fread(fid, 1, 'int16');		% 10 + 2
   hdr.horiz_axis_type = fread(fid, 1, 'int16');	% 12 + 2
   hdr.curr_channel = fread(fid, 1, 'int16');		% 14 + 2
   hdr.sample_time = fread(fid, 1, 'double');		% 16 + 8
   hdr.time_offset = fread(fid, 1, 'double');		% 24 + 8
   hdr.time_scale = fread(fid, 1, 'double');		% 32 + 8
   hdr.time_cursor1 = fread(fid, 1, 'double');		% 40 + 8
   hdr.time_cursor2 = fread(fid, 1, 'double');		% 48 + 8
   hdr.chart_window = fread(fid, 4, 'int16');		% 56 + 8
   hdr.measurement = fread(fid, 6, 'int16');		% 64 + 12
   hdr.hilite = fread(fid, 1, 'int16');			% 76 + 2
   hdr.first_time_offset = fread(fid, 1, 'double');	% 78 + 8
   hdr.rescale = fread(fid, 1, 'int16');		% 86 + 2
   hdr.horiz_units1 = deblank(fread(fid, 40, '*char')'); % 88 + 40
   hdr.horiz_units2 = deblank(fread(fid, 10, '*char')'); % 128 + 10
   hdr.in_memory = fread(fid, 1, 'int16');		% 138 + 2
   hdr.grid = fread(fid, 1, 'int16');			% 140 + 2
   hdr.markers = fread(fid, 1, 'int16');		% 142 + 2
   hdr.plot_draft = fread(fid, 1, 'int16');		% 144 + 2
   hdr.display_mode = fread(fid, 1, 'int16');		% 146 + 2
   hdr.reserved = fread(fid, 1, 'int16');		% 148 + 2

   if hdr.file_version > 33 & hdr.file_version < 68
      hdr.show_toolbar = fread(fid, 1, 'int16');	% 150 + 2
      hdr.show_chan_butt = fread(fid, 1, 'int16');	% 152 + 2
      hdr.show_measurement = fread(fid, 1, 'int16');	% 154 + 2
      hdr.show_marker = fread(fid, 1, 'int16');		% 156 + 2
      hdr.show_journal = fread(fid, 1, 'int16');	% 158 + 2
      hdr.cur_x_channel = fread(fid, 1, 'int16');	% 160 + 2
      hdr.mmt_precision = fread(fid, 1, 'int16');	% 162 + 2
   end

   if hdr.file_version > 34 & hdr.file_version < 68
      hdr.measurement_row = fread(fid, 1, 'int16');	% 164 + 2
      hdr.mmt = fread(fid, 40, 'int16');		% 166 + 80
      hdr.mmt_chan = fread(fid, 40, 'int16');		% 246 + 80
   end

   if hdr.file_version > 35 & hdr.file_version < 68
      hdr.mmt_calc_opnd1 = fread(fid, 40, 'int16');	% 326 + 80
      hdr.mmt_calc_opnd2 = fread(fid, 40, 'int16');	% 406 + 80
      hdr.mmt_calc_op = fread(fid, 40, 'int16');	% 486 + 80
      hdr.mmt_calc_constant = fread(fid, 40, 'double');	% 566 + 320
   end

   if hdr.file_version > 37 & hdr.file_version < 68
      hdr.new_grid_minor = fread(fid, 1, 'int32');	% 886 + 4
      hdr.color_major_grid = fread(fid, 1, 'int32');	% 890 + 4
      hdr.color_minor_grid = fread(fid, 1, 'int32');	% 894 + 4
      hdr.major_grid_style = fread(fid, 1, 'int16');	% 898 + 2
      hdr.minor_grid_style = fread(fid, 1, 'int16');	% 900 + 2
      hdr.major_grid_width = fread(fid, 1, 'int16');	% 902 + 2
      hdr.minor_grid_width = fread(fid, 1, 'int16');	% 904 + 2
      hdr.fixed_units_div = fread(fid, 1, 'int32');	% 906 + 4
      hdr.mid_range_show = fread(fid, 1, 'int32');	% 910 + 4
      hdr.start_middle_point = fread(fid, 1, 'double');	% 914 + 8
      hdr.offset_point = fread(fid, 60, 'double');	% 922 + 480
      hdr.h_grid = fread(fid, 1, 'double');		% 1402 + 8
      hdr.v_grid = fread(fid, 60, 'double');		% 1410 + 480
      hdr.enable_wave_tool = fread(fid, 1, 'int32');	% 1890 + 4

      %  interpret color_major_grid
      %
      color_str = sprintf('%06s', dec2hex(hdr.color_major_grid));
      hdr.color_major_grid = ...
         [hex2dec(color_str(5:6)) hex2dec(color_str(3:4)) hex2dec(color_str(1:2))]/255;

      %  interpret color_minor_grid
      %
      color_str = sprintf('%06s', dec2hex(hdr.color_minor_grid));
      hdr.color_minor_grid = ...
         [hex2dec(color_str(5:6)) hex2dec(color_str(3:4)) hex2dec(color_str(1:2))]/255;

   end

   if hdr.file_version > 38 & hdr.file_version < 68
      hdr.horiz_precision = fread(fid, 1, 'int16');	% 1894 + 2
   end

   if hdr.file_version > 40 & hdr.file_version < 68
      hdr.reserved2 = fread(fid, 20, 'int8');		% 1896 + 20
      hdr.overlap_mode = fread(fid, 1, 'int32');	% 1916 + 4
      hdr.show_hardware = fread(fid, 1, 'int32');	% 1920 + 4
      hdr.x_auto_plot = fread(fid, 1, 'int32');	% 1924 + 4
      hdr.x_auto_scroll = fread(fid, 1, 'int32');	% 1928 + 4
      hdr.start_butt_visible = fread(fid, 1, 'int32');	% 1932 + 4
      hdr.compressed = fread(fid, 1, 'int32');		% 1936 + 4
      hdr.always_start_butt_visible = fread(fid, 1, 'int32');	% 1940 + 4
   end

   if hdr.file_version > 42 & hdr.file_version < 68
      hdr.path_video = deblank(fread(fid, 260, '*char')'); % 1944 + 260
      hdr.opt_sync_delay = fread(fid, 1, 'int32');	% 2204 + 4
      hdr.sync_delay = fread(fid, 1, 'double');		% 2208 + 8
      hdr.hrp_paste_measurement = fread(fid, 1, 'int32'); % 2216 + 4
   end

   if hdr.file_version > 44 & hdr.file_version < 68
      hdr.graph_type = fread(fid, 1, 'int32');		% 2220 + 4
      hdr.mmt_calc_expr = fread(fid, [40 256], '*char'); % 2224 + 10240
      hdr.mmt_moment_order = fread(fid, 40, 'int32');	% 12464 + 160
      hdr.mmt_time_delay = fread(fid, 40, 'int32');	% 12624 + 160
      hdr.mmt_embed_dim = fread(fid, 40, 'int32');	% 12784 + 160
      hdr.mmt_mi_delay = fread(fid, 40, 'int32');	% 12944 + 160
   end

   if hdr.file_version >= 68
      unused = fread(fid, 411, 'int16')';		% 150 + 822
      hdr.compressed = fread(fid, 1, 'int32');		% 972 + 4
   end

   return;						% graph


%---------------------------------------------------------------------
function hdr = per_chan_data(fid, file_version)

   %  Struct						% off + size
   hdr.chan_header_len = fread(fid, 1, 'int32')';	% 0 + 4
   hdr.num = fread(fid, 1, 'int16')';			% 4 + 2
   hdr.comment_text = deblank(fread(fid, 40, '*char')'); % 6 + 40
   hdr.rgb_color = fread(fid, 1, 'int32')';		% 46 + 4
   hdr.disp_chan = fread(fid, 1, 'int16')';		% 50 + 2
   hdr.volt_offset = fread(fid, 1, 'double')';		% 52 + 8
   hdr.volt_scale = fread(fid, 1, 'double')';		% 60 + 8
   hdr.units_text = deblank(fread(fid, 20, '*char')');	% 68 + 20
   hdr.buf_length = fread(fid, 1, 'int32')';		% 88 + 4
   hdr.ampl_scale = fread(fid, 1, 'double')';		% 92 + 8
   hdr.ampl_offset = fread(fid, 1, 'double')';		% 100 + 8
   hdr.chan_order = fread(fid, 1, 'int16')';		% 108 + 2
   hdr.disp_size = fread(fid, 1, 'int16')';		% 110 + 2


   if file_version >= 68

      unused = fread(fid, 5, 'double')';		% 112 + 40
      hdr.var_sample_divider = fread(fid, 1, 'int16')'; % 152 + 2

      % FIX according to comments on fileexchange
      hdr.var_sample_divider = 1;

   else

   hdr.plot_mode = fread(fid, 1, 'int16')';		% 112 + 2
   hdr.mid = fread(fid, 1, 'double')';			% 114 + 8

   %  interpret rbg_color
   %
   color_str = sprintf('%06s', dec2hex(hdr.rgb_color));
   hdr.rgb_color = ...
      [hex2dec(color_str(5:6)) hex2dec(color_str(3:4)) hex2dec(color_str(1:2))]/255;

   if file_version > 37
      hdr.description = deblank(fread(fid, 128, '*char')'); % 122 + 128
      hdr.var_sample_divider = fread(fid, 1, 'int16')'; % 250 + 2
      
      % AM does not make sense to have a zero divider, set those to 1. No
      % warranties this is standard, but .ACQ files generated with biopac
      % programs produce files with var_sample_divider = 0 although this is
      % not documented. This fix works for the files seen so far.
      if hdr.var_sample_divider==0
          hdr.var_sample_divider = 1; 
      end
   else
      hdr.var_sample_divider = 1;
   end

   if file_version > 38
      hdr.vert_precision = fread(fid, 1, 'int16')';	% 252 + 2
   end

   if file_version > 42
      hdr.active_seg_color = fread(fid, 1, 'int32')';	% 254 + 4
      hdr.active_seg_style = fread(fid, 1, 'int32')';	% 258 + 4

      %  interpret active_seg_color
      %
      color_str = sprintf('%06s', dec2hex(hdr.active_seg_color));
      hdr.active_seg_color = ...
         [hex2dec(color_str(5:6)) hex2dec(color_str(3:4)) hex2dec(color_str(1:2))]/255;
   end

   end						% if file_version < 68

   return;						% per_chan_data


%---------------------------------------------------------------------
function hdr = foreign(fid, file_version)

   %  Struct						% off + size
   hdr.length = fread(fid, 1, 'short')';		% 0 + 2
   hdr.id = fread(fid, 1, 'short')';			% 2 + 2

   if file_version < 68
      hdr.by_foreign_data = fread(fid, hdr.length-4, 'int8')'; % 4 + x
      hdr.length2 = hdr.length;
   elseif file_version < 80
      unused = fread(fid, 1, 'int32')';			% 4 + 4
      hdr.length2 = fread(fid, 1, 'int32')' + 8;	% 8 + 4
   else
      unused = fread(fid, 1, 'int32')';			% 4 + 4
      hdr.length2 = hdr.length + 8;
   end

   return;						% foreign


%---------------------------------------------------------------------
function hdr = per_chan_type(fid)

   %  Struct						% off + size
   hdr.size = fread(fid, 1, 'short')';			% 0 + 2
   hdr.type = fread(fid, 1, 'short')';			% 2 + 2

   return;						% per_chan_type


%---------------------------------------------------------------------
function data = read_acq(fid, hdr, chan_by_chan)

   start_real_chan = hdr.graph.ext_item_header_len;

   for i = 1:hdr.graph.num_channels
      start_real_chan = start_real_chan + hdr.per_chan_data(i).chan_header_len;

      if hdr.per_chan_type(i).type ~= 1			% if integer
         hdr.per_chan_type(i).size = 2;			% only use int16
      end
   end

   start_real_chan = start_real_chan + hdr.foreign.length2 + 4*hdr.graph.num_channels;
   
   sample_divider = [hdr.per_chan_data.var_sample_divider];

   if length(unique(sample_divider))==1 & unique(sample_divider)==1
      min_len = min([hdr.per_chan_data.buf_length]);
   else
      min_len = min([hdr.per_chan_data.buf_length].*sample_divider);
   end

   if length(unique([hdr.per_chan_type.type])) & ~chan_by_chan & ...
	length(unique(sample_divider))==1 & unique(sample_divider)==1

      half_chan = round(hdr.graph.num_channels/2);

      for i = 1:half_chan
         fprintf('.');
      end

      if hdr.per_chan_type(i).type == 1			% double

         data=fread(fid,[hdr.graph.num_channels min_len],'double');
         data=data';
      else						% int

         data=fread(fid,[hdr.graph.num_channels min_len],'int16');
         data=data'.*(ones(min_len,1)*[hdr.per_chan_data.ampl_scale]) ...
		    + ones(min_len,1)*[hdr.per_chan_data.ampl_offset];
      end

      for i = 1:(hdr.graph.num_channels-half_chan)
         fprintf('.');
      end

   else				% if we have to do it chan_by_chan

      if length(unique(sample_divider))==1 & unique(sample_divider)==1

         %  Since data are arranged like: "channel in sample"
         %  { s1 of {ch1 ch2 ...}, s2 of {ch1 ch2 ...} ... }
         %  We can either read data sample by sample (all chan at once),
         %  or, channel by channel, but need to skip the rest of chan
         %
         size_all_chan_per_sample = 0;

         for i = 1:hdr.graph.num_channels
            size_all_chan_per_sample = size_all_chan_per_sample + ...
					hdr.per_chan_type(i).size;
         end

         for i = 1:hdr.graph.num_channels

            fprintf('.');
            fseek(fid, start_real_chan, 'bof');

            %  First jump to the start point of the right channel
            %
            start_chan = 0;

            for j = 1 : (i-1)
               start_chan = start_chan + hdr.per_chan_type(j).size;
            end

            fseek(fid, start_chan, 'cof');

            %  Need to skip the rest of chan, in order to read each sample point
            %
            skip_chan = size_all_chan_per_sample - hdr.per_chan_type(i).size;

            if hdr.per_chan_type(i).type == 1		% double

               tmp=fread(fid,hdr.per_chan_data(i).buf_length,'double',skip_chan);

            else					% int

               %  Need to be scaled & shifted by ampl_scale & ampl_offset
               %  for integer data
               %
               tmp=fread(fid,hdr.per_chan_data(i).buf_length,'int16',skip_chan) ...
			* hdr.per_chan_data(i).ampl_scale ...
			+ hdr.per_chan_data(i).ampl_offset;
            end

            data(:,i) = tmp(1:min_len);
         end

      else				% sample_divider>1 or different

         %  If there is a sample_divider>1 in any channel, we have to
         %  read data sample by sample (very slow!)
         %
         data = zeros(min_len, hdr.graph.num_channels);
         mask=zeros(min_len, hdr.graph.num_channels);
         sample_divider = abs(sample_divider);

         for j = 1:hdr.graph.num_channels
            mask(1:sample_divider(j):min_len, j)=1;
         end

         for i = 1:min_len

            if mod(i-1, ceil(min_len/hdr.graph.num_channels))==0
               fprintf('.');
            end

            for j = 1:hdr.graph.num_channels
               if mask(i,j) | sample_divider(j)==1
                  if hdr.per_chan_type(j).type == 1		% double
                     data(i,j) = fread(fid,1,'double');
                  else
                     data(i,j) = fread(fid,1,'int16');		% int

                     %  Need to be scaled & shifted by ampl_scale & ampl_offset
                     %  for integer data
                     %
                     data(i,j) = data(i,j) * hdr.per_chan_data(j).ampl_scale ...
					+ hdr.per_chan_data(j).ampl_offset;
                  end
               else
                  data(i,j) = data(i-1,j);
               end		% if mask ... read sample_divider

            end		% for j
         end		% for i

      end	% if length(unique(sample_divider
   end	% length(unique([hdr.per_chan_type.type

   return;					% read_acq


 %---------------------------------------------------------------------
 function info = read_markers(fid, hdr, data_size)
     
     % AM foolproof way of getting to the markers: see where data starts and
     % add the size of the data..
     
     start_real_chan = hdr.graph.ext_item_header_len;

     for i = 1:hdr.graph.num_channels
         start_real_chan = start_real_chan + hdr.per_chan_data(i).chan_header_len;

         if hdr.per_chan_type(i).type ~= 1			% if integer
             hdr.per_chan_type(i).size = 2;			% only use int16
         end
     end

     start_real_chan = start_real_chan + hdr.foreign.length2 + 4*hdr.graph.num_channels;

     size_channel_data = sum([hdr.per_chan_type.size])*data_size(1);
     start_markers = start_real_chan + size_channel_data;

     % AM place the file position in the position where data ends
     fseek(fid,start_markers,'bof');

     % AM then follow specifications; borrowed from code of the
     % non-functioning (at least for me) ACQREAD,
     %    ACQREAD, version 2.0 (2007-08-21)
     %    Copyright (C) 2006-2007  Sebastien Authier and Vincent Finnerty

     info.lLength = fread(fid,1,'*int32');
     info.lMarkers = fread(fid,1,'*int32');		% Number of markers
     if (info.lLength > 0) & (info.lMarkers > 0)
         for n = 1:double(info.lMarkers)
             info.lSample(n) = fread(fid,1,'*int32');	% Location of marker
             info.fSelected(n) = fread(fid,1,'*int16');
             info.fTextLocked(n) = fread(fid,1,'*int16');
             info.fPositionLocked(n) = fread(fid,1,'*int16');
             info.nTextLength(n) = fread(fid,1,'*int16');  % Length of marker text string
             info.szText{n} = deblank(fread(fid,double(info.nTextLength(n))+1,'*char')');  % Marker text string
         end
     else
         info.lSample = [];
         info.fSelected = [];
         info.fTextLocked = [];
         info.fPositionLocked = [];
         info.nTextLength = [];
         info.szText = [];
     end

     return;

