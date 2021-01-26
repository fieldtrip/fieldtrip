% pop_editeventvals() - Edit events contained in an EEG dataset structure. 
%               If the dataset is the only input, a window pops up 
%               allowing the user to insert the relevant parameter values.
%
% Usage: >> EEGOUT = pop_editeventvals( EEG, 'key1', value1, ...
%                                    'key2', value2, ... );
% Input:
%   EEG  - EEG dataset
%
% Optional inputs:
%   'sort'        - { field1 dir1 field2 dir2 } Sort events based on field1
%                   then on optional field2. Arg dir1 indicates the sort 
%                   direction (0 = increasing, 1 = decreasing).
%   'changefield' - {num field value} Insert the given value into the specified 
%                   field in event num. (Ex: {34 'latency' 320.4})
%   'changeevent' - {num value1 value2 value3 ...} Change the values of
%                   all fields in event num.
%   'add','append','insert' - {num value1 value2 value3 ...} Insert event
%                   before or at event num, and assign value to structure 
%                   fields. Note that the latency field must be in second 
%                   and will be converted to data sample. Note also that 
%                   the index of the event is often irrelevant, as events
%                   will be automatically resorted by latencies.
%   'delete'      - vector of indices of events to delete
%
% Outputs:
%   EEGOUT        - EEG dataset with the selected events only
%
% Ex:  EEG = pop_editeventvals(EEG,'changefield', { 1 'type' 'target'});
%        % set field type of event number 1 to 'target'
%
% Author: Arnaud Delorme & Hilit Serby, SCCN, UCSD, 15 March 2002
%
% See also: pop_selectevent(), pop_importevent()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 15 March 2002, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 03-16-02 text interface editing -sm & ad 
% 03-18-02 automatic latency switching display (epoch/continuous) - ad & sm 
% 03-18-02 debug soring order - ad
% 03-18-02 put latencies in ms - ad, lf & sm
% 03-29-02 debug latencies in ms - ad & sm
% 04-02-02 debuging test - ad & sm

function [EEG, com] = pop_editeventvals(EEG, varargin);

com ='';
if nargin < 1
   help pop_editeventvals;
   return;
end;	

if nargin >= 2 | isstr(EEG) % interpreting command from GUI or command line
    
    if isstr(EEG) % GUI
        gui    = 1;
        varargin = { EEG varargin{:} };
        
        % user data
        % ---------
        userdata  = get(gcf, 'userdata');
        EEG       = userdata{1};
        oldcom    = userdata{2};
        allfields = fieldnames(EEG.event);
        tmpind    = strmatch('urevent', allfields);
        allfields(tmpind) = [];
        
        % current event
        % -------------
        objevent  = findobj('parent', gcf, 'tag', 'numval');
        valnum    = str2num(get(objevent, 'string'));
        shift     = 0;
    
    else % command line    
        
        gui = 0;
        if isempty(EEG.event)
            disp('Getevent: cannot deal with empty event structure');
            return;
        end;   
        
        allfields = fieldnames(EEG.event);
        tmpind = strmatch('urevent', allfields);
        allfields(tmpind) = [];
        
    end;
    
    % retinterpret inputs (fix BUG 454)
    % --------------------------------
    if ~gui
        newvararg = {};
        for indfield = 1:2:length(varargin)
            
            com     = varargin{ indfield   };
            tmpargs = varargin{ indfield+1 };
            newvararg = { newvararg{:} com };
            
            if any(strcmpi({'add','insert','append'},com))
                evtind     = tmpargs{1};
                fields     = fieldnames(EEG.event);
                emptycells = cell(1,length(fields)-1);
                newvararg  = { newvararg{:}, { evtind emptycells{:} } };
                if strcmpi(com, 'append'), evtind = evtind+1; end;
                
                for ind = 2:length( tmpargs )
                    if ~strcmpi(fields{ind-1}, 'urevent')
                        newvararg = { newvararg{:},'changefield',{ evtind fields{ind-1} tmpargs{ind} } };
                    end;
                end;
            else
                newvararg = { newvararg{:} tmpargs };
            end;
        end;
        varargin = newvararg;
    end;
    
    % scan inputs
    % -----------
    for indfield = 1:2:length(varargin)
        if length(varargin) >= indfield+1
            tmparg = varargin{ indfield+1 };
        end;
        
        switch lower(varargin{indfield})
        
     case 'goto', % ******************** GUI ONLY ***********************
      
      % shift time
      % ----------
      shift     = tmparg;
      valnum    = valnum + shift;
      if valnum < 1,                valnum = 1;                end;
      if valnum > length(EEG.event), valnum = length(EEG.event); end;
      set(objevent, 'string', num2str(valnum,5));

      % update fields
      % -------------
      for index = 1:length(allfields) 
          
          enable = 'on';
          if isfield(EEG.event, 'type') 
              if strcmpi(EEG.event(valnum).type, 'boundary'), enable = 'off'; end;
          end;
          if strcmp( allfields{index}, 'latency') & ~isempty(EEG.event(valnum).latency)
              if isfield(EEG.event, 'epoch')
                   value = eeg_point2lat( EEG.event(valnum).latency, EEG.event(valnum).epoch, ...
                                          EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
              else value = (EEG.event(valnum).latency-1)/EEG.srate+EEG.xmin;
              end;
          elseif strcmp( allfields{index}, 'duration') & ~isempty(EEG.event(valnum).duration)
              if isfield(EEG.event, 'epoch')
                   value = EEG.event(valnum).duration/EEG.srate*1000; % milliseconds
              else value = EEG.event(valnum).duration/EEG.srate;      % seconds
              end;
          else
              value = getfield( EEG.event(valnum), allfields{index});
          end;
          
          % update interface
          % ----------------
          tmpobj = findobj('parent', gcf, 'tag', allfields{index});
          set(tmpobj, 'string', num2str(value,5), 'enable', enable);
      end;
      
      % update original
      % --------------- 
      tmpobj = findobj('parent', gcf, 'tag', 'original');
      if isfield(EEG.event, 'urevent') & EEG.event(valnum).urevent ~= valnum
           set(tmpobj, 'string', [ 'originally ' int2str(EEG.event(valnum).urevent)], ...
                       'horizontalalignment', 'center');
      else set(tmpobj, 'string', ' '); 
      end;
      return; % NO NEED TO SAVE ANYTHING
          
     case { 'append' 'insert' 'add' }, % **********************************************************

      if gui
          shift     = tmparg; % shift is for adding before or after the event

          % add epoch number if data epoch
          % ------------------------------
          tmpcell    = cell(1,1+length(fieldnames(EEG.event))); 
          tmpcell{1} = valnum;
          if EEG.trials > 1
              indepoch = strmatch('epoch', fieldnames(EEG.event), 'exact');
              if valnum > 1, tmpprevval = valnum-1;
              else           tmpprevval = valnum+1;
              end;
              if tmpprevval <= length(EEG.event)
                  tmpcell{indepoch+1} = EEG.event(tmpprevval).epoch;
              end;
          end;
          
          % update commands
          % ---------------
          if shift
              oldcom     = { oldcom{:} 'append', tmpcell };
          else
              oldcom     = { oldcom{:} 'insert', tmpcell };
          end;
      else
          if strcmpi(lower(varargin{indfield}), 'append') % not 'add' for backward compatibility
               shift = 1;
          else shift = 0;
          end;
          valnum = tmparg{1};
      end;
      
      % find ur index
      % -------------
      valnum    = valnum   + shift;
      if isfield(EEG.event, 'epoch'), curepoch = EEG.event(valnum).epoch; end;
      if isfield(EEG, 'urevent') & isfield(EEG.event, 'urevent')
          % find non empty urvalnum
          urvalnum = [];
          count = 0;
          while isempty(urvalnum)
              tmpindex = mod(valnum+count-1, length(EEG.event)+1)+1;
              urvalnum = EEG.event(valnum+count).urevent;
              count = count+1;
          end;
          if isfield(EEG.urevent, 'epoch'), urcurepoch = EEG.urevent(urvalnum).epoch; end;
          urvalnum = urvalnum;
      end;
      
      % update urevents
      % ---------------
      if isfield(EEG, 'urevent') & isfield(EEG.event, 'urevent')
          EEG.urevent(end+3)            = EEG.urevent(end);
          EEG.urevent(urvalnum+1:end-2) = EEG.urevent(urvalnum:end-3);
          EEG.urevent(urvalnum)         = EEG.urevent(end-1);
          EEG.urevent                   = EEG.urevent(1:end-2);
          if isfield(EEG.urevent, 'epoch'),  EEG.urevent(urvalnum).epoch = urcurepoch; end;
      end;
      
      % update events
      % -------------
      EEG.event(end+3)          = EEG.event(end);
      EEG.event(valnum+1:end-2) = EEG.event(valnum:end-3);
      EEG.event(valnum)         = EEG.event(end-1);
      EEG.event                 = EEG.event(1:end-2);
      if isfield(EEG.event, 'epoch'), EEG.event(valnum).epoch = curepoch; end;      
      if isfield(EEG, 'urevent') & isfield(EEG.event, 'urevent')
          EEG.event(valnum).urevent = urvalnum;
          for index = valnum+1:length(EEG.event)
              EEG.event(index).urevent = EEG.event(index).urevent+1;
          end;
      end;
      
      % update type field
      % -----------------
      for tmpind = 1:length(allfields)
          EEG.event = checkconsistency(EEG.event, valnum, allfields{tmpind});
      end;
      
     case 'delete', % **********************************************************
      
      if ~gui
          valnum = tmparg;
      end
                
      EEG.event(valnum) = []; 
      
      if gui, 
          
          valnum           = min(valnum,length(EEG.event));
          set(objevent, 'string', num2str(valnum)); 
          
          % update commands
          % ---------------
          oldcom = { oldcom{:} 'delete', valnum };
      
      end;
    
     case { 'assign' 'changefield' }, % **********************************************************
      
      if gui, % GUI case
          
          field    = tmparg;
          objfield = findobj('parent', gcf, 'tag', field);
          editval     = get(objfield, 'string');
          if ~isempty(editval) & ~isempty(str2num(editval)), editval = str2num(editval); end;

          % update history
          % --------------
          oldcom = { oldcom{:},'changefield',{ valnum field editval }};
      
      else % command line case
          
          valnum  = tmparg{1};
          field   = tmparg{2};
          editval = tmparg{3};
          
      end;
          
      % latency and duration case
      % -------------------------
      if strcmp( field, 'latency') & ~isempty(editval)
          if isfield(EEG.event, 'epoch')
               editval = eeg_lat2point( editval, EEG.event(valnum).epoch, ...
                                       EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
          else editval = (editval- EEG.xmin)*EEG.srate+1;
          end;
      end;
      if strcmp( field, 'duration') & ~isempty(editval)
          if isfield(EEG.event, 'epoch')
               editval = editval/1000*EEG.srate; % milliseconds
          else editval = editval*EEG.srate;      % seconds
          end;
      end;
      
      % adapt to other formats
      % ----------------------
      EEG.event(valnum) = setfield(EEG.event(valnum), field, editval);
      EEG.event = checkconsistency(EEG.event, valnum, field);
      
      % update urevents
      % ---------------
      if isfield(EEG, 'urevent') & isfield(EEG.event, 'urevent')
          urvalnum  = EEG.event(valnum).urevent;
          
          % latency case
          % ------------
          if strcmp( field, 'latency') & ~isempty(editval)
              if isfield(EEG.urevent, 'epoch')
                  urepoch = EEG.urevent(urvalnum).epoch;
                  
                  
                  % find closest event latency
                  % --------------------------
                  if valnum<length(EEG.event)
                      if EEG.event(valnum+1).epoch == urepoch
                          urlatency = EEG.urevent(EEG.event(valnum+1).urevent).latency;
                          latency   = EEG.event(valnum+1).latency;
                      end;
                  end;
                  if valnum>1
                      if EEG.event(valnum-1).epoch == urepoch
                          urlatency = EEG.urevent(EEG.event(valnum-1).urevent).latency;
                          latency   = EEG.event(valnum-1).latency;
                      end;
                  end;
                  
                  % update event
                  % ------------
                  if exist('urlatency') ~=1
                      disp('Urevent not updated: could not find other event in the epoch');
                  else
                      editval = urlatency - ( latency - editval ); % new latency value
                  end;
              else 
                  editval = eeg_urlatency(EEG.event, EEG.event(valnum).latency);
              end;
          elseif strcmp( field, 'latency') % empty editval
              EEG.event(valnum).latency = NaN;
          end;
          
          % duration case
          % ------------
          if strcmp( field, 'duration') & ~isempty(editval)
              if isfield(EEG.event, 'epoch')
                   editval = editval/1000*EEG.srate; % milliseconds -> point
              else editval = editval*EEG.srate;      % seconds -> point
              end;
          end;
          EEG.urevent = setfield(EEG.urevent, {urvalnum}, field, editval);
      end;
      
     case 'sort', % **********************************************************
      
      if gui % retrieve data
          
          field1 = get(findobj('parent', gcf, 'tag', 'listbox1'), 'value');
          field2 = get(findobj('parent', gcf, 'tag', 'listbox2'), 'value');
          dir1   = get(findobj('parent', gcf, 'tag', 'order1'),   'value');
          dir2   = get(findobj('parent', gcf, 'tag', 'order2'),   'value');
          
          if field1 > 1, field1 = allfields{field1-1}; else return; end;
          if field2 > 1, field1 = allfields{field2-1}; else field2 = []; end;
          
          % update history
          % --------------
          oldcom = { oldcom{:},'sort',{ field1 dir1 field2 dir2 } };

      else % command line
          field1 = tmparg{1};
          if length(tmparg) < 2, dir1 = 0;
          else                   dir1 = tmparg{2}; 
          end;
          if length(tmparg) < 3, field2 = [];
          else                   field2 = tmparg{3}; 
          end;
          if length(tmparg) < 4, dir2 = 0;
          else                   dir2 = tmparg{4}; 
          end;
      end;
      
      % Toby edit 11/16/2005 This section is scrambling the eeg.event
      % fields. Requires further investigation.
 
      if ~isempty(field2)
          tmpevent = EEG.event;
          if ~ischar(EEG.event(1).(field2))
              tmparray = [ tmpevent.(field2) ];              
          else tmparray = { tmpevent.(field2) };
          end
          % Commented out 11/18/2005, Toby
          % These lines were incorrectly sorting the event.latency field in
          % units of time (seconds) relevant to each event's relative
          % latency time as measured from the start of each epoch. It is
          % possible that there are occasions when it is desirable to do
          % so, but in the case of pop_mergset() it is not. 
          
          %if strcmp(field2, 'latency') & EEG.trials > 1
          %    tmparray = eeg_point2lat(tmparray, {EEG.event.epoch}, EEG.srate, [EEG.xmin EEG.xmax], 1);
          %end;
          [X I] = mysort( tmparray );
          if dir2 == 1, I = I(end:-1:1); end;
          events = EEG.event(I);
      else
          events = EEG.event;
      end;  
      tmpevent = EEG.event;
      if ~ischar(EEG.event(1).(field1))
          tmparray = [ tmpevent.(field1) ];
      else tmparray = { tmpevent.(field1) };
      end
      % Commented out 11/18/2005, Toby
      %if strcmp( field1, 'latency') & EEG.trials > 1
      %    tmparray = eeg_point2lat(tmparray, {events.epoch}, EEG.srate, [EEG.xmin EEG.xmax], 1);
      %end;
      [X I] = mysort( tmparray );
      if dir1 == 1, I = I(end:-1:1); end;
      EEG.event = events(I);

      
      if gui
          % warn user
          % ---------
          warndlg2('Sorting done');
      else 
          noeventcheck  = 1; % otherwise infinite recursion with eeg_checkset
      end;
      
    end; % end switch
    end; % end loop

    % save userdata
    % -------------
    if gui
        userdata{1} = EEG;
        userdata{2} = oldcom;
        set(gcf, 'userdata', userdata);
        pop_editeventvals('goto', shift);              
    else
        if ~exist('noeventcheck','var')
            EEG = eeg_checkset(EEG, 'eventconsistency');
            EEG = eeg_checkset(EEG, 'checkur');
        end;
    end;
    return;
end;

% ----------------------
% graphic interface part
% ----------------------

if isempty(EEG.event)
    disp('Getevent: cannot deal with empty event structure');
    return;
end;   

allfields = fieldnames(EEG.event);
tmpind = strmatch('urevent', allfields);
allfields(tmpind) = [];

if nargin<2
    % add field values
    % ----------------
    geometry = { [2 0.5] };
    tmpstr = sprintf('Edit event field values (currently %d events)',length(EEG.event));
    uilist = { { 'Style', 'text', 'string', tmpstr, 'fontweight', 'bold'  } ...
               { 'Style', 'pushbutton', 'string', 'Delete event', 'callback', 'pop_editeventvals(''delete'');'  }};

    for index = 1:length(allfields) 

        geometry = { geometry{:} [1 1 1 1] };
        
        % input string
        % ------------
        if strcmp( allfields{index}, 'latency') | strcmp( allfields{index}, 'duration') 
            if EEG.trials > 1
                 inputstr =  [ allfields{index} ' (ms)'];
            else inputstr =  [ allfields{index} ' (sec)'];
            end;   
		else inputstr =  allfields{index};
		end;
        
		% callback for displaying help
		% ----------------------------
        if index <= length( EEG.eventdescription )
             tmptext = EEG.eventdescription{ index };
			 if ~isempty(tmptext)
				 if size(tmptext,1) > 15,    stringtext = [ tmptext(1,1:15) '...' ]; 
				 else                        stringtext = tmptext(1,:); 
				 end;
			 else stringtext = 'no-description'; tmptext = 'no-description';
			 end;
        else stringtext = 'no-description'; tmptext = 'no-description';
        end;
		cbbutton = ['questdlg2(' vararg2str(tmptext) ...
					',''Description of field ''''' allfields{index} ''''''', ''OK'', ''OK'');' ];

        % create control
        % --------------
        cbedit = [ 'pop_editeventvals(''assign'', ''' allfields{index} ''');' ]; 
		uilist   = { uilist{:}, { }, ...
					 { 'Style', 'pushbutton', 'string', inputstr, 'callback',cbbutton  }, ...
					 { 'Style', 'edit', 'tag', allfields{index}, 'string', '', 'callback', cbedit } ...
                     { } };
    end;

    % add buttons
    % -----------
    geometry = { geometry{:} [1] [1.2 0.6 0.6 1 0.6 0.6 1.2] [1.2 0.6 0.6 1 0.6 0.6 1.2] [2 1 2] };
    
    tpappend = 'Append event after the current event';
    tpinsert = 'Insert event before the current event';
    tporigin = 'Original index of the event (in EEG.urevent table)';
    uilist   = { uilist{:}, ...
          { }, ...
          { },{ },{ }, {'Style', 'text', 'string', 'Event Num', 'fontweight', 'bold' }, { },{ },{ }, ...
          { 'Style', 'pushbutton', 'string', 'Insert event',  'callback', 'pop_editeventvals(''append'', 0);', 'tooltipstring', tpinsert }, ...
          { 'Style', 'pushbutton', 'string', '<<',            'callback', 'pop_editeventvals(''goto'', -10);' }, ...
          { 'Style', 'pushbutton', 'string', '<',             'callback', 'pop_editeventvals(''goto'', -1);' }, ...
          { 'Style', 'edit',       'string', '1',             'callback', 'pop_editeventvals(''goto'', 0);', 'tag', 'numval' }, ...
          { 'Style', 'pushbutton', 'string', '>',             'callback', 'pop_editeventvals(''goto'', 1);' }, ...
          { 'Style', 'pushbutton', 'string', '>>',            'callback', 'pop_editeventvals(''goto'', 10);' }, ...
          { 'Style', 'pushbutton', 'string', 'Append event',  'callback', 'pop_editeventvals(''append'', 1);', 'tooltipstring', tpappend }, ...
          { }, { 'Style', 'text',  'string', ' ', 'tag', 'original' 'horizontalalignment' 'center' 'tooltipstring' tporigin } { } };

    % add sorting options
    % -------------------
    listboxtext = 'No field selected';  
    for index = 1:length(allfields) 
         listboxtext = [ listboxtext '|' allfields{index} ]; 
    end;
    geometry = { geometry{:} [1] [1 1 1] [1 1 1] [1 1 1] };
    uilist = {  uilist{:}, ...
         { 'Style', 'text',       'string', 'Re-order events (for review only)', 'fontweight', 'bold'  }, ...
         { 'Style', 'text',       'string', 'Main sorting field:'  }, ...
         { 'Style', 'popupmenu',  'string', listboxtext, 'tag', 'listbox1' }, ...
         { 'Style', 'checkbox',   'string', 'Click for decreasing order', 'tag', 'order1' } ...
         { 'Style', 'text',       'string', 'Secondary sorting field:'  }, ...
         { 'Style', 'popupmenu',  'string', listboxtext, 'tag', 'listbox2' }, ...
         { 'Style', 'checkbox',   'string', 'Click for decreasing order', 'tag', 'order2' }, ...
         { } { 'Style', 'pushbutton', 'string', 'Re-sort', 'callback', 'pop_editeventvals(''sort'');' }, ...
         { } };
   
    userdata = { EEG {} };
    inputgui( geometry, uilist, 'pophelp(''pop_editeventvals'');', ...
                                  'Edit event values -- pop_editeventvals()', userdata, 'plot');
    pop_editeventvals('goto', 0);
    
    % wait for figure
    % ---------------
    fig = gcf;
    waitfor( findobj('parent', fig, 'tag', 'ok'), 'userdata');
    try, userdata = get(fig, 'userdata'); close(fig); % figure still exist ?
    catch, return; end;
    
    % transfer events
    % ---------------
    if ~isempty(userdata{2})
        com = sprintf('%s = pop_editeventvals(%s,%s);', inputname(1), inputname(1), vararg2str(userdata{2}));
    end;
    if isempty(findstr('''sort''', com))
        if ~isempty(userdata{2}) % some modification have been done
            EEG       = userdata{1};
            disp('Checking event consistency...');
            EEG = eeg_checkset(EEG, 'eventconsistency');
            EEG = eeg_checkset(EEG, 'checkur');
        end;
    else 
        com = '';
        disp('WARNING: all edits discarded because of event resorting. The EEGLAB event structure');
        disp('            must contain events sorted by latency (you may obtain an event structure');
        disp('            with resorted event by calling this function from the command line).');
    end;
    return;
    
end;
return;

% format the output field
% -----------------------
function strval = reformat( val, latencycondition, trialcondition, eventnum)
    if latencycondition
        if trialcondition
            strval = ['eeg_lat2point(' num2str(val) ', EEG.event(' int2str(eventnum) ').epoch, EEG.srate,[EEG.xmin EEG.xmax]*1000, 1E-3);' ];
        else    
            strval = [ '(' num2str(val) '-EEG.xmin)*EEG.srate+1;' ]; 
        end;
    else
        if isstr(val), strval = [ '''' val '''' ];
        else           strval = num2str(val);
        end;
    end;

% sort also empty values
% ----------------------
function [X, I] = mysort(tmparray);
    if iscell(tmparray)
        if all(cellfun('isreal', tmparray))
            tmpempty = cellfun('isempty', tmparray);
            tmparray(tmpempty) = { 0 };
            tmparray = [ tmparray{:} ];
        end;
    end;
    try, 
        [X I] = sort(tmparray);
    catch,
        disp('Sorting failed. Check that selected fields contain uniform value format.');
        X = tmparray;
        I = 1:length(X);
    end;
    
% checkconsistency of new event
% -----------------------------
function eventtmp = checkconsistency(eventtmp, valnum, field)
    
    otherval = mod(valnum+1, length(eventtmp))+1;
    
    if isstr(getfield(eventtmp(valnum), field)) & ~isstr(getfield(eventtmp(otherval), field))
        eventtmp(valnum) = setfield(eventtmp(valnum), field, str2num(getfield(eventtmp(valnum), field)));
    end;
    if ~isstr(getfield(eventtmp(valnum), field)) & isstr(getfield(eventtmp(otherval), field))
        eventtmp(valnum) = setfield(eventtmp(valnum), field, num2str(getfield(eventtmp(valnum), field)));
    end;
    if strcmpi(field, 'latency') & isempty(getfield(eventtmp(valnum), field))
        eventtmp(valnum).latency = NaN;
    end;
