function [data,ds]=setCTFDataBalance(data,ds,balance1,unit,chanList,messages);

%  Version 1.1  13 April 2007  Mod to chanList:  If chanList is omitted and
%                              size(data,2)<ds.res4.no_channels, then setCTFDataBalance
%                              sets chanList=1:size(data,2) and prints a warning the first
%                              time it happens.
%               7 Dec. 2006.   Fixed a bug.  Changed calls to getCTFBalanceCoefs so 
%               setCTFDataBalance would balance and unbalance reference gradiometers.
%  Version 1.0   27 October 2006

%  Adjusts the gradient order balance from balance0 to balance1.
%  Inputs : 
%    data : Array of data to be balanced.  precision= single or double.
%           size(data)=[Npt Nchan Ntrial].  If chanList is not supplied, it is assumed that
%           array data has dataset channels 1:Nchan (i.e. data is not just SQUID sensors
%           although only SQUID channels are affected by setCTFDataBalance.)
%    ds : Structure produced by readCTFds that describes the dataset.
%           balance1 : Character array defining balance of the output data.
%                   Parameter balance1 must be specified, even if it is only balance1=[];
%                   balance1=[],' ' or 'NONE' : Output is unbalanced.
%                   size(balance,1)=1 : MEG balancing. Assume no reference balancing.
%                   size(balance,1)=2 : MEG balance = balance(1,:).
%                                       Gref balance = balance(2,:)
%    unit : Options : [],' ', 'fT', or 'T' : Physical unit
%                     'int' or 'phi0' : Raw units.
%                     If the unit argument is not included, unit is set to 'fT'.

%    chanList : Optional list of dataset channels included in array data.
%               Channel numbering is referred to structure ds (the output of readCTFds).
%               If length(chanList)=Nchan.
%                        data(:,k,q)=trial q of channel chanList(k) of the dataset.
%               If chanList=[],<0 or not defined : Array data has dataset channels 1:Nchan.
%               Otherwise, length(chanList)=Nchan=size(data,2).
%               If length(chanList)~=size(data,2) : Error.

%    messages=0: Don't print tracing message
%                  =1: Print a message when there is a change in requested balance.
%                  =2: Always print a message

%  Outputs: data: Data with balanced SQUID signals.  Set data=[] on failure.
%             ds: ds structure for balanced data.    Set ds=-1 on failure.

%  In any event, print error messages.if exist('messages')~=1;messages=1;end
%  The balance and unbalance is made complicated because this code allows for situations 
%  where the reference gradiometers are also balanced (although this option seems never
%  to be used).

% Function calls:  Most functions are called only once, and are intended to be 
%       called only by setCTFDataBalance or its subprograms.

%   - getCTFBalanceCoefs.  Called by balance_data and unbalance_data.  Gets the balance
%                          tables from structure ds.

%    Included in this listing:
%   - getDsBalance: Gets the current balance state of the data and checks for error in
%                     the grad_order_no field of ds.rs4.senres.
%   - check_balance : Makes balance1 into [2 4] character arrays.
%   - unbalance_data : Converts data from balance0 to unbalanced.
%   - balance_data: Converts data from unbalanced to balance1.
%   - reindex : Called by functions unbalance_data and balance_data to re-order indices 
%               when chanList does not include all the channels, or the channels have 
%               been reordered in array data.

% If virtual memory is a problem, consider bringing functions balance_data and
% unbalance_data into the code.  In the present configuration, a copy of array data is
% made when these functions are called.

%  In the event of an error, returns data=[], ds1=-1  This should force an error in the 
%  calling program.
persistent chanListWarning

if nargout==0 & nargin==0
  fprintf(['\nsetCTFDataBalance: Version 1.1   13 April 2007   ',...
      'Balances MEG data.\n',...
      '\t\tCall : ',...
      '[data1,ds1]=setCTFDataBalance(data0,ds0,balance1,unit,chanList,messages);\n\n']);
  return
end

buffsize=6e6;   % process blocks of buffsize words (=8*buffsize bytes) at a time
balanceOptions=strvcat('NONE','G1BR','G2BR','G3BR');
unit_options=strvcat('fT','T','phi0','int');
default_unit='fT';

% Check that the inputs are sensible.
if nargin<3
  fprintf(['\nsetCTFDataBalance: Only %d input arguments?  ',...
    'Must specify at least data, ds and balance.\n\n'],nargin);
  data=[];ds=-1;
  return
elseif ~any(strmatch(class(data),strvcat('single','double'))) | isempty(data) | ...
    ~isstruct(ds) | isempty(ds) | isempty(balance1) |~ischar(balance1)
  fprintf('\nsetCTFDataBalance: Wrong argument types or sizes.\n\n');
  whos data ds balance1
  data=[];ds=-1;
  return
elseif ~isfield(ds,'res4')
  fprintf('\nsetCTFDataBalance: Field res4 is missing from structure ds.\n\n');
  ds
  data=[];ds=-1;
  return
end

if nargout~=2
  fprintf(['\nsetCTFDataBalance: %d output arguments specified.  ',...
      'Must have 2 outputs (data,ds).\n\n'],nargout);
  data=[];ds=-1;
  return
end

%  Single or double precision data.
prec=class(data);
% Force upper-case characters
balance1=upper(balance1);

%  Check that arrays data and chanList are consistent and sensible
if size(data,2)>ds.res4.no_channels
  fprintf('\nsetCTFDataBalance: ds.res4.no_channels=%d, but size data=[',...
    ds.res4.no_channels);
  fprintf(' %d',size(data));fprintf(']\n\n');
  data=[];ds=-1;
  return
end
if exist('chanList')~=1;chanList=[];end
if ~isempty(chanList)  % chanList~=[].  Check chanList
  if ~isnumeric(chanList)
    fprintf('\nsetCTFDataBalance: chanList has the wrong type.\n\n');
    whos data ds balance1 chanList
    data=[];ds=-1;
    return
  else  % chanList is numeric and non-empty
    %  Check that chanList is a vector
    if sum(size(chanList)>1)>1
      fprintf('\nsetCTFDataBalance: size(chanList)=[');fprintf(' %d',size(chanList));
      fprintf(']   Not a vector.\n\n');
      data=[];ds=-1;  % Force an error in the calling program
      return
    else
      chanList=reshape(chanList,1,length(chanList)); % chanList is a row vector.
      if any(chanList<=0) | any(chanList>ds.res4.no_channels) | ...
          length(chanList)~=size(data,2)
        fprintf('setCTFDataBalance: Error in input chanList:\n');
        fprintf('   min(chanList)=%d  max(chanList)=%d  Dataset has %d channels.\n',...
          min(chanList),max(chanList),ds.res4.no_channels);
        fprintf('length(chanList)=%d   size(data)=[',length(chanList));
        fprintf(' %d',size(data));
        fprintf(']    Must have length(chanList)=size(data,2)\n');
        data=[];ds=-1;  % Force an error in the calling program
        return
      end
    end
  end
  
else  % chanList=[].  Array data must include all of the channels.
  if size(data,2)<ds.res4.no_channels
    if isempty(chanListWarning)
      fprintf('setCTFDataBalance: No chanList supplied and size(data)=[');
      fprintf(' %d',size(data));fprintf(']\n');
      fprintf('                   Set chanList=1:%d\n',size(data,2));
      fprintf('                   No. of channels=ds.res4.no_channels=%d\n',...
        ds.res4.no_channels);
      chanListWarning=1;
    end
    chanList=1:size(data,2);
  end
end

%  Check the data units, convert to lower case and flag incorrect unit specification
if ~exist('unit');
  unit=default_unit;
elseif isempty(unit);
  unit=default_unit;
elseif ~ischar(unit)
  fprintf('\nsetCTFDataBalance: Input unit has the wrong type: class(unit)=%s\n\n',...
    class(unit));
  data=[];ds=-1;
  return
end
unit=lower(unit);
if isempty(strmatch(unit,lower(unit_options)))
  fprintf('\nsetCTFDataBalance: unit=%s.   Must be one of',unit);
  for k=1:size(unit_options,1);fprintf(' %s',unit_options(k,:));end;fprintf('\n');
  data=[];ds=-1;  % Force an error in the calling program
  return
end

%  Make sure message request is sensible.  This is important only as a way of checking
%  that the input list is sensible.  The value of messages does not affect the output of
%  setCTFDataBalance.
if exist('messages')~=1;messages=1;end
if ~isnumeric(messages)
  fprintf(['\nsetCTFDataBalance: Input messages has the wrong type: ',...
      'class(messages)=%s\n\n'],class(messages));
  data=[];ds=-1;
  return
elseif isempty(messages);
  messages=1;
elseif ~isequal(size(messages),[1 1])
  fprintf('\nsetCTFDataBalance: Input messages doesn''t make sense.  size(messages)=[');
  fprintf(' %d',size(messages));fprintf(']   (Must be [1 1])\n\n');
  data=[];ds=-1;
  return
elseif messages<0 | messages>2
  fprintf('\nsetCTFDataBalance: Input messages=%d.  Must be [], 0, 1 or 2.\n\n',messages);
  data=[];ds=-1;
  return
end

% Process data in blocks of nChBlock channels to reduce virtual memory demand.
nChBlock=max(1,floor(buffsize/size(data,1)));

% badGref = list of bad reference gradiometers referred to dataset channel numbering.
% badGMEG = list of bad MEG gradiometers referred to dataset channel numbering.
%  These gradiometers are considered bad because their grad_order_no is wrong, not
%  because thay are listed in ds.BadChannels 
[balance0,badGref,badGMEG,Greflist,GMEGlist]=getDsBalance(ds);
if strmatch('BAD',balance0)
  fprintf(['\nsetCTFDataBalance: The balance order of the Grefs and/or the MEG ',...
      'channels cannot be determined.\n',...
      '                   Check ds.res4.senres().grad_order_no.\n\n']);
  data=[];ds=-1;  % Force an error in the calling program
  return
end

%  Check balance1 and put in a standard form.
balance1=check_balance(balance1);
if isempty(balance1) | isempty(balance0)
  fprintf(['\nsetCTFDataBalance: size(ds balance) =[%d %d] (after call to ',...
      'getDsBalance.)\n',...
      '                size(new balance)=[%d %d] (after call to check_balance.)\n\n'],...
    size(balance0),size(balance1));
  data=[];ds=-1;  % Force an error in the calling program
  return
end

%  Print a message?
chanset=strvcat('MEGs','Grefs');
for k=1:2
  if messages==2 | (messages==1 & ~strcmp(balance0(k,:),balance1(k,:)))
    fprintf('setCTFDataBalance:  Adjusting %s from balance=%s to balance=%s.\n',...
      chanset(k,:),balance0(k,:),balance1(k,:));
  end
end
clear k chanset;

if isequal(balance0,balance1)
  return
end


% Convert from balance0 to unbalanced.
if ~isequal(balance0,strvcat('NONE','NONE'));
  data=unbalance_data(data,ds,balance0,unit,nChBlock,chanList,badGref);
  % unbalance_data returns data=[] when list badGref includes a Gref that is required for
  % balancing.
  if isempty(data);ds=-1;return;end
end

% Convert from unbalanced to balance1.
if ~isequal(balance1,strvcat('NONE','NONE'));
  % balance_data returns data=[] when badGref includes a Gref that is required for
  % balancing.
  data=balance_data(data,ds,balance1,unit,nChBlock,chanList,badGref,messages);
  if isempty(data);ds=-1;return;end
end

%  Mark the balance order in structure ds.
%  Make sure that the bad channels are marked with a different order.  This will force
%  errors if these data are used later.  Set the data in the bad channels to zero, and add
%  the bad channels to ds.BadChannels.
Gorder=strmatch(balance1(2,:),balanceOptions)-1;
MEGorder=strmatch(balance1(1,:),balanceOptions)-1;
for k=setdiff(Greflist,badGref)
  ds.res4.senres(k).grad_order_no=Gorder;
end
for k=setdiff(GMEGlist,badGMEG)
  ds.res4.senres(k).grad_order_no=MEGorder;
end
for k=badGref
  ds.res4.senres(k).grad_order_no=0;
end
for k=badGMEG
  ds.res4.senres(k).grad_order_no=round(3-MEGorder);
end
%  Set channels with bad balance order parameter to zero.
kBad=[badGref badGMEG];
if ~isempty(chanList)
  kBad=reindex(chanList,intersect(chanList,kBad));
end
for k=kBad
  if strcmp(prec,'single')
    data(:,k,:)=single(zeros(size(data,1),1,size(data,3)));
  else
    data(:,k,:)=zeros(size(data,1),1,size(data,3));
  end
end
clear k Gorder MEGorder;

if isempty([badGref badGMEG]);return;end
% Add bad channels to ds.BadChannels

if ~isfield(ds,'BadChannels');ds.BadChannels=char([]);end
BadList=union(badGref,badGMEG);
if ~isempty(ds.BadChannels)
  for k=1:size(ds.BadChannels,1)
    BadList=[BadList strmatch(deblank(ds.BadChannels(k,:)),ds.res4.chanNames)];
  end
  ds.BadChannels=char([]);
end
for k=sort(BadList)
  ds.BadChannels=strvcat(ds.BadChannels,strtok(ds.res4.chanNames(k,:),'- '));
end
return
% ************** End of function setCTFDataBalance **********************************
% **********************************************************************************

%  ********************************************************************************
% ******************* function getDsBalance **************************************
function [dsbalance,badGref,badGMEG,Greflist,GMEGlist]=getDsBalance(ds);
%  Get the current balance state from structure ds.  Check for channels that are marked
%  as bad by having a different grad_order_no.
%  Check whether things make sense.

%  Return the balance and lists of bad gradiometers.
%  Also return the complete list of reference gradiometers and MEG gradiometers since
%  setCTFDataBalance uses them.

%  In normal operation, getDsBalance will return badGref=[], badGMEG=[];
balanceOptions=strvcat('NONE','G1BR','G2BR','G3BR');
Greflist=[];  % List of reference gradiometers (sensorTypeIndex=1)
Greforder=[];
GMEGlist=[];  % List of MEG gradiometers (sensorTypeIndex=5)
GMEGorder=[];

% Make lists of sensors and gradient balancing.
for k=1:ds.res4.no_channels
  if ds.res4.senres(k).sensorTypeIndex==1
    Greflist=[Greflist k];
    Greforder=[Greforder ds.res4.senres(k).grad_order_no];
  elseif ds.res4.senres(k).sensorTypeIndex==5
    GMEGlist=[GMEGlist k];
    GMEGorder=[GMEGorder ds.res4.senres(k).grad_order_no];
  end   
end

% Reference balance OK?
dsbalance=char(zeros(2,4));
if any(Greforder<0) | any(Greforder>1)
  fprintf(['\nsetCTFDataBalance:  The reference gradiometer balance orders make no ',...
      'sense.\n            %2d Gref channels.  %2d have grad_order_no=0\n'],...
    length(Greflist),sum(Greforder==0));
  fprintf('                               %2d have grad_order_no=1\n',sum(Greforder==1));
  fprintf('                               %2d have grad_order_no<0 or >1\n',...
    sum(Greforder>1)+sum(Greforder<0));
  dsbalance=strvcat(dsbalance(1,:),'BAD');
  badGref=Greflist;  % Mark everything bad
elseif sum(Greforder==0)>sum(Greforder==1)
  dsbalance=strvcat(dsbalance(1,:),balanceOptions(1,:));
  badGref=Greflist(find(Greforder~=0));
else
  dsbalance=strvcat(dsbalance(1,:),balanceOptions(2,:));
  badGref=Greflist(find(Greforder~=1));
end

%  Sort MEG gradiometer balance.  In correct operation, all MEG channels should have the
%  same balance (0,1,2 or 3).
for bal=0:3
  nMEGbal(bal+1)=sum(GMEGorder==bal);
end

[maxbal MEGorder]=max(nMEGbal);
MEGorder=MEGorder-1;

if maxbal>ceil(0.5*length(GMEGlist)) & all(GMEGorder>=0) & all(GMEGorder<=3)
  dsbalance(1,:)=balanceOptions(MEGorder+1,:);
  badGMEG=GMEGlist(GMEGorder~=MEGorder);
else
  fprintf('\nsetCTFDataBalance:  The MEG-sensor balance orders make no sense.\n');
  fprintf('            %3d MEG gradiometer channels.\n',length(GMEGlist));
  fprintf('            %3d have grad_order_no=%d\n',[nMEGbal;0:3]);
  if sum(nMEGbal)<length(GMEGlist)
    fprintf('            %3d have grad_order_no<0 or >3\n\n',...
      length(GMEGlist)-sum(nMEGbal));
  end
  dsbalance=strvcat('BAD',dsbalance(2,:));
  badGMEG=GMEGlist;  % Mark everything bad`
end

return
% ************** End of function getDsBalance **********************************
% **********************************************************************************

%  ********************************************************************************
% ******************* function check_balance **************************************
function balance=check_balance(balance);

% make sure that character array balance has the correct format.

balance_options=strvcat('NONE','G1BR','G2BR','G3BR');

% Make balance into a [4 2] character array
if ~ischar(balance) & ~isempty(balance)
  fprintf('setCTFDataBalance (check_balance): balance is not a character array.\n');
  balance=char([]);     % Force an error.
  return;
elseif isempty(deblank(balance))
  balance=strvcat('NONE','NONE');
  return
elseif size(balance,1)==1
  balance=strvcat(balance,'NONE');  % No Gref balancing
elseif size(balance,1)>2
  fprintf('setCTFDataBalance (check_balance): size(balance)=[%d %d]?\n',size(balance));
  balance=char([]);  % Force an error.
  return
end

balance=upper(balance);
if size(balance,2)>4;balance=balance(:,1:4);end

for k=1:2  % k=1:MEGs, k=2:Grefs
  if isempty(deblank(balance(k,:)));
    balance(k,:)='NONE';
  elseif isempty(strmatch(balance(k,:),balance_options(1:(6-2*k),:)))
    fprintf('check_balance: balance(%d,:)=%s   Not an allowed option.\n',k,balance(k,:));
    balance=char([]);  % Force an error?
    return
  end
end

return
% ************** End of function check_balance **********************************
% **********************************************************************************



% ************************************************************************
%  ****** Convert from balanced to unbalanced data  **********************
function data=unbalance_data(data,ds,balance0,unit,nChBlock,chanList,badGref);
%  Inputs: data:  Data array.  size(data)=[Npt Nchan Ntrial]
%          ds:    ds strtcuture for the data. 
%          balance0: balance0(1,:) = MEG  balance before unbalancing
%                    balance0(2,:) = Gref balance before unbalancing
%          unit : 'ft','t' (physical units)
%                 'int', 'phi0'
%          chanList: List if channels included in array data (referred to the channel
%                   numbering in structure ds)
%          badGref: A list of gRef channels markes as bad because they seen to have the
%                   wrong grad_order_no.

if ~exist('chanList');chanList=[];end
prec=class(data);

%  Read the balance tables for both MEG and reference gradiometers.
%  Usually balance0(2,:)='NONE', and it returns alphaGref=[].
[alphaMEG,MEGindex,MEGbalanceindex,alphaGref,Grefindex,Gbalanceindex]=...
  getCTFBalanceCoefs(ds,balance0,unit);

%  Unbalance the Grefs first.  In almost all cases, balance0(2,:)='NONE',
%  so the program will skip this section.
if ~strcmp(upper(balance0(2,:)),'NONE')
  %  The index lists give channel numbers referred to data set channel numbering

  if isempty(Grefindex) | isempty(Gbalanceindex)
    fprintf('setCTFDataBalance (unbalance_data): balance=%s\n',...
      '\t\tsize(Grefindex)=[%d %d]  size(Gbalanceindex)=[%d %d] ??\n',balance0(2,:),...
      size(Grefindex),size(Gbalanceindex));
    data=[];  % Force an error in the calling program
    return
  end
  
  % Are there any badGref in the Gbalanceindex list?
  if ~isempty(intersect(Gbalanceindex,badGref))
    fprintf(['setCTFDataBalance (unbalance_data): A reference gradiometer marked bad ',...
        'because of its grad_order_no\n',...
        '                 appears in the list of reference sensors required for ',...
        '%s Gref balancing.\n'],balance0(2,:));
    fprintf('        Gbalanceindex=');fprintf(' %d',Gbalanceindex);fprintf('\n');
    fprintf('        badGref=');fprintf(' %d',badGref);fprintf('\n');
    data=[];  % Force an error in the calling program
    return
  end
  
  if isequal(chanList,1:max(max(Gbalanceindex),max(Grefindex)));chanList=[];end
  if ~isempty(chanList)
    if ~isempty(setdiff(Gbalanceindex,chanList))
      fprintf(['setCTFDataBalance (unbalance_data): List chanList does not include ',...
          'all of the reference sensors in the %s table for reference gradiometers.\n'],...
        balance0(2,:));
      data=[];
      return
    end
    [Gbalanceindex,Grefindex,alphaGref]=...
      reindex(chanList,Gbalanceindex,Grefindex,alphaGref);
  end
  
  %  Balancing reference gradiometers: 
  %    balanced_data(:,Grefindex)
  %                 =raw_data(:,Grefindex)-raw_data(:,Gbalanceindex)*alphaGref
  %     size(alphaGref)=[length(Gbalanceindex) length(Grefindex)]
  %  Note : not all of the sensors in list Gbalanceindex are gradiometers.
  
  %  Rewrite as a matrix equation:
  %      balanced_data(:,Gindex)=raw_data(:,Gindex)*Gbalmat
  %      Gindex=list of all the sensors involved in balancing the reference
  %             gradiometers.
  %      size(Gbalmat)=[length(Gindex) length(Gindex)]
  
  [Gindex]=union(Gbalanceindex,Grefindex);
  [Cx,Fbal]=intersect(Gindex,Gbalanceindex);  % Gindex(Fbal)=Gbalanceindex
  [Cx,Fref]=intersect(Gindex,Grefindex);      % Gindex(Fref)=Grefindex
  
  Gbalmat=eye(length(Gindex));
  Gbalmat(Fbal,Fref)=Gbalmat(Fbal,Fref)-alphaGref;
  clear Fbal Fref Cx Grefindex Gbalanceindex alphaGref;
  %  Convert to unbalanced reference gradiometers
  for pt=1:size(data,3)
    if strcmp(prec,'single')
      data(:,Gindex,pt)=single(double(data(:,Gindex,pt))/Gbalmat);
    else
      data(:,Gindex,pt)=data(:,Gindex,pt)/Gbalmat;
    end
  end
  clear Gbalmat pt;
end  %  Finished unbalancing the reference gradiometers

if strcmp(upper(balance0(1,:)),'NONE');
  return;  % No unbalancing required for MEG gradiometers
end

% Are there any badGref in the MEGbalanceindex list?
if ~isempty(intersect(MEGbalanceindex,badGref))
  fprintf(['setCTFDataBalance (unbalance_data): A reference gradiometer marked bad ',...
      'because of its grad_order_no\n',...
      '                 appears in the list of reference sensors required for ',...
      '%s MEG balancing.\n'],balance0(1,:));
  fprintf('        MEGbalanceindex=');fprintf(' %d',MEGbalanceindex);fprintf('\n');
  fprintf('        badGref=');fprintf(' %d',badGref);fprintf('\n');
  data=[];  % Force an error in the calling program
  return
end

%  Don't bother with chanList if it obviously includes all the SQUID channels and in
%  the dataset ordering.
if isequal(chanList,1:max(max(MEGindex),max(MEGbalanceindex)));chanList=[];end

if isempty(alphaMEG) | isempty(MEGindex) | isempty(MEGbalanceindex)
  return;
elseif ~isequal(size(alphaMEG),[length(MEGbalanceindex) length(MEGindex)])
  fprintf(['setCTFDataBalance (unbalance_data): size(alphaMEG)=[%d %d] ',...
      'length(MEGbalanceindex)=%d length(MEGindex)=%d\n'],...
    size(alphaMEG),length(MEGbalanceindex),length(MEGindex));
  data=[];
  return
elseif all(reshape(alphaMEG,1,length(MEGbalanceindex)*length(MEGindex))==0)
  return
end

%  Make the index lists refer to the entries in chanList.
if ~isempty(chanList)
  if ~isempty(setdiff(MEGbalanceindex,chanList))
    fprintf(['setCTFDataBalance (unbalance_data): List chanList does not include ',...
        'all of the reference sensors in the %s table for MEG gradiometers.\n'],...
      balance0(1,:));
    data=[];
    return
  end
  [MEGbalanceindex,MEGindex,alphaMEG]=...
    reindex(chanList,MEGbalanceindex,MEGindex,alphaMEG);
end

%  Reverse the balance applied to the MEG channels. Needs unbalanced reference data.
if ~isempty(alphaMEG) & ...
    ~isequal(alphaMEG,zeros(length(MEGbalanceindex),length(MEGindex)))
  %  Unbalance MEG data trial-by-trial and in blocks of nChBlock channels
  for pt=1:size(data,3)
    for m=1:nChBlock:length(MEGindex)
      mptr=m:min(m+nChBlock-1,length(MEGindex));
      MEGblock=MEGindex(mptr);
      if strcmp(prec,'single');
        data(:,MEGblock,pt)=single(double(data(:,MEGblock,pt))+...
          double(data(:,MEGbalanceindex,pt))*alphaMEG(:,mptr));
      else
        data(:,MEGblock,pt)=data(:,MEGblock,pt)+...
          data(:,MEGbalanceindex,pt)*alphaMEG(:,mptr);
      end
    end
  end
  clear pt m mptr MEGblock;
end

return
% *******************   End of function unbalance_data   ****************************
% ***********************************************************************************



% ***********************************************************************************
% ***************  Convert from unbalanced data to balanced data  *******************
function  data=balance_data(data,ds,balance,unit,nChBlock,chanList,badGref,messages);

% Extracts the balance table from structure ds and balances data as specified in balance.
% If the balancing requires a reference gradiometer marked in the badGref list, an error 
% message is printed and the data array is set to [].

%  If a MEG channel is missing from the balance table, then the channel is set to zero.

%  Inputs:  data: The aray of data to be balanced. 
%           ds:   The ds structure read by readCTFds that contains the balance table.
%           balance: The desired output balance state.
%                      balance(1,:) = MEG  balance
%                      balance(2,:) = Gref balance
%           unit : Data units ('fT','T','phio','int')
%          nChBlock : rad the data in blocks of nChBlock to reduce memory requirements.
%          chanList: List of the channels that actual;ly appear in data.  referred to the
%                     dataset channel numbering (not to te list of SQUID channels)
%          badGref: List of reference gradiometer channels that cannot be used for
%                   balancing.

if ~exist('chanList');chanList=[];end
prec=class(data);
[alphaMEG,MEGindex,MEGbalanceindex,alphaGref,Grefindex,Gbalanceindex]=...
  getCTFBalanceCoefs(ds,balance,unit);
if ~strcmp(lower(balance(1,:)),'NONE');
  %  Are there any MEG channels missing from list MEGindex (i.e. missing from the balance table)?
  %  Make a list of MEG channels referred to dataset channel numbering.
  MEGlist=[];
  for k=1:ds.res4.no_channels
    if ds.res4.senres(k).sensorTypeIndex==5;MEGlist=[MEGlist k];end
  end
  missingMEG=setdiff(MEGlist,MEGindex);
  if ~isempty(missingMEG) & messages~=0 & ...
      (isempty(chanList) | intersect(missingMEG,chanList))
    fprintf('setCTFDataBalance (balance_data): Channel(s) missing from the balance table:\n');
    for k=missingMEG
      fprintf('\t%3d: %s\n',k,strtok(ds.res4.chanNames(k,:),'- '));
    end
  end
  clear MEGlist k;
  
  % Are there any badGref in the MEGbalanceindex list?
  if ~isempty(intersect(MEGbalanceindex,badGref))
    fprintf(['setCTFDataBalance (unbalance_data): A reference gradiometer marked bad ',...
        'because of its grad_order_no\n',...
        '                 appears in the list of reference sensors required for ',...
        '%s MEG balancing.\n'],balance(1,:));
    fprintf('        MEGbalanceindex=');fprintf(' %d',MEGbalanceindex);fprintf('\n');
    fprintf('        badGref=');fprintf(' %d',badGref);fprintf('\n');
    data=[];  % Force an error in the calling program
    return
  end
  
  if isequal(chanList,1:max(max(MEGindex),max(MEGbalanceindex)));
    chanList=[];
  end
  %  Check if alphaMEG is all zeros or is improperly defined.
  do_balance=~isempty(alphaMEG) & ~isempty(MEGindex) & ~isempty(MEGbalanceindex);
  if do_balance
    if ~isequal(size(alphaMEG),[length(MEGbalanceindex) length(MEGindex)]);
      fprintf(['setCTFDataBalance (balance_data): size(alphaMEG)=[%d %d] ',...
          'length(MEGbalanceindex)=%d length(MEGindex)=%d\n'],...
        size(alphaMEG),length(MEGbalanceindex),length(MEGindex));
      data=[];
      return
    elseif isempty(chanList) & size(data,2)<max(MEGindex)
      fprintf(['setCTFDataBalance (balance_data): chanList=[], but size(data,2)=%d  ',...
          'max(MEGindex)=%d ?\n'],size(data,2),max(MEGindex));
      data=[];
      return
    end
    do_balance=~all(reshape(alphaMEG,1,length(MEGbalanceindex)*length(MEGindex))==0);
  end
  
  if do_balance
    if ~isempty(chanList)
      % Re-sort the reference sensor indices and make them refer to chanList
      [MEGbalanceindex,MEGindex,alphaMEG]=...
        reindex(chanList,MEGbalanceindex,MEGindex,alphaMEG);
      %  Refer the missing MEG channels to chanList
      if ~isempty(missingMEG);missingMEG=reindex(chanList,missingMEG);end
    end
    %  Balance MEGs data trial-by-trial and in blocks of nChBlock channels
    for pt=1:size(data,3)
      for m=1:nChBlock:length(MEGindex)
        mptr=m:min(m+nChBlock-1,length(MEGindex));
        MEGblock=MEGindex(mptr);
        if strcmp(prec,'single')
          data(:,MEGblock,pt)=single(double(data(:,MEGblock,pt))-...
            double(data(:,MEGbalanceindex,pt))*alphaMEG(:,mptr));
        else
          data(:,MEGblock,pt)=data(:,MEGblock,pt)-...
            data(:,MEGbalanceindex,pt)*alphaMEG(:,mptr);
        end
      end
      %  Zero the channel with missing coefficients.  They cannot possibly be correct.
      for m=missingMEG
        if strcmp(prec,'single')
          data(:,m,pt)=single(zeros(size(data,1),1,1));
        else
          data(:,m,pt)=zeros(size(data,1),1,1);
        end
      end
    end
  end
  clear alphaMEG MEGindex MEGbalanceindex Ix pt m mptr MEGblock do_balance missingMEG;
end    %  Finished balancing the MEG data

%  Is Gref balancing requested?  In most cases, balance(2,:)='NONE'
if size(balance,1)==1;return;end
if strcmp(balance(2,:),'NONE');return;end

% Are there any bad Gref sensors in the Gbalanceindex list?
if ~isempty(intersect(Gbalanceindex,badGref))
  fprintf(['setCTFDataBalance (balance_data): A reference gradiometer marked bad ',...
      'because of its grad_order_no\n',...
      '                 appears in the list of reference sensors required for ',...
      '%s Gref balancing.\n'],balance(2,:));
  fprintf('        Gbalanceindex=');fprintf(' %d',Gbalanceindex);fprintf('\n');
  fprintf('        badGref=');fprintf(' %d',badGref);fprintf('\n');
  data=[];  % Force an error in the calling program
  return
end

%  Make the index lists refer to the entries in chanList.
if ~isempty(chanList)
  [Gbalanceindex,Grefindex,alphaGref]=...
    reindex(chanList,Gbalanceindex,Grefindex,alphaGref);
end

%  Balance the reference channels
for pt=1:size(data,3)
  for m=1:nChBlock:length(Grefindex)
    mptr=m:min(m+nChBlock-1,length(Grefindex));
    Grefblock=Grefindex(mptr);
    if strcmp(prec,'single');
      data(:,Grefblock,pt)=single(double(data(:,Grefblock,pt))-...
        double(data(:,Gbalanceindex,pt))*alphaGref(:,mptr));
    else
      data(:,Grefblock,pt)=data(:,Grefblock,pt)-...
        data(:,Gbalanceindex,pt)*alphaGref(:,mptr);
    end
  end
end
clear pt m mptr Grefblock Grefindex Gbalanceindex;

return
% ********************  End of function balance_data   ****************************
% ***********************************************************************************


% ***********************************************************************************
% **************  Reassign indices in sensor and reference lists.  ******************

function [refindex,sensorindex,alfa2]=reindex(chanlist,reflist,sensorlist,alfa);

%  Reindex lists reflist and sensor list, and reorders 2nd index of matrix alfa.

%  Inputs : Lists of sensor and references for balancing and matrix of balance
%           coefficients.
%  Outputs : index lists: chanlist(refindex)=reflist
%                         chanlist(sensorindex)=sensorlist
%                         alfa2=The part of alfa that refers to the part of sensorlist
%                               that appears in chanlist.

refindex=[];sensorindex=[];alfa2=[];  % Force an error on an early return.

% Sensible argument lists?
if nargout>nargin-1 | nargout==0
  fprintf('reindex: %d outputs, but only %d inputs.\n',nargout,nargin);
  return
end

[X c1 c2]=intersect(chanlist,reflist);
[Y ylist]=sort(c2);
refindex=c1(ylist);

if length(refindex)~=length(reflist)
  fprintf('setCTFDataBalance (reindex): length(refindex)=%d  length(reflist)=%d\n',...
    length(refindex),length(reflist));
  fprintf('                     Some references are not included in chanlist.\n');
  return
end

if nargin>=3 & nargout>=2
  [X c1 c2]=intersect(chanlist,sensorlist);
  [Y ylist]=sort(c2);
  sensorindex=c1(ylist);
end

if nargin>=4 & nargout>=3
  alfa2=alfa(:,Y);
end
return

% ***********************************************************************************
