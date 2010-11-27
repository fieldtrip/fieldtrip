function [alphaMEG,MEGindex,MEGbalanceindex,alphaGref,Grefindex,Gbalanceindex]=...
  getCTFBalanceCoefs(ds,balanceType,unit);

%  Reads balance coefficients for SQUID data in phi0's, and converts to coefficients
%  for data in physical units.

%  Input : ds : (1) ds structure from readCTFds
%           balanceType : 'NONE'  : No balancing
%                          'G1BR'  : 1st order balancing
%                          'G2BR'  : 2nd order balancing
%                          'G3BR'  : 3rd order balancing
%                          'G3AR'  : 3rd order balancing + adaptive
%                  If only MEG balance table is requested, size(balanceType)=[1 5]
%                  If Gref balance table is requested also, size(balanceType)=[2 4]
%          unit:  Type of units.  Option: 'fT','T','phi0','int'.
%                 If unit is not entered, of unit=[], the default is unit='fT'.

%  Outputs : 
%     alphaMEG,MEGindex,MEGbalanceindex : MEG balancing coefficients for data in fT or T.
%            Suppose data array returned by getTrial2 has size(data)=[npt nchan].
%            nMEG = number of MEG channels in the data set.   
%            MEG channels have sensorTypeIndex==5.
%            
%            MEGindex = list of MEG channels referred to the channel numbering in
%                       the complete dataset (i.e. it is not referred to only the list of
%                       SQUID sensors).
%            by looking for sensors with ds.res4.senres.sensorTypeIndex==5
%                      size(MEGindex)=[1 nMEG]
%                      MEG data = data(:,MEGindex)

%            MEGbalanceindex = list of reference channels for MEG balancing.   
%                size(MEGbalanceindex)=[1 nRef]
%                Reference data for MEG balancing = data(:,MEGbalanceindex)

%            alphaMEG = balance coefficients.  size(alphaMEG)=[nRef nMEG]
%
%            Balancing MEG data : 
%              - Start with data AFTER converting to physical units.
%                (I.e. SQUID data must be in fT or T.)
%              - balanced_data(:,MEGindex)=data(:,MEGindex)-data(:,MEGbalanceindex)*alphaMEG

%            If user specifies balanceType=[], ' ' or 'none', then getCTFBalanceCoefs returns 
%            alphaMEG=zeros(0,nMEG), MEGindex=list of MEG sensor channels,  
%            MEGbalanceindex=[], alphaGref=zeros(0,nGef), Grefindex=list of Gref channels
%            and Gbalanceindex=[];

%   alphaGref,Grefindex,Gbalanceindex : Reference-gradiometer balancing coeficients.
%            These output arrays are optional.  If the calling statement includes
%            them, then this program extracts a table of balance coefficients for
%            the reference gradiometers from the G1BR table in the coefficient files.
%               nGref = no. of reference gradiometers channels in the data set
%                       (sensorTypeIndex=1)
%               Grefindex = list of reference channels.   size(Grefindex)=[1 nGref]
%                           Gradient reference data = data(:,Grefindex)
%               Gbalanceindex = list of channels for balancing reference 
%                                 gradiometers.  size(Gbalanceindex)=[1 nGbalcoef]
%               Balancing reference data = data(:,Gbalanceindex)

%               alphaGref = balance coefficients for ref. grads.  
%                          size(alphaGref)=[nGbalcoef nGref]
%
%               Balancing ref. gradiometer data :
%              - Use data AFTER converting to physical units. (I.e. data in fT or T.)
%                  balanced_data(:,Grefindex)=data(:,Grefindex)-data(:,Gbalanceindex)*alphaGref

%  Calls function getRawCTFBalanceCoefs (included in this listing).

if nargout==0 & nargin==0
  fprintf(['\ngetCTFBalanceCoefs: Version 1.1   17 April 2007   ',...
      'Reads balance coefficients from ds.res4.scrr.\n\n',...
      '\tMEG balancing : [alphaMEG,MEGindex,MEGbalanceindex]=',...
      'getCTFBalanceCoefs(ds,balanceType,unit);\n\n',...
      '\tMEG & Gref balancing : [alphaMEG,MEGindex,MEGbalanceindex,',...
      'alphaGref,Grefindex,Gbalanceindex]=\n',...
      '\t                                                    ',...
      'getCTFBalanceCoefs(ds,balanceType,unit);\n\n']);
  return
end

balanceOptions=strvcat('NONE','G1BR','G2BR','G3BR', 'G3AR');
balanceOptionsEnds = [size(balanceOptions, 1), 2]; % Which options are available for MEGs and Grefs.
physical_options=strvcat('fT','T');
raw_options=strvcat('phi0','int');
unit_options=strvcat(physical_options,raw_options);
default_unit='fT';

%  Specify outputs in case of an early return due to an error.
MEGindex=[];
Grefindex=[];
alphaMEG=[];
MEGbalanceindex=[];
alphaGref=[];
Gbalanceindex=[];

% Check that the inputs are sensible.
if nargin<2
  fprintf(['\ngetCTFBalanceCoefs: Only %d input arguments?  ',...
      'Must specify at least ds and balanceType.\n\n'],nargin);
  return
elseif ~isstruct(ds) | isempty(ds) | ~ischar(balanceType) | isempty(balanceType)
  fprintf('\ngetCTFBalanceCoefs: Wrong argument types or sizes.\n\n');
  whos ds balanceType
  return
elseif ~isfield(ds,'res4')
  fprintf('\ngetCTFBalanceCoefs: Field res4 is missing from structure ds.\n\n');
  ds
  return
elseif size(balanceType,1)>2
  fprintf('\ngetCTFBalanceCoefs: size(balanceType)=[');fprintf(' %d',size(balanceType));...
    fprintf(']   Must be[1 4] or [2 4].\n\n');
  return
end
%  Must have 3 or 6 output arguments.  Set balanceType(2,:)='NONE' if necessary.
if ~any(nargout==[3 6]);
  fprintf(['\ngetCTFBalanceCoefs: Called with %d output arguments.  ',...
    'Must be 3 or 6.\n\n'],nargout);
  return
elseif (nargout==3 & size(balanceType,1)>1) | (nargout==6 & size(balanceType,1)==1)
  balanceType=strvcat(deblank(balanceType(1,:)),'NONE');
end

%  At this point, size(balanceType,1)=2.

%  Check that balanceType has allowed values
for k=1:size(balanceType,1)  % k=1:MEGs, k=2:Grefs
  if isempty(strmatch(balanceType(k,:),balanceOptions(1:balanceOptionsEnds(k),:)))
    fprintf('\ngetCTFBalanceCoefs: balanceType(%d,:)=%s   Not an allowed option.\n\n',...
      k,balanceType(k,:));
    return
  end
end

%  Check the data units, convert to lower case and flag incorrect unit specification
if ~exist('unit');
  unit=default_unit;
elseif isempty(unit);
  unit=default_unit;
elseif ~ischar(unit)
  fprintf('\ngetCTFBalanceCoefs: Input unit has the wrong type: class(unit)=%s\n\n',...
    class(unit));
  return
end
unit=lower(unit);
if isempty(strmatch(unit,lower(strvcat(physical_options,raw_options))))
  fprintf('\ngetCTFBalanceCoefs: unit=%s.   Must be one of ',unit);
  for k=1:size(unit_options,1);fprintf(' %s,',deblank(unit_options(k,:)));end;
  fprintf('\n\n');
  return
end
physical=~isempty(strmatch(unit,lower(physical_options)));

balanceType=upper(deblank(balanceType));

[betaMEG,MEGindex,MEGbalanceindex,betaGref,Grefindex,Gbalanceindex]=...
  getRawCTFBalanceCoefs(ds,balanceType);

% No balancing is requested, just return lists of MEGindex and Grefindex.
if isempty(MEGbalanceindex)
  alphaMEG=zeros(0,length(MEGindex));
  MEGbalanceindex=[]; % force size=[0 0]
end
if isempty(Gbalanceindex)
  alphaGref=zeros(0,length(Grefindex));
  Gbalanceindex=[]; % force size=[0 0]
end
if isempty(MEGbalanceindex) & isempty(Gbalanceindex)
  return
end

%  betaMEG and betaGref are the balance coefficients when signals are in phi0's.
%  Convert to balance coefficients when signals are in physical units (T or fT).
%  invproperGain is introduced to take care of situations where bad channels 
%  are labelled by setting their gain to zero.

if physical
  properGain=zeros(1,ds.res4.no_channels);
  invproperGain=zeros(1,ds.res4.no_channels);
  for k=1:ds.res4.no_channels
    if any(ds.res4.senres(k).sensorTypeIndex==[0 1 5:7])  % SQUIDs only
      properGain(k)=ds.res4.senres(k).properGain;     % properGain = phi0/T
      if properGain(k)~=0
        invproperGain(k)=1/properGain(k);
      end
    end
  end
end

if ~isempty(MEGbalanceindex)
  if physical
    alphaMEG=betaMEG.*(properGain(MEGbalanceindex)'*invproperGain(MEGindex));
  else
    alphaMEG=betaMEG;
  end
end

if ~isempty(Gbalanceindex)
  if physical
    alphaGref=betaGref.*(properGain(Gbalanceindex)'*invproperGain(Grefindex));
  else
    alphaGref=betaGref;
  end
end
return
%%%%%%%%%%% End of  getBalanceCoef  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Function getRawCTFBalanceCoefs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [betaMEG,MEGindex,Refindex,betaGref,Grefindex,Gbalanceindex]=...
  getRawCTFBalanceCoefs(ds,balanceType);

%  getRawCTFBalanceCoefs.  Extracts raw-data (i.e. integer data) gradiometer balancing
%                          coefficients from structure array ds.res4.scrr (ds is produced
%                          by readCTFds).

%  Date : 24 OCt 2006
%  Author : Harold Wilson

%  Inputs : ds   : Structure returned by readCTFds.
%           balanceType : 'NONE'  : No balancing
%                          'G1BR'  : 1st order balancing
%                          'G2BR'  : 2nd order balancing
%                          'G3BR'  : 3rd order balancing
%                          'G3AR'  : 3rd order balancing + adaptive
%                  If only MEG balance table is requested, size(balanceType)=[1 4]
%                  If Gref balance table is requested also, size(balanceType)=[2 4]
%           If balancing 'NONE' is specified, getRawCTFBalanceCoefs returns lists MEGindex and
%           Grefindex, and sets Refindex=[], Gbalanceindex=[].
%           The reference gradiometers have only G1BR balancing coefficients.

%  Outputs : 
%     betaMEG,MEGindex,Refindex : MEG balancing coefficients.
%            Suppose data array returned by getCTFdata has size(data)=[npt nchan].
%            nMEG = number of MEG channels in the data set.   
%            MEG channels have sensorTypeIndex==5.
%            
%            MEGindex = list of MEG channels in the data set.  It is assembled
%            by looking for sensors with ds.res4.senres.sensorTypeIndex==5
%                      size(MEGindex)=[1 nMEG]
%                      MEG data = data(:,MEGindex)

%            Refindex = list of reference channels for MEG balancing.   
%                size(Refindex)=[1 nRef]
%                Reference data = data(:,Refindex)

%            betaMEG = balance coefficients.  size(betaMEG)=[nRef nMEG]
%
%            Balancing MEG data : 
%              - Start with data BEFORE converting to physical units.
%                (I.e. SQUID data must be in phi0's or raw integers.)
%              - balanced_data(:,MEGindex)=data(:,MEGindex)-data(:,Refindex)*betaMEG

%            If user specifies balanceType=[] or ' ', then getRawCTFBalanceCoefs returns 
%            betaMEG=zeros(1,nMEG),Refindex=3 and nRef=1;


%   betaGref,Grefindex,Gbalanceindex : Reference-gradiometer balancing coeficients.
%            These output arrays are optional.  If the calling statement includes
%            them, then this program extracts a table of balance coefficients for
%            the reference gradiometers from the G1BR table in the coefficient files.
%               nGref = no. of reference gradiometers channels in the data set
%                       (sensorTypeIndex=1)
%               Grefindex = list of reference channels.   size(Grefindex)=[1 nGref]
%                           Gradient reference data = data(:,Grefindex)
%               Gbalanceindex = list of channels for balancing reference 
%                                 gradiometers.  size(Gbalanceindex)=[1 nGbalcoef]
%               Balancing reference data = data(:,Gbalanceindex)

%               betaGref = balance coefficients for ref. grads.  
%                          size(betaGref)=[nGbalcoef nGref]
%
%               Balancing ref. gradiometer data : 
%              balanced_data(:,Grefindex)=data(:,Grefindex)-data(:,Gbalanceindex)*betaGref

%  No function calls.
missingMEGMessage=0;
missingGrefMessage=0;
balanceOptions=strvcat('NONE','G1BR','G2BR','G3BR', 'G3AR');
balanceOptionsEnds = [size(balanceOptions, 1), 2]; % Which options are available for MEGs and Grefs.
common_mode_only=0;

Brefindex=[];  %  Index list of reference magnetometers in the data arrays
Grefindex=[];  %  Index list of reference gradiometers in the data arrays
MEGindex=[];   %  Index list of MEG sensors in the data arrays
Gbalanceindex=[]; % Index list of sensors used as references to balance the reference
%                   gradiometers 
Refindex=[];   %  Index list of sensors used to balance the MEG sensors 
betaMEG=[];
betaGref=[];

% Check that the inputs are sensible.
if nargin<2
  fprintf(['\ngetRawCTFBalanceCoefs: Only %d input arguments?  ',...
      'Must specify at least ds and balance.\n\n'],nargin);
  return
elseif ~isstruct(ds) | isempty(ds) | ~ischar(balanceType) | isempty(balanceType)
  fprintf('\ngetRawCTFBalanceCoefs: Wrong argument types or sizes.\n\n');
  whos ds balanceType
  return
elseif ~isfield(ds,'res4')
  fprintf('\ngetRawCTFBalanceCoefs: Field res4 is missing from structure ds.\n\n');
  ds
  return
end

%  Check that the output list is OK.
if nargout~=3 & nargout~=6
  fprintf('\ngetRawCTFBalanceCoefs : Call specifies %d output arguments.\n',nargout);
  fprintf('\t\tMust have 3 output arguments (MEG balance coefficients)\n');
  fprintf('\t\tor 6 output arguments (reference gradiometer balancing also).\n\n');
  return
elseif nargout==3
  if size(balanceType,1)>1;balanceType=balanceType(1,:);end
else
  if size(balanceType,1)==1;balanceType=strvcat(balanceType,'NONE');end
end

%  Check that balanceType has allowed values
for k=1:size(balanceType,1)  % k=1:MEGs, k=2:Grefs
  if isempty(strmatch(balanceType(k,:),balanceOptions(1:balanceOptionsEnds(k),:)))
    fprintf('\ngetRawCTFBalanceCoefs: balance(%d,:)=%s   Not an allowed option.\n\n',...
      k,balanceType(k,:));
    return
  end
end

% Make lists of reference magnetometers, reference gradiometers and MEG sensors.
for q=1:length(ds.res4.senres)
  if ds.res4.senres(q).sensorTypeIndex==0
    Brefindex=[Brefindex q];
  elseif ds.res4.senres(q).sensorTypeIndex==1
    Grefindex=[Grefindex q];
  elseif ds.res4.senres(q).sensorTypeIndex==5  % Allow other MEG sensors?
    MEGindex=[MEGindex q];
  end
end
nBref=length(Brefindex);  % Don't currently use Brefindex or nBref
nGref=length(Grefindex);
nMEG=length(MEGindex);

nGbalcoef=0;  % Set to zero until we know the number by examining ds.res4.scrr
if nargout==6 & strcmp(balanceType(2,:),'NONE')
  Gbalanceindex=[];
  betaGref=zeros(0,nGref);
elseif nargout==6 & ~strcmp(balanceType(2,:),'NONE') 
  m1=1;         %  Get coefficients for balancing the reference gradiometers.
  mtot=size(ds.res4.scrr,2);
  for n=1:nGref
    Gname=strtok(ds.res4.chanNames(Grefindex(n),:),['- ',char(0)]);
    nGchar=length(Gname);
    for m=[m1:mtot 1:(m1-1)]
      if strncmp(Gname,char(ds.res4.scrr(m).sensorName),nGchar) & ...
          strcmp('G1BR',char(ds.res4.scrr(m).coefType));
        if nGbalcoef<=0  %  1st match.  Initialize table and get list of references
          nGbalcoef=ds.res4.scrr(m).numcoefs;
          betaGref=zeros(nGbalcoef,nGref);
          %  Assemble index array for balancing the reference gradiometers
          for q=1:nGbalcoef
            Refname=strtok(char(ds.res4.scrr(m).sensor(:,q))',['- ',char(0)]);
            pRef=strmatch(Refname,ds.res4.chanNames);
            if isempty(pRef)
              fprintf(['getRawCTFBalanceCoefs : Sensor %s appears in ',...
                  'ds.res4.scrr, but not in ds.res4.chanNames\n'],Refname);
              return
            end
            Gbalanceindex=[Gbalanceindex pRef];
          end
        end  %  end setup of balancing table for Ref. gradiometers
        if ds.res4.scrr(m).numcoefs~=nGbalcoef
          fprintf('\ngetRawCTFBalanceCoefs : %s has %d coefficients\n',...
            ds.res4.chanNames(Grefindex(1),1:nGchar),nGbalcoef);
          fprintf('                    %s  "  %d       "       ????\n\n',...
            ds.res4.chanNames(Grefindex(n),1:nGchar),ds.res4.scrr(m).numcoefs);
          betaGref=[];   % Force useless output.
          return
        end
        betaGref(:,n)=reshape(ds.res4.scrr(m).coefs(1:nGbalcoef),nGbalcoef,1);
        m1=m+1;
        break;  % Break out of m-loop.  Go to nect n value.
      end
      if (m==m1-1 & m1>1) | (m==mtot & m1==1)
        if missingGrefMessage==0
          if strncmp(balanceType(2,:), balanceOptions(5,:), 4)
            % Avoid warning for all sensors for adaptive coefficients.
            fprintf('\ngetRawCTFBalanceCoefs: Failed to find %s balance coefficients for reference sensors.\n',...
              balanceType(2,:));
          else
            fprintf(['\ngetRawCTFBalanceCoefs: Failed to find %s balance coefficients',...
              ' for sensor(s)'],balanceType(2,:));
            fprintf('\n\t\t\t\t');
          end
        end
        missingGrefMessage=missingGrefMessage+1;
        if ~strncmp(balanceType(2,:), balanceOptions(5,:), 4)
          if missingGrefMessage==10*round(missingGrefMessage/10)
            fprintf('\n\t\t\t\t');
          end
          fprintf(' %s',Gname);
        end
        betaGRef(:,n)=zeros(nGbalcoef,1);
        return
      end
    end  % End loop over m (searching for scrr(m).sensorName)
  end  % End loop over n (list of reference gradiometers)
end  % End of section getting coefficients to balance the reference gradiometers.
if missingGrefMessage>0;fprintf('\n');end
  
if strcmp(balanceType(1,:),'NONE')
  Refindex=[];
  betaMEG=zeros(0,nMEG);
  return
end

%  Get balance coefficients for the MEG sensors
nRef=0;
%  Pointers for search through ds.res4.scrr structure array.
m1=1;
mtot=size(ds.res4.scrr,2);

for n=1:nMEG
  MEGname=strtok(ds.res4.chanNames(MEGindex(n),:),['- ',char(0)]);
  nChar=length(MEGname);
  for m=[m1:mtot 1:(m1-1)]
    if strncmp(MEGname,char(ds.res4.scrr(m).sensorName),nChar) & ...
        strcmp(balanceType(1,:),char(ds.res4.scrr(m).coefType));
      if nRef<=0
        nRef=ds.res4.scrr(m).numcoefs;
        betaMEG=zeros(nRef,nMEG);
        for q=1:nRef   %  Assemble index array for balancing the MEG sensors
          Refname=strtok(char(ds.res4.scrr(m).sensor(:,q))',['- ',char(0)]);
          pRef=strmatch(Refname,ds.res4.chanNames);
          if isempty(pRef)
            fprintf(['\ngetRawCTFBalanceCoefs : Sensor %s appears in ',...
                'ds.res4.scrr, but not in ds.res4.chanNames\n\n'],Refname);
            return
          end
          Refindex=[Refindex pRef];
        end
      end   %  end setup of balancing table for MEG sensors
      if ds.res4.scrr(m).numcoefs~=nRef
        fprintf('\ngetRawCTFBalanceCoefs : %s - %s has %d coefficients\n',...
          ds.res4.chanNames(MEGindex(1),1:nChar),balanceType,nRef);
        fprintf('                    %s - %s "  %d       "       ????\n\n',...
          ds.res4.chanNames(MEGindex(n),1:nChar),balanceType,...
          ds.res4.scrr(m).numcoefs);
        betaMEG=[];   % An output that will force an error in the calling program.
        return
      end
      betaMEG(:,n)=reshape(ds.res4.scrr(m).coefs(1:nRef),nRef,1);
      m1=m+1;
      break;
    end
    if (m==m1-1 & m1>1) | (m==mtot & m1==1)
      if missingMEGMessage==0
        if strncmp(balanceType(2,:), balanceOptions(5,:), 4)
          % Avoid warning for all sensors for adaptive coefficients.
          fprintf('\ngetRawCTFBalanceCoefs: Failed to find %s balance coefficients for sensors.\n',...
            balanceType(2,:));
        else
          fprintf(['\ngetRawCTFBalanceCoefs: Failed to find %s balance coefficients',...
            ' for sensor(s)'],balanceType(1,:));
          fprintf('\n\t\t\t\t');
        end
      end
      missingMEGMessage=missingMEGMessage+1;
      if ~strncmp(balanceType(2,:), balanceOptions(5,:), 4)
        if missingMEGMessage==10*round(missingMEGMessage/10)
          fprintf('\n\t\t\t\t');
        end
        fprintf(' %s',MEGname);
      end
      betaMEG(:,n)=zeros(nRef,1);
    end
  end  % End of loop over m (ds.res4.scrr table)
end  % End of loop over MEG sensors
if missingMEGMessage>0;fprintf('\n');end

if common_mode_only
  if size(betaMEG,1)>3 & nMEG>0
    betaMEG=betaMEG(1:3,:);
    Refindex=Refindex(1:3);
  end
  if size(betaGref,1)>3 & nGref>0
    betaGref=betaGref(1:3,:);
    Gbalanceindex=Gbalanceindex(1:3);
  end
end

return
