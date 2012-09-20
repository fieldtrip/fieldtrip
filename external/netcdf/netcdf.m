function S = netcdf(File,varargin)
% Function to read NetCDF files
%   S = netcdf(File)
% Input Arguments
%   File = NetCDF file to read
% Optional Input Arguments:
%   'Var',Var - Read data for VarArray(Var), default [1:length(S.VarArray)]
%   'Rec',Rec - Read data for Record(Rec), default [1:S.NumRecs]
% Output Arguments:
%   S    = Structure of NetCDF data organised as per NetCDF definition
% Notes:
%   Only version 1, classic 32bit, NetCDF files are supported. By default
% data are extracted into the S.VarArray().Data field for all variables.
% To read the header only call S = netcdf(File,'Var',[]);
%
% SEE ALSO
% ---------------------------------------------------------------------------
S = [];

try
   if exist(File,'file') fp = fopen(File,'r','b');
   else fp = []; error('File not found'); end
   if fp == -1   error('Unable to open file'); end

% Read header
   Magic = fread(fp,4,'uint8=>char');
   if strcmp(Magic(1:3),'CDF') error('Not a NetCDF file'); end
   if uint8(Magic(4))~=1       error('Version not supported'); end
   S.NumRecs  = fread(fp,1,'uint32=>uint32');
   S.DimArray = DimArray(fp);
   S.AttArray = AttArray(fp);
   S.VarArray = VarArray(fp);

% Setup indexing to arrays and records
   Var = ones(1,length(S.VarArray));
   Rec = ones(1,S.NumRecs);
   for i = 1:2:length(varargin)
      if     strcmp(upper(varargin{i}),'VAR') Var=Var*0; Var(varargin{i+1})=1;
      elseif strcmp(upper(varargin{i}),'REC') Rec=Rec*0; Rec(varargin{i+1})=1;
      else error('Optional input argument not recognised'); end
   end
   if sum(Var)==0 fclose(fp); return; end

% Read non-record variables
   Dim = double(cat(2,S.DimArray.Dim));
   ID  = double(cat(2,S.VarArray.Type));

   for i = 1:length(S.VarArray)
      D = Dim(S.VarArray(i).DimID+1); N = prod(D); RecID{i}=find(D==0);
      if isempty(RecID{i})
         if length(D)==0 D = [1,1]; N = 1; elseif length(D)==1 D=[D,1]; end
         if Var(i)
            S.VarArray(i).Data = ReOrder(fread(fp,N,[Type(ID(i)),'=>',Type(ID(i))]),D);
            fread(fp,(Pad(N,ID(i))-N)*Size(ID(i)),'uint8=>uint8');
         else fseek(fp,Pad(N,ID(i))*Size(ID(i)),'cof'); end
      else S.VarArray(i).Data = []; end
   end

% Read record variables
   for k = 1:S.NumRecs
      for i = 1:length(S.VarArray)
         if ~isempty(RecID{i})
            D = Dim(S.VarArray(i).DimID+1); D(RecID{i}) = 1; N = prod(D);
            if length(D)==1 D=[D,1]; end
            if Var(i) & Rec(k)
               S.VarArray(i).Data = cat(RecID{i},S.VarArray(i).Data,...
                  ReOrder(fread(fp,N,[Type(ID(i)),'=>',Type(ID(i))]),D));
               if N > 1 fread(fp,(Pad(N,ID(i))-N)*Size(ID(i)),'uint8=>uint8'); end
            else fseek(fp,Pad(N,ID(i))*Size(ID(i)),'cof'); end
         end
      end
   end

   fclose(fp);
catch
   Err = lasterror; fprintf('%s\n',Err.message);
   if ~isempty(fp) && fp ~= -1 fclose(fp); end
end

% ---------------------------------------------------------------------------------------
% Utility functions

function S = Size(ID)
% Size of NetCDF data type, ID, in bytes
   S = subsref([1,1,2,4,4,8],struct('type','()','subs',{{ID}}));

function T = Type(ID)
% Matlab string for CDF data type, ID
   T = subsref({'int8','char','int16','int32','single','double'},...
               struct('type','{}','subs',{{ID}}));

function N = Pad(Num,ID)
% Number of elements to read after padding to 4 bytes for type ID
   N = (double(Num) + mod(4-double(Num)*Size(ID),4)/Size(ID)).*(Num~=0);

function S = String(fp)
% Read a CDF string; Size,[String,[Padding]]
   S = fread(fp,Pad(fread(fp,1,'uint32=>uint32'),1),'uint8=>char').';

function A = ReOrder(A,S)
% Rearrange CDF array A to size S with matlab ordering
   A = permute(reshape(A,fliplr(S)),fliplr(1:length(S)));

function S = DimArray(fp)
% Read DimArray into structure
   if fread(fp,1,'uint32=>uint32') == 10 % NC_DIMENSION
      for i = 1:fread(fp,1,'uint32=>uint32')
         S(i).Str = String(fp);
         S(i).Dim = fread(fp,1,'uint32=>uint32');
      end
   else fread(fp,1,'uint32=>uint32'); S = []; end

function S = AttArray(fp)
% Read AttArray into structure
   if fread(fp,1,'uint32=>uint32') == 12 % NC_ATTRIBUTE
      for i = 1:fread(fp,1,'uint32=>uint32')
         S(i).Str = String(fp);
         ID       = fread(fp,1,'uint32=>uint32');
         Num      = fread(fp,1,'uint32=>uint32');
         S(i).Val = fread(fp,Pad(Num,ID),[Type(ID),'=>',Type(ID)]).';
      end
   else fread(fp,1,'uint32=>uint32'); S = []; end

function S = VarArray(fp)
% Read VarArray into structure
   if fread(fp,1,'uint32=>uint32') == 11 % NC_VARIABLE
      for i = 1:fread(fp,1,'uint32=>uint32')
         S(i).Str      = String(fp);
         Num           = double(fread(fp,1,'uint32=>uint32'));
         S(i).DimID    = double(fread(fp,Num,'uint32=>uint32'));
         S(i).AttArray = AttArray(fp);
         S(i).Type     = fread(fp,1,'uint32=>uint32');
         S(i).VSize    = fread(fp,1,'uint32=>uint32');
         S(i).Begin    = fread(fp,1,'uint32=>uint32'); % Classic 32 bit format only
      end
   else fread(fp,1,'uint32=>uint32'); S = []; end
