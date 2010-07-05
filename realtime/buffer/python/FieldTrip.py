"""
FieldTrip buffer (V1) client in pure Python

Not complete yet - only simple GET_xxx calls implemented
(C) 2010 S. Klanke
"""

# We need socket, struct, and numpy 
import socket
import struct
import numpy

VERSION = 1
PUT_HDR = 0x101
PUT_DAT = 0x102
PUT_EVT = 0x103
PUT_OK  = 0x104
PUT_ERR = 0x105
GET_HDR = 0x201
GET_DAT = 0x202
GET_EVT = 0x203
GET_OK  = 0x204
GET_ERR = 0x205
FLUSH_HDR = 0x301
FLUSH_DAT = 0x302
FLUSH_EVT = 0x303
FLUSH_OK  = 0x304
FLUSH_ERR = 0x305
WAIT_POLL = 0x401

DATATYPE_CHAR = 0
DATATYPE_UINT8 = 1
DATATYPE_UINT16 = 2
DATATYPE_UINT32 = 3
DATATYPE_UINT64 = 4
DATATYPE_INT8 = 5
DATATYPE_INT16 = 6
DATATYPE_INT32 = 7
DATATYPE_INT64 = 8
DATATYPE_FLOAT32 = 9
DATATYPE_FLOAT64 = 10
DATATYPE_UNKNOWN = 0xFFFFFFFF

# List for converting FieldTrip datatypes to Numpy datatypes
numpyType = ['int8', 'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16', 'int32', 'int64', 'float32', 'float64']
# Corresponding word sizes
wordSize = [1,1,2,4,8,1,2,4,8,4,8]
# FieldTrip data type as indexed by numpy dtype.num
# this goes  0 => nothing, 1..4 => int8, uint8, int16, uint16, 7..10 => int32, uint32, int64, uint64  11..12 => float32, float64
dataType = [-1, 5, 1, 6, 2, -1, -1, 7, 3, 8, 4, 9, 10]

def serialize(A):
	"""Returns Fieldtrip data type and string representation of the given object, if possible."""
	if isinstance(A, str):
		return (0,A)
	
	if isinstance(A, numpy.ndarray):
		dt = A.dtype
		if not(dt.isnative) or dt.num<1 or dt.num>=len(dataType):
			return (DATATYPE_UNKNOWN, None)
			
		ft = dataType[dt.num]
		if ft == -1:
			return (DATATYPE_UNKNOWN, None)
			
		try:
			buf = A.data
		except:
		    buf = A.flatten().data
			
		return (ft, str(buf))
	
	if isinstance(A, int):
		return (DATATYPE_INT32, struct.pack('i', A))
	
	if isinstance(A, float):
		return (DATATYPE_FLOAT64, struct.pack('d', A))

	return (DATATYPE_UNKNOWN, None)


def recv_all(s, N):
	"""Receive N bytes from a socket 's' and return it as a string."""
	A = s.recv(N)
	while len(A)<N:
		A += s.recv(N-len(A))
	return A

	
def send_all(s, A):
	"""Send all bytes of the string 'A' out to socket 's'."""
	N = len(A);
	nw = s.send(A)
	while nw<N:
		nw += s.send(A[nw:])

	
class Chunk:
	def __init__(self):
		self.type = 0
		self.size = 0
		self.buf = ''

class Header:
	"""Class for storing header information in the FieldTrip buffer format"""
	def __init__(self):
		self.nChannels = 0
		self.nSamples = 0
		self.nEvents = 0
		self.fSample = 0.0
		self.dataType = 0
		self.chunks = []
		
	def __str__(self):
		return 'Channels.: %i\nSamples..: %i\nEvents...: %i\nSampFreq.: %f\nDataType.: %i\n'%(self.nChannels, self.nSamples, self.nEvents, self.fSample, self.dataType)
		
class Event:
	"""Class for storing events in the FieldTrip buffer format"""
	def __init__(self, S = None):
		if S is None:
			self.type = ''
			self.value = ''
			self.sample = 0
			self.offset = 0
			self.duration = 0
		else:
			self.deserialize(S)
	
	def __str__(self):
		return 'Type.....: %s\nValue....: %s\nSample...: %i\nOffset...: %i\nDuration.: %i\n'%(str(self.type),str(self.value), self.sample, self.offset, self.duration)
	
	def deserialize(self, buf):
		bufsize = len(buf)
		if bufsize < 32:
			return 0
		
		(type_type, type_numel, value_type, value_numel, sample, offset, duration, bsiz) = struct.unpack('IIIIIiiI', buf[0:32])
		
		self.sample = sample
		self.offset = offset
		self.duration = duration
		
		st = type_numel * wordSize[type_type]
		sv = value_numel * wordSize[value_type]
		
		if bsiz+32 > bufsize or st+sv > bsiz:
			raise IOError('Invalid event definition -- does not fit in given buffer')
			
		raw_type = buf[32:32+st]
		raw_value = buf[32+st:32+st+sv]
				
		if type_type == 0:
			self.type = raw_type
		else:
			self.type = numpy.ndarray((type_numel), dtype=numpyType[type_type], buffer=raw_type)

		if value_type == 0:
			self.value = raw_value
		else:
			self.value = numpy.ndarray((value_numel), dtype=numpyType[value_type], buffer=raw_value)

		return bsiz + 32
		
	def serialize(self):
		"""Returns the contents of this event as a string, ready to send over the network, 
		   or None in case of conversion problems.
		"""
		type_type, type_buf = serialize(self.type)
		if type_type == DATATYPE_UNKNOWN:
			return None
		type_size = len(type_buf)
		type_numel = type_size / wordSize[type_type]
			
		value_type, value_buf = serialize(self.value)
		if value_type == DATATYPE_UNKNOWN:
			return None
		value_size = len(value_buf)
		value_numel = value_size / wordSize[value_type]
		
		bufsize = type_size + value_size
		
		S = struct.pack('IIIIIiiI', type_type, type_numel, value_type, value_numel, int(self.sample), int(self.offset), int(self.duration), bufsize)
		return S + type_buf + value_buf 
		
class Client:
	"""Class for managing a client connection to a FieldTrip buffer."""
	def __init__(self):
		self.isConnected = False
		self.sock = []
	
	def connect(self, hostname, port=1972):
		"""connect(hostname [, port]) -- make a connection, default port is 1972."""
		self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		self.sock.connect((hostname, port))
		self.sock.setblocking(True)
		self.isConnected = True

	def disconnect(self):
		"""disconnect() -- close a connection."""
		if self.isConnected:
			self.sock.close()
			self.sock = []
			self.isConnected = False

	def getHeader(self):
		"""getHeader() -- grabs header information (no chunks yet) from the buffer an returns it as a Header object."""
		if not(self.isConnected):
			raise IOError('Not connected to FieldTrip buffer')
		
		request = struct.pack('HHI', VERSION, GET_HDR, 0)
		nWr = send_all(self.sock, request)
		resp_hdr = recv_all(self.sock, 8)
		
		(version, command, bufsize) = struct.unpack('HHI', resp_hdr)
		
		if command==GET_ERR:
			return None
		
		if version!=VERSION or command!=GET_OK:
			self.disconnect()
			raise IOError('Bad response from buffer server - disconnecting')
		
		if bufsize > 0:
			resp_buf = self.sock.recv(bufsize)
			if len(resp_buf) < bufsize:
				self.disconnect()
				raise IOError('Error during socket read operation')
				
		if bufsize < 24:
			self.disconnect()
			raise IOError('Invalid HEADER packet received (too few bytes) - disconnecting')
				
		(nchans, nsamp, nevt, fsamp, dtype, bfsiz) = struct.unpack('IIIfII', resp_buf[0:24])
		
		H = Header()
		H.nChannels = nchans
		H.nSamples = nsamp
		H.nEvents = nevt
		H.fSample = fsamp
		H.dataType = dtype
		
		# TODO: chunks
		
		return H
		

	def getData(self, index = None):
		"""getData([indices]) -- retrieve data samples and return them as a Numpy array, samples in rows(!).
			The 'indices' argument is optional, and if given, must be a tuple or list with inclusive, zero-based 
			start/end indices.
		"""
		if not(self.isConnected):
			raise IOError('Not connected to FieldTrip buffer')
			
		if index is None:
			request = struct.pack('HHI', VERSION, GET_DAT, 0)
			nWr = send_all(self.sock, request)
		else:
			indS = int(index[0])
			indE = int(index[1])
			request = struct.pack('HHIII', VERSION, GET_DAT, 8, indS, indE)
			nWr = send_all(self.sock, request)

		resp_hdr = recv_all(self.sock, 8)
		
		(version, command, bufsize) = struct.unpack('HHI', resp_hdr)
		
		if command == GET_ERR:
			return None
		
		if version!=VERSION or command!=GET_OK:
			self.disconnect()
			raise IOError('Bad response from buffer server - disconnecting')
		
		if bufsize > 0:
			resp_buf = recv_all(self.sock, bufsize)
				
		if bufsize < 16:
			self.disconnect()
			raise IOError('Invalid DATA packet received (too few bytes)')
				
		(nchans, nsamp, datype, bfsiz) = struct.unpack('IIII', resp_buf[0:16])
		
		if bfsiz < bufsize - 16 or datype >= len(numpyType):
			raise IOError('Invalid DATA packet received')
			
		raw = resp_buf[16:bfsiz+16]
		D = numpy.ndarray((nsamp, nchans), dtype=numpyType[datype], buffer=raw)
		
		return D
		
		
	def getEvents(self, index = None):
		"""getEvents([indices]) -- retrieve events and return them as a list of Event objects.
			The 'indices' argument is optional, and if given, must be a tuple or list with 
			inclusive, zero-based start/end indices. The 'type' and 'value' fields of the event
			will be converted to strings or Numpy arrays.
		"""
		if not(self.isConnected):
			raise IOError('Not connected to FieldTrip buffer')
			
		if index is None:
			request = struct.pack('HHI', VERSION, GET_EVT, 0)
			nWr = send_all(self.sock, request)
		else:
			indS = int(index[0])
			indE = int(index[1])
			request = struct.pack('HHIII', VERSION, GET_EVT, 8, indS, indE)
			nWr = send_all(self.sock, request)

		resp_hdr = recv_all(self.sock, 8)
		
		(version, command, bufsize) = struct.unpack('HHI', resp_hdr)
		
		if command == GET_ERR:
			return []
		
		if version!=VERSION or command!=GET_OK:
			self.disconnect()
			raise IOError('Bad response from buffer server - disconnecting')
		
		if bufsize > 0:
			resp_buf = recv_all(self.sock, bufsize)
			
		offset = 0
		E = []
		while 1:
			e = Event()
			nextOffset = e.deserialize(resp_buf[offset:])
			if nextOffset == 0:
				break
			E.append(e)
			offset = offset + nextOffset
		
		return E
		

# Just a small demo for testing purposes...
# This should be moved to a separate file at some point		
ftc = Client()
ftc.connect('localhost',1972)
print '\nConnected - trying to read header...'
print ftc.getHeader()	

print '\nTrying to read (all) data...'
D = ftc.getData()
print D.shape
print D

print '\nTrying to read (all) events...'
E = ftc.getEvents()
for e in E:
	print e

ftc.disconnect()