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
GET_HDR = 513
GET_DAT = 514
GET_EVT = 515
GET_OK  = 516
GET_ERR = 517

numpyType = ['int8', 'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16', 'int32', 'int64', 'float32', 'float64']
wordSize = [1,1,2,4,8,1,2,4,8,4,8]

def recv_all(s, N):
	A = s.recv(N)
	while len(A)<N:
		A += s.recv(N-len(A))
	return A

	
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
	def __init__(self):
		self.type = ''
		self.value = ''
		self.sample = 0
		self.offset = 0
		self.duration = 0
	
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
		nWr = self.sock.send(request)
		resp_hdr = self.sock.recv(8)
		nRd = len(resp_hdr)
		
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
			nWr = self.sock.send(request)
		else:
			indS = int(index[0])
			indE = int(index[1])
			request = struct.pack('HHIII', VERSION, GET_DAT, 0, indS, indE)
			nWr = self.sock.send(request)

		resp_hdr = self.sock.recv(8)
		nRd = len(resp_hdr)
		
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
			nWr = self.sock.send(request)
		else:
			indS = int(index[0])
			indE = int(index[1])
			request = struct.pack('HHIII', VERSION, GET_EVT, 8, indS, indE)
			nWr = self.sock.send(request)

		resp_hdr = self.sock.recv(8)
		nRd = len(resp_hdr)
		
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
ftc.connect('mentat069',1972)
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