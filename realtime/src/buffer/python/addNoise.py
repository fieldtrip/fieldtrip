"""
Demo for FieldTrip buffer Python interface, 
simulating an asynchronous pipeline element.

This script will read data and events from one
FT buffer, add noise, and write the results to another
buffer. By default, the buffer addresses are localhost:1972
and localhost:1973, but the port numbers can be changed by
command line arguments (or in the script).

Also, this script will poll for events, and when finding an
event with type = 'noise' and a value that can be turned into a 
floating point number, that value will be used as the new noise
level.

(C) 2010 S. Klanke
"""


import FieldTrip
import sys
import time
import numpy

hostname = 'localhost'
portIn = 1972
portOut = 1973
	
if len(sys.argv)>1:
	hostname = sys.argv[1]
if len(sys.argv)>2:
	try:
		portIn = int(sys.argv[2])
	except:
		print 'Error: second argument (%s) must be a valid (=integer) port number'%sys.argv[2]
		sys.exit(1)
if len(sys.argv)>3:
	try:
		portOut = int(sys.argv[3])
	except:
		print 'Error: second argument (%s) must be a valid (=integer) port number'%sys.argv[2]
		sys.exit(1)		
	
ftIn  = FieldTrip.Client()		
ftOut = FieldTrip.Client()		
		
print 'Trying to connect to buffer on %s:%i ...'%(hostname,portIn)
ftIn.connect(hostname, portIn)
print 'Trying to connect to buffer on %s:%i ...'%(hostname,portOut)
ftOut.connect(hostname, portOut)

noiseLevel = 0.1

while 1:
	H = ftIn.getHeader()
	if H is None:
		print 'Header could not be read.'
		time.sleep(0.5)
		continue
		
	ftOut.putHeader(H.nChannels, H.fSample, H.dataType) # H.labels
	
	numSmp = H.nSamples
	numEvt = H.nEvents
	
	
	while 1:
		newSmp, newEvt = ftIn.wait(numSmp, numEvt, 500)
		
		if newSmp < numSmp:
			print 'Number of samples decreased - re-reading header.'
			break
		
		if newEvt > numEvt:
			E = ftIn.getEvents([numEvt, newEvt-1])
			for e in E:
				print e
				if e.type == 'noise':
					try:
						noiseLevel = float(e.value)
					except:
						pass
			numEvt = newEvt
						
		if newSmp == numSmp:
			continue
		
		D = ftIn.getData([numSmp, newSmp-1])
		
		noise = numpy.float32(numpy.random.randn(*D.shape))
		
		Dn = D + noiseLevel*noise
		
		ftOut.putData(Dn)
		
		print 'Send out samples %i - %i'%(numSmp,newSmp-1)
		
		numSmp = newSmp

