"""
Demo for FieldTrip buffer Python interface, 
showing off the putEvents function.

This is intended to be used in interactive mode.
Type setNoise(X) with some number X to send an
event as described in addNoise.py.

(C) 2010 S. Klanke
"""

import FieldTrip
import sys
import numpy

hostname = 'localhost'
port = 1972
	
if len(sys.argv)>1:
	hostname = sys.argv[1]
if len(sys.argv)>2:
	try:
		port = int(sys.argv[2])
	except:
		print 'Error: second argument (%s) must be a valid (=integer) port number'%sys.argv[2]
		sys.exit(1)
	
ftc  = FieldTrip.Client()		
		
print 'Trying to connect to buffer on %s:%i ...'%(hostname,port)
ftc.connect(hostname, port)

def setNoise(val):
	E = FieldTrip.Event()
	ns, ne = ftc.poll()
	E.type = 'noise'
	E.sample = ns
	E.value = val
	ftc.putEvents(E)

