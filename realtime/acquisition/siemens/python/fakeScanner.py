import time
import dircache
import shutil
import random
import os, errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


print '\nPress CTRL-C to stop\n'

random.seed()

srcdir = 'D:/fmri_data/46p-4064/'
#srcdir = 'D:/fmri_data/49p-4064/'
#srcdir = 'D:/fmri_data/52p-4064/'
#srcdir = 'D:/fmri_data/47p-4064/'
destdir = 'D:/watch/image/%ip-%i'%(random.randint(10,100), random.randint(1000,10000))
mkdir_p(destdir)

FN = dircache.listdir(srcdir)

PixFN = []
for fn in FN:
	if '.PixelData' in fn:
		PixFN.append(fn)

numSrc = len(PixFN)
numDest = 25

idxSrc  = 0
idxDest = 1

logFile = file(destdir + '-logging.txt','w')

shutil.copy(srcdir + 'mrprot.txt', 'D:/watch/image/mrprot.txt')

echoTimes = [1.5, 0.01, 0.03, 0.25]
#echoTimes = [1];
curEcho = 0;

while 1:
	time.sleep(echoTimes[curEcho])
	curEcho = curEcho + 1
	if curEcho == len(echoTimes):
		curEcho = 0
	
	T = time.time()
	shutil.copy(srcdir + PixFN[idxSrc], '%s/Scan%04i.PixelData'%(destdir,idxDest))
	info = '%16.4f  %04i  %s\n'%(T,idxDest,PixFN[idxSrc])
	print info
	
	#shutil.copy(srcdir + PixFN[idxSrc], '%s/Scan%04ip.PixelData'%(destdir,idxDest))
	#info = '%16.4f  %04ip  %s\n'%(T,idxDest,PixFN[idxSrc])
	#print info
	
	logFile.write(info)
	logFile.flush()
	
	idxSrc = idxSrc+1
	idxDest = idxDest+1
	if idxSrc>=numSrc:
		idxSrc=0
	if idxDest>numDest:
		idxDest=1
	
