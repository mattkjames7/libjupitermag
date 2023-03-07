import numpy as np
import os

headerfiles = [ 'src/coordconv.h',
				'src/footprint.h',
				'lib/libspline/include/spline.h',
				'src/interptraceclosestpos.h',
				'lib/libcon2020/include/con2020.h',
				'lib/libinternalfield/include/internalfield.h',
				'src/model.h',
				'src/trace.h',
				'src/libjupitermag.h']

def ReadFile(fname):
	'''
	Simply read in the lines from an ASCII file

	Input
	=====
	fname : str
		Name of the file to read

	Returns
	=======
	lines : np.array
		Array of strings - one for each line of the file.

	'''

	with open(fname,'r') as f:
		lines = f.readlines()

	return np.array(lines)

def _RemoveDirectives(lines):
	'''
	Remove compiler directives and includes
	'''
	lines = np.array(lines)
	nl = lines.size

	use = np.ones(nl,dtype='bool')
	for i in range(0,nl):
		if lines[i].strip().startswith('#'):
			use[i] = False

	use = np.where(use)

	return lines[use]

def _SplitHeaderDefs(lines):
	'''
	split code into C and C++ code

	'''
	lines = np.array(lines)
	ltype = np.zeros(lines.size,dtype='int')
	isC = False
	for i in range(0,lines.size):

		if isC and lines[i].strip() == '}':
			isC = False
			ltype[i] = 0
		elif isC:
			ltype[i] = 1
		else:
			ltype[i] = 2
		if 'extern "C"' in lines[i].strip():
			isC = True
			ltype[i] = 0

	usec = np.where(ltype == 1)[0]
	usecc = np.where(ltype == 2)[0]

	c = lines[usec]
	cc = lines[usecc]

	return c,cc

def _ReadHeader(fname):
	'''
	This will read a header file in and  remove
	any compiler directives and split into C/C++
	
	'''

	lines = ReadFile(fname)


	hasextC = False
	for l in lines:
		if 'extern "C"' in l:
			hasextC = True
			break

	code = _RemoveDirectives(lines)

	if hasextC:
		c,cc = _SplitHeaderDefs(code)
	else:
		c = []
		cc = code

	return c,cc

def GenerateHeader():


	#get the template file first
	top = ReadFile('include/jupitermag.template')
	
	#add version 
	ver = Version()

	top = np.append(top,np.array(ver))

	#a few more lines
	extc = np.array(['#ifdef __cplusplus\n','extern "C" {\n','#endif\n'])
	top = np.append(top,extc)



	c = []
	cc = []


	for f in headerfiles:
		fc,fcc = _ReadHeader(f)

		if len(fc) > 0:
			c = c + fc.tolist()
		if len(fcc) > 0:
			cc = cc + fcc.tolist()
	
	cc = np.array(cc)
	keep = np.ones(cc.size,dtype='bool')
	for i in range(0,cc.size):
		if 'typedef void (*modelFieldPtr)(double,double,double,double*,double*,double*);' in cc[i]:
			keep[i] = False
	use = np.where(keep)[0]
	cc = cc[use]

	cpp = '#ifdef __cplusplus\n}\n'

	bottom = '''
#endif
#endif
	'''
	print('Saving header file: include/jupitermag.h')
	f = open('include/jupitermag.h','w')
	f.writelines(top)
	f.writelines(c)
	f.write(cpp)
	f.writelines(cc)
	f.write(bottom)
	f.close()

	
	
def Version():

	with open('VERSION','r') as f:
		line = f.readline()
	
	line = line.strip()
	mj,mn,pa = line.split('.') 

	out = [ '#define LIBJUPITERMAG_VERSION_MAJOR '+mj+'\n',
			'#define LIBJUPITERMAG_VERSION_MINOR '+mn+'\n',
			'#define LIBJUPITERMAG_VERSION_PATCH '+pa+'\n',]
	return out


if __name__ == '__main__':

	GenerateHeader()