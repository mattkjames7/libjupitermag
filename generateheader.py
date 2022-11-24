import numpy as np
import os

headerfiles = [ 'src/libjupitermag.h',
                'src/trace.h',
                'src/model.h',
                'src/traceclosestpos.h',
                'lib/libcon2020/include/con2020.h',
                'lib/libspline/include/spline.h',
                'lib/libinternalfield/include/internalfield.h']

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

def HeaderDef(fname):
    '''
    Return the string used as header define


    '''
    bn = os.path.basename(fname)
    n,e = bn.split('.')
    out = '__{:s}_{:s}__'.format(n.upper(),e.upper())
    return out

def ExtractHeaderLines(fname):
    '''
    Exctract the lines we want in our 
    new header file.

    Inputs
    ======
    fname : str
        Name of the header file

    Returns
    =======
    include : np.array
        List of include lines
    defines : np.array
        List of defines
    externs : np.array
        List of extern prototypes
    others : np.array
        List of everything else to include
    

    '''


    lines = ReadFile(fname)

    #remove header define etc
    remstr = [HeaderDef(fname),'#ifndef','#endif']
    keep = np.ones(lines.size,dtype='bool')
    for i,l in enumerate(lines):
        for r in remstr:
            if r in l:
                keep[i] = False
                break
    lines = lines[keep]


    #extract include lines first
    isinclude = np.zeros(lines.size,dtype='bool')

    for i,l in enumerate(lines):
        if l.startswith('#include'):
            isinclude[i] = True

    include = lines[isinclude]
    lines = lines[isinclude == False]

    #extract defines
    isdefine = np.zeros(lines.size,dtype='bool')
    keep = []
    for i,l in enumerate(lines):
        if l.startswith('#define'):
            isdefine[i] = True
        else:
            keep.append(i)
    keep = np.array(keep)
    defines = lines[isdefine]
    lines = lines[keep]

    #locate extern C 
    br_cnt = 0
    externs = []
    keep = []
    for i,l in enumerate(lines):
        if 'extern "C"' in l:
            br_cnt = 1
            externs.append(l)
            continue
        if br_cnt > 0:
            br_cnt += l.count('{')
            br_cnt -= l.count('}')
            externs.append(l)
        else:
            keep.append(i)
    keep = np.array(keep)
    lines = lines[keep]
    externs = np.array(externs)[1:-1]

    #remove extra empty lines
    isempty = np.zeros(lines.size,dtype='bool')
    for i,l in enumerate(lines):
        if len(l.strip()) == 0:
            isempty[i] = True
    remove = np.append(isempty[:-1] & isempty[1:],False)
    keep = np.where(remove == False)[0]
    others = lines[keep]

    return include,defines,externs,others
    
def Version():

    with open('VERSION','r') as f:
        line = f.readline()
    
    line = line.strip()
    mj,mn,pa = line.split('.') 

    out = [ '#define LIBJUPITERMAG_VERSION_MAJOR '+mj+'\n',
            '#define LIBJUPITERMAG_VERSION_MINOR '+mn+'\n',
            '#define LIBJUPITERMAG_VERSION_PATCH '+pa+'\n',]
    return out

def ListFiles(start,ReturnNames=False):
	'''
	Should list the files that exist within a folder.
	'''
	
	FileOut = []
	NameOut = []
	for root,dirs,files in os.walk(start,topdown=False,followlinks=True):
		for name in files:
			FileOut.append(root+'/'+name)
			NameOut.append(name)
	
	FileOut = np.array(FileOut)
	NameOut = np.array(NameOut)
	
	if ReturnNames:
		return FileOut,NameOut
	else:
		return FileOut



def CombineHeaders():



    # read everything in
    include = np.array([])
    defines = np.array([])
    externs = np.array([])
    others = np.array([])
    for hf in headerfiles:
        i,d,e,o = ExtractHeaderLines(hf)
        include = np.append(include,i)
        defines = np.append(defines,d)
        externs = np.append(externs,e)
        others = np.append(others,o)

    #remove duplicates
    include = np.unique(include)
    defines = np.unique(defines)

    #list the names of headers
    files = ListFiles('.')
    hbn = []
    for f in files:
        _,ext = os.path.splitext(f)
        if ext in ['.h','.hpp']:
            hbn.append(os.path.basename(f))

    #remove unwanted bits from the includes
    keep = np.ones(include.size,dtype='bool')
    for i,inc in enumerate(include):
        for h in hbn:
            if h in inc:
                keep[i] = False
                break
    include = include[keep]

    #frame the externs
    externs = np.append('extern "C" {\n',externs)
    externs = np.append(externs,'}\n')

    #get the current version
    ver = Version()

    #combine everything
    out = np.array(['#ifndef __LIBJUPITERMAG_H__\n',
                    '#define __LIBJUPITERMAG_H__\n'])
    out = np.append(out,include)
    out = np.append(out,'\n')
    out = np.append(out,ver)
    out = np.append(out,'\n')
    out = np.append(out,defines)
    out = np.append(out,'\n')
    out = np.append(out,externs)
    out = np.append(out,'\n')
    out = np.append(out,others)
    out = np.append(out,'#endif\n')


    #create the output file
    with open('include/jupitermag.h','w') as f:
        f.writelines(out)
        print('Saved header file')

if __name__ == '__main__':

    CombineHeaders()