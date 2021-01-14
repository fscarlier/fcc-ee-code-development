import numpy as np
import re
import os


def save_seq(fname,fmask='./save_seq.mask',mexec='/machfs/swhite/madx/madx'):
    #save a madx sequence
    fin = fname.replace('.seq','_py.seq')
    fout = fname.replace('.seq','_for_at.seq')
    mfile = fmask.replace('.mask','.madx')
    fhandle = fhandle = open(fmask,'r')
    lines = fhandle.readlines()
    fhandle.close()
    fhandle = open(mfile,'w')
    for line in lines:
        line = line.replace('$FIN',fin)
        line = line.replace('$FOUT',fout)
        fhandle.write(line)
    fhandle.close()
    cmd = mexec+' < '+mfile
    os.system(cmd)
    os.system('rm -f ./'+fin)
    #step needed for at converter
    lines = open(fout, 'r').readlines()
    new_last_line = (lines[-1].rstrip() + '\nreturn;')
    lines[-1] = new_last_line
    open(fout, 'w').writelines(lines)
    
    

def reformat_seq(fname,nskip = 5):
    #some pre-format needed for the AT converter
    fout = fname.replace('.seq','_py.seq')
    fhandle = open(fname,'r')
    lines = fhandle.readlines()
    fhandle.close()
    fhandle = open(fout,'w')
    #read the variables
    n=nskip
    dvars={}
    line=lines[n]
    while line !='\n':
        fhandle.write(lines[n]) 
        sl = line.split('=')
        dvars[sl[0]]=sl[1]
        n=n+1
        line = lines[n]
    #read element definition
    n=n+1
    delems={}
    line=lines[n]
    while line !='\n': 
        sl = line.replace(';','').replace('\n','').split(' : ')
        delems[sl[0]]=sl[1]
        n=n+1
        line = lines[n]
    #read the sequence, replace '.' and write elements definition using scalars
    n=n+1
    dseq = []
    header=lines[n]
    for i in range(n+1,len(lines)-1):
        sl = re.split(': |,',lines[i])
        dseq.append(sl[0].replace('.','_')+' : '+sl[0].replace('.','_')+' , '+sl[2])
        typen = sl[0].split('.')[0]
        fhandle.write(sl[0].replace('.','_')+' : '+delems[typen]+';\n')
    #write the sequence
    fhandle.write(header)
    for l in dseq:
        fhandle.write(l)
    fhandle.write(lines[len(lines)-1])
    fhandle.close()

if __name__ == "__main__":
    fname = 'Optics_FCCee_h_217_nosol_3.seq'
    reformat_seq(fname,nskip = 5)
    save_seq(fname,fmask='./save_seq.mask',mexec='/machfs/swhite/madx/madx')
    
