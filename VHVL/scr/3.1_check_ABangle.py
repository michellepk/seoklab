import os, time

un=[]
with open('../list/unable_Fv_apo') as l:
    lines=l.readlines()
    for i in range(len(lines)):
        line=lines[i]
        un.append(line[:6])

if not os.path.exists('./3.0_run'):
    os.mkdir('./3.0_run')

with open('../list/human_Fv_apo') as h:
    lines=h.readlines()
    for i in range(len(lines)):
        line=lines[i]
        pdb= line[:6]
        if line[:6] in un: continue
        if not os.path.exists('./3.0_run/%s'%pdb):
            os.mkdir('./3.0_run/%s'%pdb)
        os.system('ABangle -i ../3.0_reorienate/%s.pdb -usernumbered -change_data_path ./3.0_run/%s -store y'%(pdb,pdb))
