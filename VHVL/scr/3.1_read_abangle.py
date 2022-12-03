import os

un = []
with open('../list/unable_Fv_apo') as l:
    lines = l.readlines()
    for i in range(len(lines)):
        line = lines[i]
        un.append(line[:6])

HL = -59.50
HC1 = 71.50
LC1 = 118.50
HC2 = 119.50
LC2 = 83.50
dc = 16.50

par = [HL, HC1, LC1, HC2, LC2, dc]
lpdb =[]
with open('../list/human_Fv_apo') as h:
    lines = h.readlines()
    for i in range(len(lines)):
        line = lines[i]
        pdb = line[:6]
        lpdb.append(pdb)

for pdb in lpdb:
    if pdb in un: continue
    if not os.path.exists('./3.0_run/%s/ABangleData/UserAngles.dat'%pdb):
        print('ERROR %s'%pdb)
    with open('./3.0_run/%s/ABangleData/UserAngles.dat'%pdb) as d:
        lines = d.readlines()
        HL = float(lines[1].split()[1])
        HC1 = float(lines[1].split()[2])
        LC1 = float(lines[1].split()[3])
        HC2 = float(lines[1].split()[4])
        LC2 = float(lines[1].split()[5])
        dc = float(lines[1].split()[6])
        check = [HL, HC1, LC1, HC2, LC2, dc]
        if check != par:
            print('wrong orientation %s'%pdb)
