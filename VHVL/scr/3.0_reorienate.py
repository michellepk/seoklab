import os, math
import numpy as np
import tempfile, subprocess

db_path='/home/minjae_pk/DB/SAbDab/human_Fv_apo/0.2_Fv/2.0_HL/'
data_path='/home/minjae_pk/GalaxyPipe/lib/abangle/data/'

# Get coreset positions
coresetL = [l.strip()[1:] for l in open(os.path.join(data_path, "Lcoresetfw.txt")).readlines()]
coresetH = [l.strip()[1:] for l in open(os.path.join(data_path, "Hcoresetfw.txt")).readlines()]

# Read in the plane vectors precalculated on the consensus structure
Lpos = map( lambda x: map(float, x), map(str.split, open(os.path.join(data_path, "pcL.txt")).readlines()))
Hpos = map( lambda x: map(float, x), map(str.split, open(os.path.join(data_path, "pcH.txt")).readlines())) 

#########################
# Calculation functions #
#########################

def create_coreset(fname):
    try:
        # Parse the input file
        Hfd, Hf = tempfile.mkstemp('.pdb','H')
        Lfd, Lf = tempfile.mkstemp('.pdb','L')

        Htmp = os.fdopen(Hfd, 'w')
        Ltmp = os.fdopen(Lfd, 'w')
        fin = open(fname, 'r').readlines()
        for line in fin:
            l = line.split()
            if l[0] != 'ATOM': continue
            elif l[4] == 'L' and l[5] in coresetL:
                Ltmp.write(line)
            elif l[4] == 'H' and l[5] in coresetH:
                Htmp.write(line)
        Htmp.close()
        Ltmp.close()

    except Exception, e:
        os.remove(Hf)
        os.remove(Lf)
        raise Exception(str(e)+"\n")
    return Hf, Lf

def mapvectors(fname, PAPS_def = False):
    # Get transformation matrices by alining the core of the domains
    Hf, Lf = create_coreset(fname)

    uL = align(os.path.join(data_path,"consensus_L.pdb"), Lf)
    uH = align(os.path.join(data_path,"consensus_H.pdb"), Hf)
    os.remove(Hf)
    os.remove(Lf)

    if PAPS_def:
        # centroid of interface residues
        cH = Hpos[2]
        cL = Lpos[2]
    else:
        # The minimally varying centroid vector is at as calculated:
        cH = [ -10*0.5*Hpos[0][i] + 1*0.5*Hpos[1][i] + Hpos[2][i] for i in range(3) ]
        cL = [ 6*0.5*Lpos[0][i] - 2*0.5*Lpos[1][i] + Lpos[2][i] for i in range(3) ]

    # Define the plane vectors from the centroid point
    # On VL domain
    L1 = [ cL[i] + Lpos[0][i] for i in range(3) ]
    L2 = [ cL[i] + Lpos[1][i] for i in range(3) ]

    # On VH domain
    H1 = [ cH[i] + Hpos[0][i] for i in range(3) ]
    H2 = [ cH[i] + Hpos[1][i] for i in range(3) ]

    # Do the transformation onto the
    Lpoints = map( lambda x: transform(x, uL), (cL, L1, L2) )
    Hpoints = map( lambda x: transform(x, uH), (cH, H1, H2) )

    return Lpoints, Hpoints

def align(file1, file2):
    """ Aligns file1 to file2 using tmalign and returns the trasnformation matrix """
    # Temp file for the matrix for latest versions of TMalign
    mtmpfd, mtmp = tempfile.mkstemp('.txt', 'matrix')
    os.close(mtmpfd)
    # Align file1 to file2 using TMalign

    try:
        subpr = subprocess.Popen(['TMalign',file1,file2,'-m',mtmp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        TMresult = subpr.communicate()
    
    except OSError:
        # If is not found, point to webpage for installation
	raise Exception('Cannot execute TMalign. Please install and ensure it is in your path.\nTMalign can be downloaded from:\n'\
	'http://zhanglab.ccmb.med.umich.edu/TM-align/\n'\
	'Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9\n' ) 

    # Parse the output of TMalign. Some versions don't output the matrix
    result = TMresult[0].split("\n")
    attempt = 0
    while 1:
        try:
            i = 0
            while 1:
		if result[i].upper().startswith(' -------- ROTATION MATRIX'):
		    # Grab transformation matrix		
		    u = []
		    u.append( map (float, result[i+2].split()[1:] ) )
		    u.append( map (float, result[i+3].split()[1:] ) )
		    u.append( map (float, result[i+4].split()[1:] ) )
		    break
		else:
		    i+=1
	    break
	except IndexError:
	    try:
		if not attempt:
		    ftmp = open( mtmp )
		    result = ftmp.readlines()
		    ftmp.close()
		    attempt=1
	    except IOError:
		break
	    
    if os.path.exists( mtmp ):
        os.remove( mtmp )

    # Return the transformation matrix
    try:
        return u
    except NameError:
	raise Exception("TMalign alignment file not in an expected format, check output gives rotation matrix (or with -m option )\n")

def transform(coords, u):
    """ Transforms coords by a matrix u. u is found using tmalign """
    # Ensure coordinates are of type float
    coords = map( float, coords )

    # Do transformation
    X = u[0][0]+u[0][1]*coords[0] +u[0][2]*coords[1] +u[0][3]*coords[2] 
    Y = u[1][0]+u[1][1]*coords[0] +u[1][2]*coords[1] +u[1][3]*coords[2] 
    Z = u[2][0]+u[2][1]*coords[0] +u[2][2]*coords[1] +u[2][3]*coords[2] 

    # Return transformed coordinates
    return [X,Y,Z]
   
def normalise(vec):
    mag = (sum(map(lambda x: x**2, vec)))**0.5
    return map(lambda x: x/mag, vec)

def length(crd):
    d = 0
    for i in range(3):
        d+=crd[i]**2
    dd = d**0.5
    return dd

# input degree to radian
def degree_radian(degree):
    rad = degree*(math.pi/180.0)
    return rad

un = [l.strip()[1:] for l in open('../list/unable_Fv_apo').readlines()]
lapo = []
with open('../list/human_Fv_apo') as apo:
    lines=apo.readlines()
    for i in range(len(lines)):
        line=lines[i]
        if line[:6] in un: continue
        lapo.append(line[:6])

def original_matrix(fname):
    Lpoints, Hpoints = mapvectors(fname)
    pcL = Lpoints[0]
    pL1 = Lpoints[1]
    pL2 = Lpoints[2]

    L1 = np.subtract(pL1, pcL)
    L2 = np.subtract(pL2, pcL)
    cL = np.cross(L1, L2)

    A = np.array([[cL[0], L1[0], L2[0]],
                    [cL[1], L1[1], L2[1]],
                    [cL[2], L1[2], L2[2]]])

    return A

def new_matrix(fname, dc, HC1, LC1, HC2, LC2, HL):
    # all degrees given as radian

    A = original_matrix(fname)
    Lpoints, Hpoints = mapvectors(fname)

    pcL = Lpoints[0]
    pL1 = Lpoints[1]
    pL2 = Lpoints[2]

    pcH = Hpoints[0]
    pH1 = Hpoints[1]
    pH2 = Hpoints[2]
   
    L1 = np.subtract(pL1, pcL)
    L2 = np.subtract(pL2, pcL)

    H1 = np.subtract(pH1, pcH)
    H2 = np.subtract(pH2, pcH)
    
    nL1 = normalise(L1)
    nL2 = normalise(L2)
    nH1 = normalise(H1)
    nH2 = normalise(H2)

    L12 = np.dot(nL1, nL2)
    L12 = math.acos(L12)

    # derive new pcL
    np.set_printoptions(precision=16)
    N = np.cross(nH1, nH2)
    d1 = -(length(nH1)*dc*np.cos(HC1) + np.dot(nH1, pcH))
    d2 = -(length(nH2)*dc*np.cos(HC2) + np.dot(nH2, pcH))
    c = 0 
    b = (d2*nH1[0]-d1*nH2[0])/(nH1[1]*nH2[0]-nH1[0]*nH2[1])
    a = -(nH1[1]*b+nH1[2]*c+d1)/nH1[0]
    P = np.array([length(N)**2, 2*(N[0]*(a-pcH[0])+N[1]*(b-pcH[1])+N[2]*(c-pcH[2])), (a-pcH[0])**2+(b-pcH[1])**2+(c-pcH[2])**2-dc**2])
    
    for i in range(2) :
        x = N[0]*np.roots(P)[i] + a
        y = N[1]*np.roots(P)[i] + b
        z = N[2]*np.roots(P)[i] + c
        cL = np.array([x,y,z])
        oC = np.subtract(pcL,pcH)
        nC = np.subtract(cL, pcH)
        if np.dot(oC, nC) > 0:
            pcL = cL

    C = np.subtract(pcH, pcL)
    Cminus = np.subtract(pcL, pcH)
    nC = normalise(C)
    nCminus = normalise(Cminus)

    # derive L1
    n_x = np.cross(nH1, nCminus)
    n_y = np.cross(nCminus, n_x)
    nnx = normalise(n_x)
    nny = normalise(n_y)
    tmpC_ = normalise([np.dot(nCminus, nCminus), np.dot(nCminus, nnx), np.dot(nCminus, nny)])
    tmpCm = normalise([np.dot(nC, nCminus), np.dot(nC, nnx), np.dot(nC, nny)])
    tmpH_ = normalise([0, np.dot(nH1, n_x), np.dot(nH1, n_y)])

    if HC1 < 90.0:
        alpha = 90.0 + HC1
    tmpH1 = normalise([np.dot(nH1, nCminus), np.dot(nH1, nnx), np.dot(nH1, nny)])
    
    tmpL_ = normalise([0, -np.sin(HL), np.cos(HL)])
    
    HL = math.acos(np.dot(tmpL_, tmpH_))
    HL = HL*(180.0/math.pi)

    a = np.cos(LC1)
    b = tmpL_[1]*(1-a**2)**0.5
    c = tmpL_[2]*(1-a**2)**0.5
    tmpL1 = normalise([-a, b, c])
    theta = math.acos(np.dot(tmpC_, tmpL1))
    theta = theta*(180.0/math.pi)

    N = np.array([[nCminus[0],nnx[0],nny[0]],
                    [nCminus[1],nnx[1],nny[1]],
                    [nCminus[2],nnx[2],nny[2]]])
    
    r = np.matmul(N, tmpL1)
    pL1 = np.add(pcL, r)
    L1 = np.subtract(pL1, pcL)
    nL1 = normalise(L1)
    
    LC1 = math.acos(np.dot(nC, nL1))
    LC1 = LC1*(180.0/math.pi)
    
    # derive L2
    N = np.cross(nC, nL1)
    d1 = -(np.dot(nC, pcL)+length(L2)*np.cos(LC2))
    d2 = -(np.dot(nL1, pcL)+length(L2)*np.cos(L12))
    
    c = 0
    b = (nC[0]*d2-nL1[0]*d1)/(nL1[0]*nC[1]-nL1[1]*nC[0])
    a = -(d1+nC[1]*b+nC[2]*c)/nC[0]
    ## problem!!##
    P = np.array([length(N)**2, 2*(N[0]*(a-pcL[0])+N[1]*(b-pcL[1])+N[2]*(c-pcL[2])), (a-pcL[0])**2+(b-pcL[1])**2+(c-pcL[2])**2-length(L2)**2])
    D = P[1]**2-4*P[0]*P[2]
    
    for i in range(2):
        x = N[0]*np.roots(P)[i] + a
        y = N[1]*np.roots(P)[i] + b
        z = N[2]*np.roots(P)[i] + c

        tL2 = np.array([x, y, z])
        nL2 = np.subtract(tL2, pcL)
        if np.dot(nL2, L2) > 0:
            pL2 = tL2

    L2 = np.subtract(pL2, pcL)
    crL = np.cross(L1, L2)

    B = np.array([[crL[0], L1[0], L2[0]],
                    [crL[1], L1[1], L2[1]],
                    [crL[2], L1[2], L2[2]]])
    
    n_x = np.cross(nL1, nC)
    n_y = np.cross(nC, n_x)

    tmpH_ = normalise([0, np.dot(nH1, n_x), np.dot(nH1, n_y)])
    tmpL_ = normalise([0, np.dot(nL1, n_x), np.dot(nL1, n_y)])

    HL = math.acos(np.dot(tmpL_, tmpH_))
    HL = HL*(180.0/math.pi)
    if np.dot(np.cross(tmpL_,tmpH_), [1,0,0]) < 0:
        HL = -HL
    
    return A, B, pcL

def transformed_pdb(fname, save_path, dc, HC1, LC1, HC2, LC2, HL):
    Lpoints, Hpoints = mapvectors(fname)
    A, B, pcL = new_matrix(fname, dc, HC1, LC1, HC2, LC2, HL)
    ocL = Lpoints[0] 
    pdb = fname.split('/')[-1]
    f = open('%s%s'%(save_path,pdb),'w')
    with open(fname) as o:
        lines = o.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if line[20:22].strip()=='H':
                f.write(line)
            if line[20:22].strip()=='L':
                x = float(line[27:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                crd = np.array([x, y, z])
                PP = np.subtract(crd, ocL)
                inv = np.linalg.inv(A)
                anw = np.matmul(inv, PP)
                rep =  np.add(pcL, np.matmul(B, anw))
                f.write('%s%11.3f%8.3f%8.3f%s\n'%(line[:27], rep[0], rep[1], rep[2], line[54:78]))
    f.close()
                

if not os.path.exists('../3.0_reorienate/'):
    os.mkdir('../3.0_reorienate/')

save_path = '../3.0_reorienate/'
dc = 16.5
HC1 = degree_radian(71.5)
LC1 = degree_radian(118.5)
HC2 = degree_radian(119.5)
LC2 = degree_radian(83.5)
HL = degree_radian(-59.5)

for pdb in lapo:
    pdb_fn = '%s.pdb'%pdb
    fn_path = '%s%s'%(db_path, pdb_fn)
    transformed_pdb(fn_path, save_path, dc, HC1, LC1, HC2, LC2, HL)
