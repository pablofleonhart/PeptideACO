import os
import sys
import shutil
import numpy as np
import random
from time import time
import multiprocessing
import datetime
import math
from scipy.stats import norm
import matplotlib.pyplot as plt
from collections import defaultdict
from operator import itemgetter
import copy
import csv
from pdbReader import PDBReader
from pdbAligner import PDBAligner
from prediction import Prediction
from aminoPhiPsi import AminoPhiPsi
import rmsd

class ACOR:
    pdbPattern = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s} {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
    NHC_ATOMS = ("N", "H", "1H", "H1", "2H", "H2", "3H", "H3", "CA")
    generations = []
    values = []
    mod = []

    '''path = os.getcwd()
    # get folder name from input
    case = str(sys.argv[1])
    if os.path.exists(path+'/'+case):
        raw_input('Folder exists! Press enter to delete it and continue')
        shutil.rmtree(path+'/'+case)
    os.mkdir(path+'/'+case)
    os.chdir(path+'/'+case)'''

    def rotate_to( self, ang, atoms, aminoAcids, posAtoms ):
        #print "ang", ang
        angles = [6.28319, 3.14159, -3.14159, 3.14159, -3.14159, -3.14159, 3.14159, 3.14159, -3.14159,  6.28319]
        n_aa = len( aminoAcids )
        for i in xrange(n_aa):
            if i + min( aminoAcids) <= max(aminoAcids):
                #ROTATE PHI
                #print atoms, aminoAcids
                n_i = zip(atoms, aminoAcids).index((" N  ", i + min(aminoAcids)))   
                ca_i = zip(atoms, aminoAcids).index((" CA ", i + min(aminoAcids)))
                current_angles = angles
                #print current_angles
                dphi = math.atan2(math.sin(ang[2*i] - current_angles[2*i]), math.cos(ang[2*i] - current_angles[2*i]))
                #print "dphi", degrees( dphi )
                n_pos = posAtoms[n_i]
                ca_pos = posAtoms[ca_i]                
                ia = 0
                for atom in zip(atoms, aminoAcids):
                    if (i > 0) and (atom[1] > i + min(aminoAcids) or (atom[1] == i + min(aminoAcids) and (atom[0].strip() not in self.NHC_ATOMS))): 
                        posAtoms[ia] = rotate_atom_around_bond(dphi, posAtoms[ia], n_pos, ca_pos)
                        #print(atom[0], atom[1])   
                    ia += 1        
                #ROTATE PSI    
                c_i  = zip(atoms, aminoAcids).index((" C  ",  i + min(aminoAcids)))  
                ca_i = zip(atoms, aminoAcids).index((" CA ", i + min(aminoAcids)))
                current_angles = angles
                #print current_angles
                dpsi = math.atan2(math.sin(ang[2*i+1] - current_angles[2*i+1]), math.cos(ang[2*i+1] - current_angles[2*i+1]))              
                c_pos = posAtoms[c_i] 
                ca_pos = posAtoms[ca_i]
                ia = 0
                for atom in zip(atoms, aminoAcids):
                    if (i+min(aminoAcids) < max(aminoAcids)) and (atom[1] > i+min(aminoAcids) or (atom[1] == i+min(aminoAcids) and (atom[0]==" O  "))): 
                        posAtoms[ia] = rotate_atom_around_bond(dpsi, posAtoms[ia], ca_pos, c_pos)
                        #print(atom[0], atom[1])          
                ia += 1
                
    def normalize( self, v ):
        norm = np.linalg.norm( v )
        if norm == 0: 
            return v
        return v/norm  

    def rotate_atom_around_bond( self, theta, atom_pos, bond_start, bond_end ):
        #https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        v = np.array( atom_pos ) - np.array( bond_start )
        k = np.array( bond_end ) - np.array( bond_start )
        k = normalize( k )
        rot_pos = v * np.cos( theta ) + ( np.cross( k, v ) ) * np.sin( theta ) + k * ( np.dot( k, v ) ) * ( 1.0 - np.cos( theta ) )
        return list( rot_pos + np.array( bond_start ) )

    def calcKabschRMSD( self, exp, mod ):
        P = np.array( exp )
        Q = np.array( mod )
        #print rmsd.kabsch_rmsd( P, Q )
        P -= rmsd.centroid( P )
        Q -= rmsd.centroid( Q )
        result = rmsd.kabsch_rmsd( P, Q )
        #print "{:15s} {:6.2f}".format( "Kabsch RMSD:", result )
        return result

    def evaluator( self, x, refPosAtoms, modPosAtoms ):
        '''rotation = []
        for i in xrange( len( x )+2 ):
            if i == 0 or i == len( x )+1:
                rotation[i] = 0.0
            else:
                rotation[i] = x[i-1]*2*math.pi'''
        rotation = [ (2*math.pi*i)-math.pi for i in x ]
        rotation.append( 0.0 )

        rt = []
        rt.append( 0.0 )
        for i in xrange( len( rotation ) ):
            rt.append( rotation[i] )

        #print rt
        mod = copy.deepcopy( modPosAtoms )

        rotate_to( rt, prediction.modified.atoms, prediction.modified.aminoAcids, mod )

        #print mod
        f = calcKabschRMSD( refPosAtoms, mod )
        #print rt

        #print "tr", translation
        #print transformation
        '''aligner = PDBAligner()
        solution = aligner.transform( copy.deepcopy( modPosAtoms ), translation, rotation )

        f = aligner.calcRMSD( refPosAtoms, solution )
        #print "f", f

        ta = tAligner()
        f = ta.rmsd( transformation, refPosAtoms, copy.deepcopy( modPosAtoms ) )'''
        # calculate values for other responses
        res = {'x1':f, 'x2':f, 'x3':f,'x4':f, 'x5':f,'x6':f,'x7':f,'x8':f,'x9':f,'x10':f }
        fitness = dict(Obj=f,**res)
        return fitness

    def mp_evaluator( self, x, refPosAtoms, modPosAtoms ):
        '''Multiprocessing evaluation'''
        # ste number of cpus
        nprocs = 4
        # create pool
        pool = multiprocessing.Pool(processes=nprocs)
        results = [pool.apply_async(evaluator,[c,refPosAtoms, modPosAtoms]) for c in x]
        pool.close()
        pool.join()
        f = [r.get()['Obj'] for r in results]
        for r in results:
            del r.get()['Obj']
        # maximization or minimization problem
        maximize = False
        return (f, [r.get() for r in results],maximize)

    def initialize( self, ants,var):
        '''Create initial solution matrix'''
        X = np.random.uniform(low=0,high=1,size=(ants,var))
        '''for i in xrange( ants ):
            X[i][0] = 0.5;
            X[i][1] = 0.5;
            X[i][2] = 0.5;
            X[i][3] = 0;
            X[i][4] = 0;
            X[i][5] = 0;
        print "X", X'''
        return X

    def evolve( self, refPosAtoms, modPosAtoms, display ):
        '''Executes the optimization'''
        self.generations = []
        self.values = []
        start_time = time()

        # number of variables
        parameters_v = ['psi1', 'phi2', 'psi2', 'phi3', 'psi3', 'phi4', 'psi4', 'phi5']
        response_v = ['x1','x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8']

        # create output file
        projdir = os.getcwd()
        ind_file_name = '{0}/results.csv'.format(projdir)
        ind_file = open(ind_file_name, 'w')

        # number of variables
        nVar = len(parameters_v)
        # size of solution archive
        nSize = 200
        # number of ants
        nAnts = 150

        # parameter q
        q = 0.0001

        # standard deviation
        qk = q*nSize

        # parameter xi (like pheromone evaporation)
        xi = 0.85

        # maximum iterations
        maxiter = 10
        # tolerance
        errormin = 0.01

        # bounds of variables
        Up = [1]*nVar
        Lo = [0]*nVar

        # initilize matrices
        S = np.zeros((nSize,nVar))
        S_f = np.zeros((nSize,1))

        plt.figure()

        # initialize the solution table with uniform random distribution and sort it
        print '-----------------------------------------'
        print 'Starting initilization of solution matrix'
        print '-----------------------------------------'

        Srand = self.initialize(nSize,nVar)
        #print "Inicial", Srand
        f,S_r,maximize = self.mp_evaluator(Srand, refPosAtoms, modPosAtoms)
        #print f

        S_responses = []

        for i in range(len(S_r)):
            S_f[i] = f[i]
            k = S_r[i]
            row = []
            for r in response_v:
                row.append(k[r])
            S_responses.append(row)

        # add responses and "fitness" column to solution
        S = np.hstack((Srand,S_responses,S_f))
        # sort according to fitness (last column)
        S = sorted(S, key=lambda row: row[-1],reverse = maximize)
        S = np.array(S)
        #init_observer(ind_file,S,parameters_v,response_v)

        # initilize weight array with pdf function
        w = np.zeros((nSize))
        for i in range(nSize):
            w[i] = ( 1/( qk * math.sqrt( 2*math.pi ) ) )*math.exp( -math.pow( i, 2 )/( 2*math.pow( q, 2 )*math.pow( nSize, 2 ) ) )
            #print i, w[i]

        if display:
            x = []
            y = []
            for i in S:
                x.append(i[0])
                y.append(i[1])

            plt.scatter(x,y)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.pause(2)
            plt.cla()

        # initialize variables
        iterations = 1
        best_par = []
        best_obj = []
        best_sol = []
        best_res = []
        worst_obj = []
        best_par.append(S[0][:nVar])
        best_obj.append(S[0][-1])
        best_sol.append(S[0][:])
        best_res.append(S[0][nVar:-1])
        worst_obj.append(S[-1][-1])

        p = w/sum(w)
        stop = 0

        # iterations
        while True:
            print '-----------------------------------------'
            print 'Iteration', iterations
            print '-----------------------------------------'

            Stemp = np.zeros((nAnts,nVar))

            # para cada formiga
            for k in range( nAnts ):
                # escolhe uma solucao
                cs = np.random.random_sample()
                #print cs
                total = 0
                for z in xrange( nSize-1, -1, -1 ):
                    total += p[z]
                    if cs <= total:
                        sol = z
                        break

                # para cada variavel
                for i in range( nVar ):
                    # calcula desvio padrao da solucao 'sol'
                    sigma = 0
                    for y in xrange( nSize ):
                        sigma += abs( S[y][i] - S[sol][i] )/( nSize-1 )

                    #print "sigma", i, "=", sigma

                    # calcula valor de 'i' com as funcoes gaussianas
                    x = np.random.random_sample()
                    gi = w[sol]*math.exp( -math.pow( x - S[sol][i], 2 ) / ( 2*math.pow( sigma, 2 ) ) )* (1/( sigma*math.pow( 2*math.pi, 2 ) ))

                    #print gi
                    #Stemp[k][i] = sigma[i] * x + S[selection][i]
                    Stemp[k][i] = gi
                    if Stemp[k][i] > Up[i]:
                        Stemp[k][i] = Up[i]
                    elif Stemp[k][i] < Lo[i]:
                        Stemp[k][i] = Lo[i]

            #print Stemp
            f,S_r,maximize = self.mp_evaluator(Stemp, refPosAtoms, modPosAtoms)
            #print f

            S_f = np.zeros((nAnts,1))
            S_responses = []

            for i in range(len(S_r)):
                S_f[i] = f[i]
                k = S_r[i]
                row = []
                for r in response_v:
                    row.append(k[r])
                S_responses.append(row)

            # add responses and "fitness" column to solution
            Ssample = np.hstack((Stemp,S_responses,S_f))

            # add new solutions in the solutions table
            Solution_temp = np.vstack((S,Ssample))

            # sort according to "fitness"
            Solution_temp = sorted(Solution_temp, key=lambda row: row[-1],reverse = maximize)
            Solution_temp = np.array(Solution_temp)

            # keep best solutions
            S = Solution_temp[:nSize][:]

            #print "S", S
            # keep best after each iteration
            best_par.append(S[0][:nVar])
            best_obj.append(S[0][-1])
            best_res.append(S[0][nVar:-1])
            best_sol.append(S[0][:])
            worst_obj.append(S[-1][-1])

            #print "Best individual:", parameters_v
            #print best_sol[0][0:len(parameters_v)]
            print "Fitness:", S[0][:][10]
            self.generations.append( iterations )
            self.values.append( S[0][:][10] )

            #print best_sol

            iterations += 1
            if iterations > maxiter or stop > 5:
                break

        ind_file.close()

        total_time_s = time() - start_time
        total_time = datetime.timedelta(seconds=total_time_s)
        #total_time = formatTD(total_time)

        # fix varibales values in output file
        #correct_par(ind_file_name,parameters_v)

        best_sol = sorted( best_sol, key=lambda row: row[-1], reverse = maximize )

        print "Best individual:", parameters_v
        print best_sol[0][0:len(parameters_v)]
        print "Fitness:"
        print best_sol[0][-1]
        print "Responses:", response_v
        print best_sol[0][len(parameters_v):-1]

        print self.generations
        print self.values

        '''plt.plot( values )
        #plt.title('RMSD over iterations: ' + reference_file.replace(".pdb", "-F"))
        plt.xlabel('Generations')
        plt.ylabel('RMSD')
        plt.xlim( 0, 1000 )
        fig = plt.gcf()
        #fig.set_size_inches(10, 10)
        fig.savefig( 'acor.png', dpi=100, bbox_inches='tight')'''

        rotation = [ (2*math.pi*i)-math.pi for i in best_sol[0][0:len(parameters_v)] ]
        rotation.append( 0.0 )

        rt = []
        rt.append( 0.0 )
        for i in xrange( len( rotation ) ):
            rt.append( rotation[i] )

        #print rt
        mod = copy.deepcopy( modPosAtoms )

        rotate_to( rt, prediction.modified.atoms, prediction.modified.aminoAcids, mod )

        pdbNew = open( "1PLX-F.pdb", "w" )
        countTotal = 1
        acid = 0
        aa = None
        for z in range( 0, len( prediction.modified.atoms ) ):
            if prediction.modified.aminoAcids[z] != aa:
                aa = prediction.modified.aminoAcids[z]
                acid += 1
            pdbNew.write( self.pdbPattern.format( "ATOM", countTotal, str( prediction.modified.atoms[z] ), " ", str( prediction.modified.aAcids[z] ), " ", \
                          acid, " ", float( mod[z][0] ), float( mod[z][1] ), float( mod[z][2] ), float( 1.00 ), float( 0.0 ) ) + "\n" )

            countTotal += 1

        pdbNew.write( "TER\n" )
        pdbNew.close()