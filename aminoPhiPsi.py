from backbone import Backbone
import copy
import math
import matplotlib.pyplot as plt
import numpy as np

class AminoPhiPsi:
	pdbPattern = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s} {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
	filename = None
	ATOM_TAG = "ATOM"
	END_TAG = "TER"
	ALPHA_TAG = "CA"
	CARBON_TAG = "C"
	NITROGEN_TAG = "N"
	OC_ATOMS = ("C", "O", "OC", "HOC", "HC", "HO")
	NH_ATOMS = ("N", "H", "1H", "H1", "2H", "H2", "3H", "H3")
	NHC_ATOMS = ("N", "H", "1H", "H1", "2H", "H2", "3H", "H3", "CA")

	dicContent = {}
	phi = []
	psi = []
	angles = []
	omega = []
	aminoAcids = []
	aAcids = []
	posAtoms = []
	atoms = []
	dicResults = {}

	def __init__( self, filename ):
		self.filename = filename
		self.readFile()
		self.calcAngles()
		self.plotRamanchandran()

	def isNumber( self, value ):
		try:
			int( value )
			return True

		except ValueError:
			return False

	def degrees( self, x ):
		return x*180/math.pi

	def rad( self, x ):
		return math.pi*x/180

	def rotaxis2m( self, theta, vector ): 
		vector=vector.copy() 
		vector.normalize() 
		c=numpy.cos(theta) 
		s=numpy.sin(theta) 
		t=1-c 
		x, y, z=vector.get_array() 
		rot=numpy.zeros((3, 3)) 
		# 1st row 
		rot[0, 0]=t*x*x+c 
		rot[0, 1]=t*x*y-s*z 
		rot[0, 2]=t*x*z+s*y 
		# 2nd row 
		rot[1, 0]=t*x*y+s*z 
		rot[1, 1]=t*y*y+c 
		rot[1, 2]=t*y*z-s*x 
		# 3rd row 
		rot[2, 0]=t*x*z-s*y 
		rot[2, 1]=t*y*z+s*x 
		rot[2, 2]=t*z*z+c 
		return rot

	def AxB( self, A, B ):   
		A_linhas = len(A)
		A_colunas = len(A[0])
		B_linhas = len(B)
		B_colunas = len(B[0])
		if A_colunas == B_linhas:
			comum = A_colunas
			M = [[sum(A[m][n] * B[n][p] for n in range(comum)) \
				for p in range(B_colunas)] for m in range(A_linhas)]
			return M
		else:
			return -1

	def readFile( self ):
		if self.filename is not None:
			pdb = open( self.filename, "r" )
			stop = False
			key = 0
			index = 0
			backbone = None
			self.aminoAcids = []
			self.posAtoms = []
			self.atoms = []
			self.aAcids = []
			ii = 0

			while not stop:
				line = pdb.readline()
				if not line:
					stop = True
				else:
					line = line.split()
					if line[0] == self.END_TAG:
						stop = True
					elif line[0] == self.ATOM_TAG:
						atom = line[2]
						if self.isNumber( line[4] ):
							aminoAcid = int( line[4] )

						elif self.isNumber( line[5] ):
							aminoAcid = int( line[5] )

						posInit = 0
						self.atoms.append( atom )
						self.aminoAcids.append( aminoAcid )

						for i in xrange( len( line ) ):
							if "." in line[i]:
								posInit = i
								break

						self.dicResults[str( ii )] = line
						ii += 1

						if aminoAcid != key:
							if key > 0:
								self.dicContent[str( index )] = None
								self.dicContent[str( index )] = backbone
								index += 1
							key = aminoAcid
							backbone = Backbone()

						pos = map( float, line[posInit:posInit + 3] )
						self.posAtoms.append( pos )
						backbone.setPosAtom( line[2], line[3], pos )
						self.aAcids.append( line[3] )

			self.dicContent[str( index )] = None
			self.dicContent[str( index )] = backbone

		#print self.dicContent

	def calcAngles( self ):
		file = open( "aminoPhiPsi.txt", "w" )
		file.write( "{:5s}  {:>7s}  {:>7s}  {:>7s}".format( "Amino", "Phi", "Psi", "Omega" ) + "\n" )
		self.phi = []
		self.psi = []
		self.angles = []
		self.omega = []

		for i in range ( 0, len( self.dicContent ) ):
			phiValue = 360.00
			psiValue = 360.00
			omegaValue = 360.00

			if i > 0:
				phiValue = self.calcDihedralAngle( self.dicContent.get( str( i-1 ) ).getPosC(), self.dicContent.get( str( i ) ).getPosN(), self.dicContent.get( str( i ) ).getPosCA(), self.dicContent.get( str( i ) ).getPosC() )
			if i < len( self.dicContent )-1:
				psiValue = self.calcDihedralAngle( self.dicContent.get( str( i ) ).getPosN(), self.dicContent.get( str( i ) ).getPosCA(), self.dicContent.get( str( i ) ).getPosC(), self.dicContent.get( str( i+1 ) ).getPosN() )
				omegaValue = self.calcDihedralAngle( self.dicContent.get( str( i ) ).getPosCA(), self.dicContent.get( str( i ) ).getPosC(), self.dicContent.get( str( i+1 ) ).getPosN(), self.dicContent.get( str( i+1 ) ).getPosCA() )

			'''if phiValue < 180.0:
				diff = 180 - phiValue
				print self.rad( diff )'''

			self.phi.append( phiValue )
			self.psi.append( psiValue )
			self.angles.append( self.rad( phiValue ) )
			self.angles.append( self.rad( psiValue ) )
			self.omega.append( omegaValue )

			file.write( "{:5s}  {:7.2f}  {:7.2f}  {:7.2f}".format( self.dicContent.get( str( i ) ).getAminoAcid(), phiValue, psiValue, omegaValue ) + "\n" )

		file.close()

	def calcDihedralAngle( self, atom1, atom2, atom3, atom4 ):
		vector1 = [( atom2[0] - atom1[0] ), ( atom2[1] - atom1[1] ), ( atom2[2] - atom1[2] )]
		vector2 = [( atom3[0] - atom2[0] ), ( atom3[1] - atom2[1] ), ( atom3[2] - atom2[2] )]
		vector3 = [( atom4[0] - atom3[0] ), ( atom4[1] - atom3[1] ), ( atom4[2] - atom3[2] )]

		normal1 = [( vector1[1] * vector2[2] - vector1[2] * vector2[1] ), ( vector1[2] * vector2[0] - vector1[0] * vector2[2] ), ( vector1[0] * vector2[1] - vector1[1] * vector2[0] )]
		normal2 = [( vector2[1] * vector3[2] - vector2[2] * vector3[1] ), ( vector2[2] * vector3[0] - vector2[0] * vector3[2] ), ( vector2[0] * vector3[1] - vector2[1] * vector3[0] )]

		n1 = math.sqrt( ( normal1[0]**2) + ( normal1[1]**2 ) + ( normal1[2]**2 ) )
		n2 = math.sqrt( ( normal2[0]**2) + ( normal2[1]**2 ) + ( normal2[2]**2 ) )

		normal1 = [( normal1[0]/n1 ), ( normal1[1]/n1 ), ( normal1[2]/n1 )]
		normal2 = [( normal2[0]/n2 ), ( normal2[1]/n2 ), ( normal2[2]/n2 )]

		v2 = math.sqrt( ( vector2[0]**2 ) + ( vector2[1]**2 ) + ( vector2[2]**2 ) )
		vector2 = [( vector2[0]/v2 ), ( vector2[1]/v2 ), ( vector2[2]/v2 )]

		m1 = [( ( normal1[1] * normal2[2] ) - ( normal1[2] * normal2[1] ) ), ( ( normal1[2] * normal2[0] ) - ( normal1[0] * normal2[2] ) ), ( ( normal1[0] * normal2[1] ) - ( normal1[1] * normal2[0] ) )]

		x = ( normal1[0] * normal2[0] ) + ( normal1[1] * normal2[1] ) + ( normal1[2] * normal2[2] )
		y = ( m1[0] * vector2[0] ) + ( m1[1] * vector2[1] ) + ( m1[2] * vector2[2] )

		angle = round( math.degrees( math.atan2( y, x ) ), 2 )
		return angle

	def rotate_omegas( self, angles=[] ):
		n_aa = len( self.aminoAcids )
		#print n_aa
		if len( angles ) == 0:
			angles = [math.pi]*n_aa
			for i in xrange( n_aa ):
				if i + min( self.aminoAcids) < max(self.aminoAcids):
					#ROTATE OMEGA
					#print 'i',i, self.dicContent.get( str( i ) )
					c_i  = zip(self.atoms, self.aminoAcids).index(("C",  i + min(self.aminoAcids))) #C from aminoacid i		
					nn_i = zip(self.atoms, self.aminoAcids).index(("N", i + 1 + min(self.aminoAcids))) #N from aminoacid i+1					
					current_omegas = self.omega
					domega = math.atan2( math.sin( angles[i] - current_omegas[i]), math.cos(angles[i] - current_omegas[i]))
					#print "domega", domega
					c_pos  = self.posAtoms[c_i]
					nn_pos = self.posAtoms[nn_i]
					#print "c_i", c_i, c_pos
					#print "nn_i", nn_i, nn_pos
					ia = 0
					for atom in zip(self.atoms, self.aminoAcids):
						if (atom[1] > i + 1 + min(self.aminoAcids) or (atom[1] == i + 1 + min(self.aminoAcids) and (atom[0] != "N"))): 
							self.posAtoms[ia] = self.rotate_atom_around_bond(domega, self.posAtoms[ia], c_pos, nn_pos)
							'''print self.atoms[ia], self.aAcids[ia], self.posAtoms[ia]
							print i
							self.dicContent.get( str( i ) ).setPosAtom( self.atoms[ia], self.aAcids[ia], self.posAtoms[ia] )'''
						ia += 1

	def rotate_to( self, ang ):
		#print "ang", ang
		n_aa = len( self.aminoAcids )
		for i in xrange(n_aa):
			if i + min( self.aminoAcids) <= max(self.aminoAcids):
			#ROTATE PHI
				#print self.atoms, self.aminoAcids
				n_i = zip(self.atoms, self.aminoAcids).index(("N", i + min(self.aminoAcids)))   
				ca_i = zip(self.atoms, self.aminoAcids).index(("CA", i + min(self.aminoAcids)))
				current_angles = self.angles
				#print current_angles
				dphi = math.atan2(math.sin(ang[2*i] - current_angles[2*i]), math.cos(ang[2*i] - current_angles[2*i]))
				#print "dphi", self.degrees( dphi )
				n_pos = self.posAtoms[n_i]
				ca_pos = self.posAtoms[ca_i]                
				ia = 0
				for atom in zip(self.atoms, self.aminoAcids):
					if (i > 0) and (atom[1] > i + min(self.aminoAcids) or (atom[1] == i + min(self.aminoAcids) and (atom[0] not in self.NHC_ATOMS))): 
						self.posAtoms[ia] = self.rotate_atom_around_bond(dphi, self.posAtoms[ia], n_pos, ca_pos)
						#print(atom[0], atom[1])   
					ia += 1        
				#ROTATE PSI    
				c_i  = zip(self.atoms, self.aminoAcids).index(("C",  i + min(self.aminoAcids)))  
				ca_i = zip(self.atoms, self.aminoAcids).index(("CA", i + min(self.aminoAcids)))
				current_angles = self.angles
				#print current_angles
				dpsi = math.atan2(math.sin(ang[2*i+1] - current_angles[2*i+1]), math.cos(ang[2*i+1] - current_angles[2*i+1]))              
				c_pos = self.posAtoms[c_i] 
				ca_pos = self.posAtoms[ca_i]
				ia = 0
				for atom in zip(self.atoms, self.aminoAcids):
					if (i+min(self.aminoAcids) < max(self.aminoAcids)) and (atom[1] > i+min(self.aminoAcids) or (atom[1] == i+min(self.aminoAcids) and (atom[0]=="O"))): 
						self.posAtoms[ia] = self.rotate_atom_around_bond(dpsi, self.posAtoms[ia], ca_pos, c_pos)
						#print(atom[0], atom[1])          
					ia += 1
            
	def normalize( self, v ):
		norm = np.linalg.norm( v )
		if norm == 0: 
			return v
		return v/norm  

	def rotate_atom_around_bond( self, theta, atom_pos, bond_start, bond_end ):
		v = np.array( atom_pos ) - np.array( bond_start )
		k = np.array( bond_end ) - np.array( bond_start )
		k = self.normalize( k )
		rot_pos = v * np.cos( theta ) + ( np.cross( k, v ) ) * np.sin( theta ) + k * ( np.dot( k, v ) ) * ( 1.0 - np.cos( theta ) )
		return list( rot_pos + np.array( bond_start ) )

	def get_ca_info(self):
		ca_info = []
		for a in zip(self.atoms, self.aminoAcids, self.posAtoms):
			if a[0] == self.ALPHA_TAG:
				ca_info.append(a)
		return copy.deepcopy(ca_info)

	def get_N_info(self):
		n_info = []
		for a in zip(self.atoms, self.aminoAcids, self.posAtoms):
			if a[0] == self.NITROGEN_TAG:
				n_info.append(a)
		return copy.deepcopy(n_info)

	def get_C_info(self):
		c_info = []
		for a in zip(self.atoms, self.aminoAcids, self.posAtoms):
			if a[0] == self.CARBON_TAG:
				c_info.append(a)
		return copy.deepcopy(c_info) 

	def get_omegas(self):
		angles = []        
		ca = self.get_ca_info()
		n  = self.get_N_info()
		c  = self.get_C_info()       
		for aa in xrange(len(ca)):
			if aa < len(ca) - 1:
				name       = ca[aa][1]
				ca_pos     = ca[aa][2]
				c_pos      =  c[aa][2]                
				nex_n_pos  =  n[aa+1][2]
				nex_ca_pos = ca[aa+1][2]
				omega = self.calcDihedralAngle(ca_pos, c_pos, nex_n_pos, nex_ca_pos)
				angles.append(omega)
			else:
				angles.append(2.0*math.pi)
		return angles

	def get_angles(self):
		angles = []   
		ca = self.get_ca_info()
		n  = self.get_N_info()
		c  = self.get_C_info()      
		print ca
		for aa in xrange(len(ca)):
			name   = ca[aa][1]
			if aa > 0:
				pre_c_pos = c[aa-1][2]
			n_pos  =  n[aa][2]
			ca_pos = ca[aa][2]
			c_pos  =  c[aa][2]
			if aa < len(ca) - 1:
				nex_n_pos = n[aa+1][2]       
			if aa > 0: 
				phi = self.calcDihedralAngle(pre_c_pos, n_pos, ca_pos, c_pos)
			else:
				phi = 2 * math.pi   
			if aa < len(ca) - 1: 
				psi = self.calcDihedralAngle(n_pos, ca_pos, c_pos, nex_n_pos)
			else:
				psi = 2 * math.pi 
			angles.append(phi)
			angles.append(psi)
		return angles

	def calc_angle_3(self, pos1, posC, pos2):
		pos1 = np.array(pos1)
		posC = np.array(posC)
		pos2 = np.array(pos2)
		bond1C = self.normalize(pos1 - posC)
		bond2C = self.normalize(pos2 - posC)
		dp = np.dot(bond1C, bond2C)
		angle = np.arccos(dp)
		return angle  

	def get_peptide_bond_angles(self):
		angles = []
		ca = self.get_ca_info()
		n  = self.get_N_info()
		c  = self.get_C_info()
		for aa in xrange(len(ca)):
			if aa < len(ca) - 1:
				name       = ca[aa][1]
				c_pos      =  c[aa][2]                
				nex_n_pos  =  n[aa+1][2]
				nex_ca_pos = ca[aa+1][2]
				alpha = self.calc_angle_3(c_pos, nex_n_pos, nex_ca_pos)
				angles.append(alpha)
			else:
				angles.append(2.0*math.pi)
		return angles 

	def set_peptide_bond_angles(self, angles=[]):
		n_aa = len( self.aminoAcids )
		if len(angles) == 0:
			angles = [math.radians(120.0)]*n_aa
		for i in xrange(n_aa):
			if i + min(self.aminoAcids) < max(self.aminoAcids):
				#ROTATE ALPHA
				c_i   = zip(self.atoms, self.aminoAcids).index(("C",  i + min(self.aminoAcids)))     # C from aminoacid i
				nn_i  = zip(self.atoms, self.aminoAcids).index(("N",  i + 1 + min(self.aminoAcids))) # N from aminoacid i+1
				nh_i  = -1.0
				for a in self.NH_ATOMS:
					if a != "N":                                                 
						nh_i  = zip(self.atoms, self.aminoAcids).index((a,  i + 1 + min(self.aminoAcids))) if (a,  i + 1 + min(self.aminoAcids)) in zip(self.atoms, self.aminoAcids) else -1
						if nh_i >= 0:
							break
					    
				nca_i = zip(self.atoms, self.aminoAcids).index(("CA", i + 1 + min(self.aminoAcids))) #CA from aminoacid i+1
				current_alphas = self.get_peptide_bond_angles()
				dalpha = math.atan2(math.sin(angles[i] - current_alphas[i]), math.cos(angles[i] - current_alphas[i]))
				c_pos   = self.posAtoms[c_i]
				nn_pos  = self.posAtoms[nn_i]
				if nh_i >= 0:
					nh_pos  = self.posAtoms[nh_i]
					current_alphaH = self.calc_angle_3(c_pos, nn_pos, nh_pos)
					dalphaH = math.atan2(math.sin(angles[i] - current_alphaH), math.cos(angles[i] - current_alphaH)) 
				nca_pos = self.posAtoms[nca_i]
				ia = 0
				for atom in zip(self.atoms, self.aminoAcids):
					if (atom[1] > i + 1 + min(self.aminoAcids) or (atom[1] == i + 1 + min(self.aminoAcids) and (atom[0] not in self.NH_ATOMS))): 
						self.posAtoms[ia] = self.bend_bonds(dalpha, self.posAtoms[ia], c_pos, nn_pos, nca_pos)
					elif nh_i >= 0 and atom[1] == i + 1 + min(self.aminoAcids) and atom[0] in self.NH_ATOMS and atom[0] != "N":
						self.posAtoms[ia] = self.bend_bonds(dalphaH, self.posAtoms[ia], c_pos, nn_pos, nh_pos)    
						#print(atom[0], atom[1])
					ia += 1                   

	def bend_bonds(self, theta, atom_pos, pos1, posC, pos2):
		#https://pt.stackoverflow.com/questions/25923/vetores-e-%C3%82ngulos-geometria-molecular
		posC = np.array(posC)
		atom_pos = np.array(atom_pos) - posC
		pos1 = np.array(pos1) - posC
		pos2 = np.array(pos2) - posC
		bond1C = self.normalize(pos1)
		bond2C = self.normalize(pos2)
		ortho = self.normalize(np.cross(bond1C, bond2C))

		c = np.cos(theta)
		s = np.sin(theta)
		t = 1.0 - c

		rot = np.matrix([[c + ortho[0] * ortho[0] * t, ortho[0] * ortho[1] * t - ortho[2] * s, ortho[0] * ortho[2] * t + ortho[1] * s], 
		                 [ortho[0] * ortho[1] * t + ortho[2] * s, c + ortho[1] * ortho[1] * t, ortho[1] * ortho[2] * t - ortho[0] * s],
		                 [ortho[2] * ortho[0] * t - ortho[1] * s, ortho[2] * ortho[1] * t + ortho[0] * s, c + ortho[2] * ortho[2] * t]])
		                 
		new_pos = np.matrix.tolist(np.matrix(atom_pos) * rot.transpose())[0]               
		new_pos = list(np.array(new_pos) + posC)
		return new_pos

	def writePDBFile( self, name ):
		pdbNew = open( name, "w" )
		countTotal = 1
		acid = 0
		aa = None
		for z in range( 0, len( self.atoms ) ):
			key = str( z )
			#print len( self.dicResults.get( key ) )
			#for y in range( 0, len( self.dicResults.get( key ) ) ):
				#self.dicResults.get( key )[y][1] = countTotal
				#self.dicResults.get( key )[y][4] = z+1

				#self.dicResults.get( key )[y][2] = self.renameAtom( self.dicResults.get( key )[y][2], self.dicResults.get( key )[y][3] )

			if self.dicResults.get( key )[4] != aa:
				aa = self.dicResults.get( key )[4]
				acid += 1
			pdbNew.write( self.pdbPattern.format( self.dicResults.get( key )[0], countTotal, self.dicResults.get( key )[2], " ", self.dicResults.get( key )[3], " ", \
						  acid, " ", float( self.posAtoms[z][0] ), float( self.posAtoms[z][1] ), float( self.posAtoms[z][2] ), float( self.dicResults.get( key )[8] ), float( self.dicResults.get( key )[9] ) ) + "\n" )

			countTotal += 1

		pdbNew.write( "TER\n" )
		pdbNew.close()

	def plotRamanchandran( self ):
		plt.plot( self.phi, self.psi, 'ro', color = "green", ms = 7.0 )
		plt.xlim( -180, 180 )
		plt.ylim( -180, 180 )

		plt.xticks( np.arange( -180.1, 180.1, 60 ) )
		plt.yticks( np.arange( -180.1, 180.1, 60 ) )
		plt.xlabel( "Phi(deg)" )
		plt.ylabel( "Psi(deg)" )
		plt.arrow( -180, 0, 360, 0 )
		plt.arrow( 0, -180, 0, 360 )

		fig = plt.gcf()
		fig.set_size_inches( 6.0, 6.0 )
		fig.savefig( 'ramachandran.png', dpi=300, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None )