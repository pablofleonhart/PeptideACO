from backbone import Backbone
import copy
import math
import matplotlib.pyplot as plt
import numpy as np

class AminoPhiPsi:
	filename = None
	ATOM_TAG = "ATOM"
	END_TAG = "TER"
	ALPHA_TAG = "CA"
	CARBON_TAG = "C"
	NITROGEN_TAG = "N"

	dicContent = {}
	phi = []
	psi = []

	def __init__( self, filename ):
		self.filename = filename
		self.readFile()
		self.calcPhiPsiAngles()
		self.plotRamanchandran()

	def isNumber( self, value ):
		try:
			int( value )
			return True

		except ValueError:
			return False

	def readFile( self ):
		if self.filename is not None:
			pdb = open( self.filename, "r" )
			stop = False
			key = 0
			index = 1
			backbone = None

			while not stop:
				line = pdb.readline()
				if not line:
					stop = True
				else:
					line = line.split()
					if line[0] == self.END_TAG:
						stop = True
					elif line[0] == self.ATOM_TAG:
						if self.isNumber( line[4] ):
							aminoAcid = int( line[4] )

						elif self.isNumber( line[5] ):
							aminoAcid = int( line[5] )

						posInit = 0

						for i in xrange( len( line ) ):
							if "." in line[i]:
								posInit = i
								break

						if aminoAcid != key:
							if key > 0:
								self.dicContent[str( index )] = None
								self.dicContent[str( index )] = backbone
								index += 1
							key = aminoAcid
							backbone = Backbone()

						backbone.setPosAtom( line[2], line[3], map( float, line[posInit:posInit + 3] ) )

			self.dicContent[str( index )] = None
			self.dicContent[str( index )] = backbone

	def calcPhiPsiAngles( self ):
		file = open( "aminoPhiPsi.txt", "w" )
		file.write( "{:5s}  {:>7s}  {:>7s}".format( "Amino", "Phi", "Psi" ) + "\n" )
		self.phi = []
		self.psi = []

		for i in range ( 1, len( self.dicContent )+1 ):
			phiValue = 360.00
			psiValue = 360.00

			if i > 1:
				phiValue = self.calcDihedralAngle( self.dicContent.get( str( i-1 ) ).getPosC(), self.dicContent.get( str( i ) ).getPosN(), self.dicContent.get( str( i ) ).getPosCA(), self.dicContent.get( str( i ) ).getPosC() )
			if i < len( self.dicContent ):
				psiValue = self.calcDihedralAngle( self.dicContent.get( str( i ) ).getPosN(), self.dicContent.get( str( i ) ).getPosCA(), self.dicContent.get( str( i ) ).getPosC(), self.dicContent.get( str( i+1 ) ).getPosN() )

			self.phi.append( phiValue )
			self.psi.append( psiValue )

			file.write( "{:5s}  {:7.2f}  {:7.2f}".format( self.dicContent.get( str( i ) ).getAminoAcid(), phiValue, psiValue ) + "\n" )

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

	def plotRamanchandran( self ):
		plt.plot( self.phi, self.psi, 'ro', color = "green", ms = 7.0 )
		plt.xlim( -180, 180 )
		plt.ylim( -180, 180 )

		plt.xticks( np.arange( -180.1, 180.1, 30 ) )
		plt.yticks( np.arange( -180.1, 180.1, 30 ) )
		plt.xlabel( "Phi(deg)" )
		plt.ylabel( "Psi(deg)" )
		plt.arrow( -180, 0, 360, 0 )
		plt.arrow( 0, -180, 0, 360 )

		fig = plt.gcf()
		fig.set_size_inches( 10.0, 10.0 )
		fig.savefig( 'rama.png', dpi=300, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None )