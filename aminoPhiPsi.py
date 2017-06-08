from backbone import Backbone
import copy

class AminoPhiPsi:
	filename = None
	ATOM_TAG = "ATOM"
	END_TAG = "TER"
	ALPHA_TAG = "CA"
	CARBON_TAG = "C"
	NITROGEN_TAG = "N"

	dicContent = {}

	def __init__( self, filename ):
		self.filename = filename
		self.readFile()
		self.calcPhiPsiAngles()

	def readFile( self ):
		if self.filename is not None:
			pdb = open( self.filename, "r" )
			stop = False
			key = 0
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
						if line[4] != key:
							if key > 0:
								self.dicContent[key] = None
								self.dicContent[key] = backbone
							key = line[4]
							backbone = Backbone()

						backbone.setPosAtom( line[2], map( float, line[5:8] ) )

			self.dicContent[key] = None
			self.dicContent[key] = backbone

	def calcPhiPsiAngles( self ):
		for i in range ( 1, len( self.dicContent ) + 1 ):
			#print i, self.dicContent.get( str( i ) ).getPosN()
			if i > 1:
				phi = self.calcDihedralAngle( 1, 2, 3, 4 )
				print i, "phi", phi
			if i < len( self.dicContent ):
				psi = self.calcDihedralAngle( 1, 2, 3, 4 )
				print i, "psi", psi

	def calcDihedralAngle( self, atom1, atom2, atom3, atom4 ):
		angle = 0
		return angle