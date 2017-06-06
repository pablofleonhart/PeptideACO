import copy

class AminoAcids( object ):

	ATOM_TAG = "ATOM"
	END_TAG = "TER"
	ATOM_REF = "CA"

	dic = { "A": "files/alanine.pdb",
			"R": "files/arginine.pdb",
			"N": "files/asparagine.pdb",
			"D": "files/aspartic_acid.pdb",
			"C": "files/cysteine.pdb",
			"E": "files/glutamic_acid.pdb",
			"Q": "files/glutamine.pdb",
			"G": "files/glycine.pdb",
			"H": "files/histidine.pdb",
			"I": "files/isoleucine.pdb",
			"L": "files/leucine.pdb",
			"K": "files/lysine.pdb",
			"M": "files/methionine.pdb",
			"F": "files/phenalalanine.pdb",
			"P": "files/proline.pdb",
			"S": "files/serine.pdb",
			"T": "files/threonine.pdb",
			"W": "files/tryptophan.pdb",
			"Y": "files/tyrosine.pdb",
			"V": "files/valine.pdb" }

	dicPositions = {}
	dicAtoms = {}
	atoms = []
	posAtoms = []
	translation = []

	def __init__( self, sequence ):
		self.sequence = sequence

	def generatePDB( self ):
		fileName = None

		for i in range( 0, len( self.sequence ) ):
			if self.sequence[i] is not None:
				fileName = self.dic.get( self.sequence[i] )

			if fileName is not None:
				self.atoms = []
				self.posAtoms = []
				self.translation = []
				atom = None
				pdb = open( fileName, "r" )

				stop = False

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
							pos = map( float, line[5:8] )

							if atom == self.ATOM_REF:
								self.translation = copy.deepcopy( pos )
								for j in xrange( 3 ):
									self.translation[j] *= -1

							self.atoms.append( atom )
							self.posAtoms.append( pos )

				pdb.close()
				print self.atoms
				print self.posAtoms
				print self.translation

				self.translateAtoms()
				print self.sequence[i]
				self.dicPositions[self.sequence[i]] = []
				self.dicPositions[self.sequence[i]].append( self.posAtoms )

				self.dicAtoms[self.sequence[i]] = []
				self.dicAtoms[self.sequence[i]].append( self.atoms )

		print self.dicPositions
		print self.dicAtoms

	def translateAtoms( self ):
		for i in range( len( self.atoms ) ):
			for j in range( 3 ):
				self.posAtoms[i][j] += self.translation[j]