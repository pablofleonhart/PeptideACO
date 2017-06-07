import copy

class AminoAcids( object ):

	pdbPattern = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s} {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
	ATOM_TAG = "ATOM"
	END_TAG = "TER"
	ATOM_REF = "CA"
	OXIGEN_CARBOXYL = "OC"
	HYDROGEN_CARBOXYL = ( "HC", "HOC" )
	HYDROGEN_AMINO = ( "H", "1H" )
	NITROGEN = "N"

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
	dicContent = {}
	dicResults = []

	def __init__( self, sequence ):
		self.sequence = sequence

	def generatePDB( self ):
		self.readSequence()
		print self.dicPositions
		print self.dicAtoms
		print self.dicContent.get( "V" )[1]
		self.combineAminoAcids()

	def readSequence( self ):
		for i in range( 0, len( self.sequence ) ):
			fileName = None
			
			if self.sequence[i] is not None:
				fileName = self.dic.get( self.sequence[i] )

			if fileName is not None:
				if self.dicPositions.get( self.sequence[i] ) is None or self.dicAtoms.get( self.sequence[i] ) is None:
					atoms = []
					posAtoms = []
					translation = []
					atom = None
					pdb = open( fileName, "r" )

					stop = False
					content = []

					while not stop:
						line = pdb.readline()
						if not line:
							stop = True
						else:
							line = line.split()							
							if line[0] == self.END_TAG:
								stop = True
								content.append( line )
							elif line[0] == self.ATOM_TAG:
								content.append( line )
								atom = line[2]				
								pos = map( float, line[5:8] )

								if i == 0 and atom == self.ATOM_REF:
									translation = copy.deepcopy( pos )
									for j in xrange( 3 ):
										translation[j] *= -1

								atoms.append( atom )
								posAtoms.append( pos )

					pdb.close()
					print atoms
					print posAtoms
					print translation

					if i == 0:
						self.translateAtoms( posAtoms, translation )

					print self.sequence[i]
					self.dicPositions[self.sequence[i]] = []
					self.dicPositions[self.sequence[i]] = copy.deepcopy( posAtoms )

					self.dicAtoms[self.sequence[i]] = []
					self.dicAtoms[self.sequence[i]] = copy.deepcopy( atoms )

					self.dicContent[self.sequence[i]] = []
					self.dicContent[self.sequence[i]] = content

			'''pdb = open( fileName, "r" )
			pdbNew = open( "newPDB.pdb", "w" )

			stop = False

			while not stop:
				line = pdb.readline()
				if not line:
					stop = True
				else:
					line = line.split()

					if line[0] == self.END_TAG:
						pdbNew.write( "TER" )
						stop = True

					elif line[0] == self.ATOM_TAG:
						index = int( line[1] ) - 1

						pdbNew.write( self.pdbPattern.format( line[0], int( line[1] ), line[2], " ", line[3], " ", \
									  int( line[4] ), " ", self.posAtoms[index][0], self.posAtoms[index][1], self.posAtoms[index][2], float( line[8] ), float( line[9] ) ) + "\n" )

			pdb.close()
			pdbNew.close()'''

	def translateAtoms( self, posAtoms, translation ):
		for i in range( len( posAtoms ) ):
			for j in range( 3 ):
				posAtoms[i][j] += translation[j]

	def matchAminoAcids( self ):
		for i in range( 0, len( self.sequence ) ):
			key = self.sequence[i]
			posOC = []
			posHC = []
			posHA = []
			posN = []
			indexOC = 0
			indexHC = 0
			indexHA = 0
			indexN = 0
			keyContent = copy.deepcopy( self.dicContent.get( key ) )

			if self.dicAtoms.get( key ) is not None:
				for ( atom, pos ) in zip( self.dicAtoms[key], self.dicPositions[key] ):
					if atom == self.OXIGEN_CARBOXYL and i < len( self.sequence )-1:
						indexOC = zip( self.dicAtoms[key], self.dicPositions[key] ).index( ( atom, pos ) )
						#keyContent.pop( index )
						posOC = pos

					elif atom in self.HYDROGEN_CARBOXYL and i < len( self.sequence )-1:
						indexHC = zip( self.dicAtoms[key], self.dicPositions[key] ).index( ( atom, pos ) )
						#keyContent.pop( index )
						posHC = pos

					if atom in self.HYDROGEN_AMINO and i > 0:
						indexHA = zip( self.dicAtoms[key], self.dicPositions[key] ).index( ( atom, pos ) )
						posHA = pos

			print posOC, indexOC
			print posHC, indexHC
			print posHA, indexHA
			#print keyContent

	def combineAminoAcids( self ):
		for i in range( 0, len( self.sequence )-1 ):
			j = i+1
			keyI = self.sequence[i]
			keyJ = self.sequence[j]
			jPosOC = []
			jPosHC = []
			jPosHA = []
			jPosN = []
			jIndexOC = 0
			jIndexHC = 0
			jIndexHA = 0
			jIndexN = 0
			iPosOC = []
			iPosHC = []
			iPosHA = []
			iPosN = []
			iIndexOC = 0
			iIndexHC = 0
			iIndexHA = 0
			iIndexN = 0

			if i == 0:
				keyContentI = copy.deepcopy( self.dicContent.get( keyI ) )

			if self.dicAtoms.get( keyI ) is not None:
				for ( atom, pos ) in zip( self.dicAtoms[keyI], self.dicPositions[keyI] ):
					if atom == self.OXIGEN_CARBOXYL and i < len( self.sequence )-1:
						iIndexOC = zip( self.dicAtoms[keyI], self.dicPositions[keyI] ).index( ( atom, pos ) )
						#keyContent.pop( index )
						iPosOC = pos

					elif atom in self.HYDROGEN_CARBOXYL and i < len( self.sequence )-1:
						iIndexHC = zip( self.dicAtoms[keyI], self.dicPositions[keyI] ).index( ( atom, pos ) )
						#keyContent.pop( index )
						iPosHC = pos

					if atom in self.HYDROGEN_AMINO and i > 0:
						iIndexHA = zip( self.dicAtoms[keyI], self.dicPositions[keyI] ).index( ( atom, pos ) )
						iPosHA = pos

			if self.dicAtoms.get( keyJ ) is not None:
				for ( atom, pos ) in zip( self.dicAtoms[keyJ], self.dicPositions[keyJ] ):
					if atom == self.OXIGEN_CARBOXYL and i < len( self.sequence )-1:
						jIndexOC = zip( self.dicAtoms[keyJ], self.dicPositions[keyJ] ).index( ( atom, pos ) )
						#keyContent.pop( index )
						jPosOC = pos

					elif atom in self.HYDROGEN_CARBOXYL and i < len( self.sequence )-1:
						jIndexHC = zip( self.dicAtoms[keyJ], self.dicPositions[keyJ] ).index( ( atom, pos ) )
						#keyContent.pop( index )
						jPosHC = pos

					if atom in self.HYDROGEN_AMINO and i > 0:
						jIndexHA = zip( self.dicAtoms[keyJ], self.dicPositions[keyJ] ).index( ( atom, pos ) )
						jPosHA = pos

					if atom == self.NITROGEN:
						jIndexN = zip( self.dicAtoms[keyJ], self.dicPositions[keyJ] ).index( ( atom, pos ) )
						jPosN = pos

			#calculate translation of jN to iOC
			if iPosOC > iPosHC:
				keyContentI.pop( iIndexOC )
				keyContentI.pop( iIndexHC )
			else:
				keyContentI.pop( iIndexHC )
				keyContentI.pop( iIndexOC )

			translation = [ 0, 0, 0 ]
			print iPosOC
			print jPosN
			for k in xrange( 3 ):
				print k
				translation[k] = iPosOC[k] - jPosN[k]

			print "translation",translation

			positionsJ = copy.deepcopy( self.dicPositions[keyJ] )
			print positionsJ
			self.translateAtoms( positionsJ, translation )
			print positionsJ

			self.dicResults.append( keyContentI )
			print self.dicResults

			keyContentJ = copy.deepcopy( self.dicContent.get( keyI ) )