from atom import Atom
from backbone import Backbone
import copy

class PDBReader:
	ALPHA_CARBON = "CA"
	BACKBONE_ATOMS = ( "N", "CA", "C", "O" )
	ATOM_TAG = "ATOM"
	END_TAG = "TER"

	fileName = ""
	atoms = []
	aminoAcids = []
	posAtoms = []
	backbone = []
	alpha = []
	content = []
	aAcids = []
	dicContent = {}

	def __init__( self, fileName ):
		self.fileName = fileName
		self.readFile()

	def readFile( self ):
		self.atoms = []
		self.aminoAcids = []
		self.posAtoms = []
		self.aAcids = []
		backbone = None
		key = 0
		index = 0

		finish = False
		file = open( self.fileName, "r" )

		while not finish:
			line = file.readline()

			if not line:
				finish = True

			else:
				atom = Atom( line )

				if atom.tag.strip() == self.END_TAG:
					finish = True

				elif atom.tag.strip() == self.ATOM_TAG:
					self.atoms.append( atom.atom )
					self.posAtoms.append( atom.getPos() )
					self.aminoAcids.append( int( atom.seqResidue ) )
					self.content.append( atom )
					self.aAcids.append( atom.residue )

					if atom.seqResidue != key:
						if key > 0:
							self.dicContent[str( index )] = None
							self.dicContent[str( index )] = backbone
							index += 1
						key = atom.seqResidue
						backbone = Backbone()

					backbone.setPosAtom( atom.getAtom(), atom.getResidue(), atom.getPos() )

			self.dicContent[str( index )] = None
			self.dicContent[str( index )] = backbone

		file.close()

	def adjustAtoms( self, refAtoms, refAminoAcids ):
		auxAtoms = []
		auxPos = []
		auxAminoAcids = []
		auxContent = []
		auxAA = []

		for i in xrange( len( refAtoms ) ):  
			if ( refAtoms[i], refAminoAcids[i] ) in zip( self.atoms, self.aminoAcids ):
				index = zip( self.atoms, self.aminoAcids ).index( ( refAtoms[i], refAminoAcids[i] ) )

				auxAtoms.append( self.atoms.pop( index ) )
				auxPos.append( self.posAtoms.pop( index ) )
				auxAminoAcids.append( self.aminoAcids.pop( index ) )
				auxContent.append( self.content.pop( index ) )
				auxAA.append( self.aAcids.pop( index ) )

		self.atoms = copy.deepcopy( auxAtoms )
		self.posAtoms = copy.deepcopy( auxPos )
		self.aminoAcids = copy.deepcopy( auxAminoAcids )
		self.aAcids = copy.deepcopy( auxAA )

	def calcBackbonePos( self ):
		self.backbone = []
		print self.atoms
		for i in range( len( self.atoms ) ):
			print self.atoms[i].strip()
			if self.atoms[i].strip() in self.BACKBONE_ATOMS:
				self.backbone.append( self.posAtoms[i] )

	def calcCaPos( self ):
		self.alpha = []
		for i in range( len( self.atoms ) ):
			if self.atoms[i].strip() == self.ALPHA_CARBON:
				self.alpha.append( self.posAtoms[i] )