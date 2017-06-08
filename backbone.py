import copy

class Backbone:
	NITROGEN_TAG = "N"
	ALPHA_TAG = "CA"
	CARBON_TAG = "C"
	posN = []
	posCA = []
	posC = []

	def setPosAtom( self, atom, positions ):
		if atom == self.NITROGEN_TAG:
			self.posN = positions
		elif atom == self.ALPHA_TAG:
			self.posCA = positions
		elif atom == self.CARBON_TAG:
			self.posC = positions

	def getPosN( self ):
		return copy.deepcopy( self.posN )

	def getPosCA( self ):
		return copy.deepcopy( self.posCA )

	def getPosC( self ):
		return copy.deepcopy( self.posC )