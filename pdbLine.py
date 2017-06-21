class PDBLine( object ):
	tag = 'ATOM'
	sequence = 0
	atom = ''
	aminoAcid = ''
	posX = 0
	posY = 0
	posZ = 0
	ns1 = 0
	ns2 = 0
	content = None

	def __init__( self, content ):
		self.content = content
		if content is not None:
			self.tag = content[0]
			self.sequence = content[1]
			self.atom = content[2]
			self.aminoAcid = content[3]
			self.posX = content[4]
			self.posY = content[5]
			self.posZ = content[6]
			self.ns1 = content[7]
			self.ns2 = content[8]