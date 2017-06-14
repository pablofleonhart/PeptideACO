from mapAminoAcids import AminoAcids
from aminoPhiPsi import AminoPhiPsi
import sys

if len( sys.argv ) <= 1:
	sequence = "VSCEDCPEHCSTQKAQAKCDNDKCVCEPI"
else:
	sequence = str( sys.argv[1:] )
	sequence = sequence.replace( "[", "" )
	sequence = sequence.replace( "]", "" )
	sequence = sequence.replace( "'", "" )

if len( sequence ) > 1:
	aminoAcids = AminoAcids( sequence )
	aminoPhiPsi = AminoPhiPsi( "results.pdb" )
else:
	print "You must need specify at least two amino acids!"