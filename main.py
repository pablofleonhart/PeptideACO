from mapAminoAcids import AminoAcids
from aminoPhiPsi import AminoPhiPsi
import sys

sequence = raw_input( "- Enter the desired aminoacid sequence or use the default (VSCEDCPEHCSTQKAQAKCDNDKCVCEPI) by pressing 'Enter': " )

if len( sequence ) == 0:
	sequence = "VSCEDCPEHCSTQKAQAKCDNDKCVCEPI"

if len( sequence ) > 1:
	filename = "results.pdb"
	aminoAcids = AminoAcids( sequence, filename )
	aminoPhiPsi = AminoPhiPsi( filename )
else:
	print "You must need specify at least two amino acids!"