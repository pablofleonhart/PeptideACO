from mapAminoAcids import AminoAcids
import sys

sequence = str( sys.argv[1:] )
sequence = sequence.replace( "[", "" )
sequence = sequence.replace( "]", "" )
sequence = sequence.replace( "'", "" )

aminoAcids = AminoAcids( sequence );
aminoAcids.generatePDB()