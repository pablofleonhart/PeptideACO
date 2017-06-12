from mapAminoAcids import AminoAcids
from aminoPhiPsi import AminoPhiPsi
from pdbReader import PDBReader
from pdbAligner import PDBAligner
import sys

def calcRMSD( fileA, fileB ):
	experimental = PDBReader( fileA )
	modified = PDBReader( fileB )

	aligner = PDBAligner()
	print "{:15s} {:6.2f}".format( "CA RMSD:", aligner.calcRMSD( experimental.alpha, modified.alpha ) )
	#print "{:15s} {:6.2f}".format( "Backbone RMSD:", aligner.calcRMSD( experimental.backbone, modified.backbone ) )
	#print "{:15s} {:6.2f}".format( "All atoms RMSD:", aligner.calcRMSD( experimental.posAtoms, modified.posAtoms ) )

sequence = raw_input( "- Enter the desired aminoacid sequence or use the default (YGGFM) by pressing 'Enter': " )

if len( sequence ) == 0:
	sequence = "YGGFM"

if len( sequence ) > 1:
	aminoAcids = AminoAcids( sequence, "1PLX-P.pdb" )
	print "OK - The file '1PLX-P.pdb' with the peptide bonds was generated."
	filename = raw_input( "- Enter the PDB filename to calc the dihedral angles: Phi and Psi, or use the '1PLX-P.pdb' by pressing 'Enter':" )

	if len( filename ) == 0:
		filename = "1PLX-P.pdb"

	aminoPhiPsi = AminoPhiPsi( filename )
	print "OK - The file 'aminoPhiPsi.txt' with the dihedral angles by amino acid was generated."
	print "OK - The file 'ramachandran.png' with the ramachandran map was generated."

	# TODO: set phi and psi in 180 degrees
	calcRMSD( "files/1PLX.pdb", "1PLX-P.pdb" )
	# TODO: get the file 1PLX-F.pdb from ACOR
	#calcRMSD( "files/1PLX.pdb", "1PLX-F.pdb" )
	# TODO: plot chart RMSD vs algorithm generation

else:
	print "ERROR - You must need specify at least two amino acids!"