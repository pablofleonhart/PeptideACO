from mapAminoAcids import AminoAcids
from aminoPhiPsi import AminoPhiPsi
from pdbReader import PDBReader
from pdbAligner import PDBAligner
import math
import numpy as np
import rmsd
import sys

class Prediction( object ):
	experimental = None
	modified = None
	mod = None

	def __init__( self ):
		sequence = raw_input( "- Enter the desired aminoacid sequence or use the default (YGGFM) by pressing 'Enter': " )

		if len( sequence ) == 0:
			sequence = "YGGFM"

		if len( sequence ) > 1:
			fileName = "1PLX-P.pdb"
			aminoAcids = AminoAcids( sequence, fileName )
			print "OK - The file '1PLX-P.pdb' with the peptide bonds was generated."
			name = raw_input( "- Enter the PDB filename to calc the dihedral angles: Phi and Psi, or use the '" + "1PLX-P.pdb" + "' by pressing 'Enter':" )

			if len( name ) == 0:
				name = fileName

			aminoPhiPsi = AminoPhiPsi( name )
			print "OK - The file 'aminoPhiPsi.txt' with the dihedral angles by amino acid was generated."
			print "OK - The file 'ramachandran.png' with the ramachandran map was generated."

			print aminoPhiPsi.get_angles()
			aminoPhiPsi.rotate_omegas()
			print aminoPhiPsi.get_omegas()
			#print aminoPhiPsi.get_angles()
			'''pis = [math.radians(360.0),math.radians(176.63),
				   math.radians(148.48),math.radians(-21.96),
				   math.radians(114.02),math.radians(29.89),
				   math.radians(-88.0),math.radians(-38.16),
				   math.radians(-74.24),math.radians(360.0)]'''
			#print pis
			aminoPhiPsi.set_peptide_bond_angles()
			#print aminoPhiPsi.get_peptide_bond_angles()		
			#aminoPhiPsi.rotate_to( pis )
			print aminoPhiPsi.get_angles()
			aminoPhiPsi.writePDBFile( "1PLX-P.pdb" )
			# TODO: set phi and psi in 180 degrees
			self.readFiles( "files/1PLX.pdb", "1PLX-P.pdb" )
			self.calcRMSD()
			self.calcKabschRMSD()
			# TODO: get the file 1PLX-F.pdb from ACOR
			#calcRMSD( "files/1PLX.pdb", "1PLX-F.pdb" )
			# TODO: plot chart RMSD VS algorithm generation

		else:
			print "ERROR - You must need specify at least two amino acids!"

	def readFiles( self, fileA, fileB ):
		self.experimental = PDBReader( fileA )
		self.modified = PDBReader( fileB )

		self.modified.adjustAtoms( self.experimental.atoms, self.experimental.aminoAcids )
		self.experimental.adjustAtoms( self.modified.atoms, self.modified.aminoAcids )

		print self.experimental.atoms
		print self.experimental.posAtoms
		print self.modified.atoms
		print self.modified.posAtoms

		self.experimental.calcBackbonePos()
		self.modified.calcBackbonePos()
		self.experimental.calcCaPos()
		self.modified.calcCaPos()

	def calcKabschRMSD( self ):
		P = np.array( self.experimental.posAtoms )
		Q = np.array( self.modified.posAtoms )
		#print rmsd.kabsch_rmsd( P, Q )
		P -= rmsd.centroid( P )
		Q -= rmsd.centroid( Q )
		print "{:15s} {:6.2f}".format( "Kabsch RMSD:", rmsd.kabsch_rmsd( P, Q ) )

	def calcRMSD( self ):
		print( len( self.experimental.atoms ), len( self.modified.atoms ) )
		#print( self.experimental.atoms )
		#print( self.modified.atoms )
		print( len( self.experimental.backbone ), len( self.modified.backbone ) )
		print( len( self.experimental.alpha ), len( self.modified.alpha ) )

		aligner = PDBAligner()
		print "{:15s} {:6.2f}".format( "CA RMSD:", aligner.calcRMSD( self.experimental.alpha, self.modified.alpha ) )
		print "{:15s} {:6.2f}".format( "Backbone RMSD:", aligner.calcRMSD( self.experimental.backbone, self.modified.backbone ) )
		print "{:15s} {:6.2f}".format( "All atoms RMSD:", aligner.calcRMSD( self.experimental.posAtoms, self.modified.posAtoms ) )

	def setModPosAtoms( self, mod ):
		self.mod = mod
