/***********************************************************************
CreateProtein - Functions to create protein structures from amino acid
sequences and secondary structure prediction sequences.
Copyright (c) 2002-2020 Oliver Kreylos
Copyright (c) 2002 Nelson Max
***********************************************************************/

#ifndef CREATEPROTEIN_INCLUDED
#define CREATEPROTEIN_INCLUDED

/* Forward declarations: */
namespace Misc {
class ConfigurationFileSection;
}
namespace MD {
class Protein;
}

namespace MD {

void initializeProteinCreation(const Misc::ConfigurationFileSection& configFileSection); // Initializes protein creation module by reading given configuration file section
void setCreateTerminatorResidues(bool newCreateTerminatorResidues); // Sets the flag for creation of terminator residues (ACE and NME)
void ReadStandards(const char* standardsDirectory); // Reads standard configuration of amino acid residues
Protein* ReadPredictionFile(const char* predictionFilename); // Reads prediction file and returns protein structure
Protein* SetDihedrals(int numResidues,const int type[],const char pred[],const double phis[],const double psis[]); // Creates protein structure from given amino acid sequence and dihedral angles
Protein* loadProtein(const char* inputFileName); // Creates protein from input file; file type is determined by extension
}

#endif
