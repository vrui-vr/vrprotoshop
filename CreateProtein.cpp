/***********************************************************************
CreateProtein - Functions to create protein structures from amino acid
sequences and secondary structure prediction sequences.
Copyright (c) 2002-2024 Oliver Kreylos
Copyright (c) 2002 Nelson Max
***********************************************************************/

#include <string.h>
#include <vector>
#include <Misc/StdError.h>
#include <Misc/StandardValueCoders.h>
#include <Misc/ConfigurationFile.h>
#include <Math/Math.h>

#include "Protein.h"
#include "ParsePdbFile.h"
#include "CreateProtein.h"

namespace MD {

double proteinCreatorHphi;
double proteinCreatorHpsi;
double proteinCreatorEphi;
double proteinCreatorEpsi;
double proteinCreatorCphi;
double proteinCreatorCpsi;
double proteinCreatorProlinePhi;
bool proteinCreatorCreateTerminatorResidues;

void initializeProteinCreation(const Misc::ConfigurationFileSection& configFileSection)
	{
	/* Read default angles for alpha helices: */
	proteinCreatorHphi=Math::rad(configFileSection.retrieveValue<double>("./alphaHelixPhi",-60.0));
	proteinCreatorHpsi=Math::rad(configFileSection.retrieveValue<double>("./alphaHelixPsi",-40.0));
	
	/* Read default angles for beta strands: */
	proteinCreatorEphi=Math::rad(configFileSection.retrieveValue<double>("./betaStrandPhi",-120.0));
	proteinCreatorEpsi=Math::rad(configFileSection.retrieveValue<double>("./betaStrandPsi",140.0));
	
	/* Read default angles for coil regions: */
	proteinCreatorCphi=Math::rad(configFileSection.retrieveValue<double>("./coilPhi",-182.0));
	proteinCreatorCpsi=Math::rad(configFileSection.retrieveValue<double>("./coilPsi",-182.0));
	
	/* Set fixed phi angle for proline: */
	proteinCreatorProlinePhi=Math::rad(configFileSection.retrieveValue<double>("./prolinePhi",-60.0));
	
	/* Read setting for creation of terminator residues: */
	proteinCreatorCreateTerminatorResidues=configFileSection.retrieveValue<bool>("./createTerminatorResidues",true);
	}

void setCreateTerminatorResidues(bool newCreateTerminatorResidues)
	{
	proteinCreatorCreateTerminatorResidues=newCreateTerminatorResidues;
	}

Protein* loadProtein(const char* inputFileName)
	{
	/* Find input file extension: */
	const char* extension=0;
	for(const char* cPtr=inputFileName;*cPtr!='\0';++cPtr)
		if(*cPtr=='.')
			extension=cPtr+1;
	
	/* Determine input file type based on extension: */
	if(extension!=0&&strcasecmp(extension,"pdb")==0)
		{
		/* Load a PDB file and return only the first polypeptide chain: */
		std::vector<Protein*> chains=parsePdbFile(inputFileName);
		while(chains.size()>1)
			{
			delete chains.back();
			chains.pop_back();
			}
		return chains.front();
		}
	else if(extension!=0&&strcasecmp(extension,"pred")==0)
		{
		/* Create protein from prediction file: */
		return ReadPredictionFile(inputFileName);
		}
	else
		throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Unknown file extension in input file %s",inputFileName);
	}

}
