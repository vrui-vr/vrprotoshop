/***********************************************************************
ParsePdbFile - Functions to read protein structures from a protein
database file.
Copyright (c) 2001-2024 Oliver Kreylos
***********************************************************************/

#include "ParsePdbFile.h"

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <Misc/StdError.h>
#include <Misc/StandardHashFunction.h>
#include <Misc/HashTable.h>

#include "Protein.h"

namespace MD {

std::vector<Protein*> parsePdbFile(const char* filename)
	{
	typedef Misc::HashTable<char,std::string> StructureMap; // Hash table mapping polypeptide chain IDs to secondary structure type strings
	
	/* Open the input file and check for validity: */
	FILE* file=fopen(filename,"rt");
	if(file==0)
		throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Cannot open input file %s",filename);
	
	/* Read the input file one polypeptide chain at a time: */
	std::vector<Protein*> result;
	char currentChain='\0';
	Protein* protein=0;
	Protein::ResidueCreator* residueCreator=0;
	int currentResidueIndex=-1;
	char currentStructureType='\0';
	char nextStructureType='\0';
	
	/* Keep track of per-chain structure maps to parse new-style PDB files: */
	StructureMap structureMap(5);
	
	bool pdbFileHasStructureMaps=false;
	char line[1024];
	int lineNumber=1;
	while(fgets(line,sizeof(line),file)!=0)
		{
		/* Parse the current line: */
		char keyword[7];
		sscanf(&line[ 0],"%6s",keyword);
		if(strcmp(keyword,"REMARK")==0)
			{
			/* Check if this is a secondary structure type notification: */
			char secondaryStructureType[20];
			secondaryStructureType[0]='\0';
			sscanf(&line[12],"%20s",secondaryStructureType);
			if(strcmp(secondaryStructureType,"ALPHA-HELIX")==0)
				nextStructureType='A';
			else if(strcmp(secondaryStructureType,"BETA-STRAND")==0)
				nextStructureType='B';
			else if(strcmp(secondaryStructureType,"COIL")==0)
				nextStructureType='C';
			}
		else if(strcmp(keyword,"HELIX")==0)
			{
			/* Parse the helix definition: */
			char beginChain=line[19];
			int beginResidueIndex;
			sscanf(&line[20],"%d",&beginResidueIndex);
			char endChain=line[31];
			int endResidueIndex;
			sscanf(&line[32],"%d",&endResidueIndex);
			
			/* Check for consistency: */
			if(beginChain!=endChain)
				{
				/* Clean up and throw an exception: */
				delete residueCreator;
				for(std::vector<Protein*>::iterator rIt=result.begin();rIt!=result.end();++rIt)
					delete *rIt;
				throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Invalid HELIX definition in line %d of input file %s",lineNumber,filename);
				}
			
			/* Update the chain's secondary structure string: */
			std::string& structure=structureMap[beginChain].getDest();
			while(int(structure.length())<=endResidueIndex)
				structure.push_back('C');
			for(int i=beginResidueIndex;i<=endResidueIndex;++i)
				structure[i]='A';
			pdbFileHasStructureMaps=true;
			}
		else if(strcmp(keyword,"SHEET")==0)
			{
			/* Parse the sheet definition: */
			char beginChain=line[21];
			int beginResidueIndex;
			sscanf(&line[22],"%d",&beginResidueIndex);
			char endChain=line[32];
			int endResidueIndex;
			sscanf(&line[33],"%d",&endResidueIndex);
			
			/* Check for consistency: */
			if(beginChain!=endChain)
				{
				/* Clean up and throw an exception: */
				delete residueCreator;
				for(std::vector<Protein*>::iterator rIt=result.begin();rIt!=result.end();++rIt)
					delete *rIt;
				throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Invalid SHEET definition in line %d of input file %s",lineNumber,filename);
				}
			
			/* Update the chain's secondary structure string: */
			std::string& structure=structureMap[beginChain].getDest();
			while(int(structure.length())<=endResidueIndex)
				structure.push_back('C');
			for(int i=beginResidueIndex;i<=endResidueIndex;++i)
				structure[i]='B';
			pdbFileHasStructureMaps=true;
			}
		else if(strcmp(keyword,"ATOM")==0)
			{
			/* Parse the atom definition: */
			int atomIndex;
			sscanf(&line[ 6],"%d",&atomIndex);
			char atomName[5];
			sscanf(&line[12],"%4s",atomName);
			char residueName[4];
			sscanf(&line[17],"%3s",residueName);
			char chainId=line[21];
			int residueIndex;
			sscanf(&line[22],"%d",&residueIndex);
			double atomPosition[3];
			sscanf(&line[30],"%lg %lg %lg",&atomPosition[0],&atomPosition[1],&atomPosition[2]);
			
			/* Check if this atom starts a new chain: */
			if(currentChain!=chainId)
				{
				/* Finish constructing the previous chain: */
				if(residueCreator!=0)
					{
					residueCreator->finishProtein();
					delete residueCreator;
					}
				
				/* Create a new chain and residue creator: */
				protein=new Protein;
				result.push_back(protein);
				residueCreator=new Protein::ResidueCreator(protein);
				
				/* Reset state trackers: */
				currentResidueIndex=-1;
				currentStructureType='\0';
				
				currentChain=chainId;
				}
			
			/* Check if this atom starts a new residue: */
			if(currentResidueIndex!=residueIndex)
				{
				/* Check if this residue starts a new secondary structure: */
				if(pdbFileHasStructureMaps)
					{
					std::string& structures=structureMap.getEntry(currentChain).getDest();
					if(residueIndex<int(structures.length()))
						nextStructureType=structures[residueIndex];
					else
						nextStructureType='C';
					}
				if(currentStructureType!=nextStructureType)
					{
					switch(nextStructureType)
						{
						case 'A':
							residueCreator->newSecondaryStructure(Protein::SecondaryStructure::ALPHA_HELIX);
							break;
						
						case 'B':
							residueCreator->newSecondaryStructure(Protein::SecondaryStructure::BETA_STRAND);
							break;
						
						case 'C':
							residueCreator->newSecondaryStructure(Protein::SecondaryStructure::COIL);
							break;
						}
					
					currentStructureType=nextStructureType;
					}
				
				/* Add the new residue to the protein: */
				residueCreator->newResidue(residueName,residueIndex);
				
				currentResidueIndex=residueIndex;
				}
			
			/* Parse the element name and structure position: */
			char elementName[2];
			char* atomNamePtr=atomName;
			while(isspace(*atomNamePtr))
				++atomNamePtr;
			elementName[0]=*atomNamePtr;
			++atomNamePtr;
			elementName[1]='\0';

			/* Add the atom to the residue: */
			residueCreator->addAtom(elementName,atomIndex,Position(atomPosition),atomNamePtr);
			}
		
		++lineNumber;
		}
	
	/* Finish constructing the final chain: */
	if(residueCreator!=0)
		{
		residueCreator->finishProtein();
		delete residueCreator;
		}
	
	/* Clean up and return the constructed polypeptide chains: */
	fclose(file);
	return result;
	}

bool writePdbFile(const Protein& protein,const char* filename,bool writeStructure)
	{
	/* Open the input file and check for validity: */
	FILE* file=fopen(filename,"wt");
	if(file==0)
		return false;
	
	#if 1
	/* Write bogus energy value into file if no structure is requested: */
	if(!writeStructure)
		fprintf(file,"%f\n",0.0);
	#endif
	
	/* Start checking for secondary structure boundaries: */
	const Protein::SecondaryStructure* currentSecondaryStructure=0;
	
	/* Write the protein's atom in ascending index order: */
	for(Protein::ConstAtomIterator aIt=protein.atomsBegin();aIt!=protein.atomsEnd();++aIt)
		{
		const Protein::Residue* rPtr=aIt->getResidue();
		if(rPtr->getSecondaryStructure()!=currentSecondaryStructure)
			{
			if(writeStructure)
				{
				fprintf(file,"REMARK      ");
				switch(rPtr->getSecondaryStructure()->getStructureType())
					{
					case Protein::SecondaryStructure::ALPHA_HELIX:
						fprintf(file,"ALPHA-HELIX\n");
						break;

					case Protein::SecondaryStructure::BETA_STRAND:
						fprintf(file,"BETA-STRAND\n");
						break;

					case Protein::SecondaryStructure::COIL:
						fprintf(file,"COIL\n");
						break;
					
					default:
						; // Just to make compiler happy
					}
				}
			currentSecondaryStructure=rPtr->getSecondaryStructure();
			}
		char atomName[5];
		strcpy(atomName,aIt->getElementName());
		strcat(atomName,aIt->getPlacement());
		Position p=aIt->getPosition();
		fprintf(file,"ATOM   %4d %4s %3s %1s %3d     ",aIt->getAtomIndex(),atomName,rPtr->getPdbResidueName(),"",rPtr->getPdbResidueIndex()+1); // Hack to get around 1-based indices
		fprintf(file,"%7.3f %7.3f %7.3f %7.2f %4.2f\n",double(p[0]),double(p[1]),double(p[2]),-99.0,0.0);
		}
	fclose(file);
	
	return true;
	}

bool updatePdbFilePositions(const Protein& protein,const char* templateFilename,const char* outputFilename)
	{
	/* Open the input file and check for validity: */
	FILE* inputFile=fopen(templateFilename,"rt");
	if(inputFile==0)
		return false;
	
	/* Create a temporary output file: */
	FILE* outputFile=fopen(outputFilename,"wt");
	if(outputFile==0)
		{
		fclose(inputFile);
		return false;
		}
	
	/* Calculate the protein's centroid to center it in the output file: */
	MD::Point centroid=protein.calcCentroid();
	
	/* Process all lines from the input file, and replace atom positions with current values: */
	bool writingSuccessful=true;
	Protein::ConstAtomIterator aIt=protein.atomsBegin();
	while(!feof(inputFile))
		{
		char line[256];
		if(fgets(line,sizeof(line),inputFile)==0)
			break;
		if(strncmp(line,"ATOM ",5)==0||strncmp(line,"atom ",5)==0)
			{
			/* Get the atom's PDB index: */
			int atomIndex;
			sscanf(&line[ 6],"%d",&atomIndex);
			
			/* Find the same atom in the protein model, and write its position into the input line: */
			for(;aIt!=protein.atomsEnd()&&aIt->getAtomIndex()!=atomIndex;++aIt)
				;
			if(aIt!=protein.atomsEnd())
				{
				Vector pos=aIt->getPosition()-centroid;
				sprintf(&line[30],"%8.3f%8.3f%8.3f",double(pos[0]),double(pos[1]),double(pos[2]));
				line[54]=' ';
				}
			else
				writingSuccessful=false;
			}
		
		/* Write the (changed) input line to the output file: */
		fputs(line,outputFile);
		}
	
	/* Clean up and return: */
	fclose(inputFile);
	fclose(outputFile);
	return writingSuccessful;
	}

}
