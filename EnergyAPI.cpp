/***********************************************************************
EnergyAPI - Classes to dynamically load energy computation libraries at
run-time.
Copyright (2) 2003 Oliver Kreylos
***********************************************************************/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dlfcn.h>
#include <string>
#include <stdexcept>

#include "Protein.h"
#include "ParsePdbFile.h"

#include "EnergyAPI.h"

/******************************
Methods of class EnergyLibrary:
******************************/

EnergyLibrary::EnergyLibrary(void* sDsoHandle)
	:dsoHandle(sDsoHandle)
	{
	}

EnergyLibrary* EnergyLibrary::load(const char* dsoName)
	{
	/* Open DSO containing class: */
	void* dsoHandle=dlopen(dsoName,RTLD_LAZY|RTLD_GLOBAL);
	if(dsoHandle==0) // Check for errors during DSO loading
		throw std::runtime_error(std::string("EnergyLibrary::load: Unable to open energy calculation DSO \"")+std::string(dsoName)+std::string("\" due to ")+std::string(dlerror()));
	
	/* Get address to library object creation function: */
	CreateEnergyLibraryFunction createEnergyLibrary=(CreateEnergyLibraryFunction)dlsym(dsoHandle,"createEnergyLibrary");
	if(createEnergyLibrary==0) // Check for errors during function lookup
		{
		dlclose(dsoHandle);
		throw std::runtime_error(std::string("EnergyLibrary::load: Unable to retrieve energy calculator creation function from DSO due to ")+std::string(dlerror()));
		}
	
	/* Create and return library object: */
	return createEnergyLibrary(dsoHandle);
	}

EnergyLibrary::~EnergyLibrary(void)
	{
	dlclose(dsoHandle);
	}

/*********************************
Methods of class EnergyCalculator:
*********************************/

int EnergyCalculator::getNumAtoms(void) const
	{
	return protein->getNumAtoms();
	}

int EnergyCalculator::getNumResidues(void) const
	{
	return protein->getNumResidues();
	}

void EnergyCalculator::getAtomCoordinates(float x[],float y[],float z[]) const
	{
	/* Copy atom positions: */
	MD::Protein::ConstAtomIterator aIt;
	int i;
	for(aIt=protein->atomsBegin(),i=0;aIt!=protein->atomsEnd();++aIt,++i)
		{
		x[i]=float(aIt->getPosition()[0]);
		y[i]=float(aIt->getPosition()[1]);
		z[i]=float(aIt->getPosition()[2]);
		}
	}

void EnergyCalculator::getAtomCoordinates(double x[],double y[],double z[]) const
	{
	/* Copy atom positions: */
	MD::Protein::ConstAtomIterator aIt;
	int i;
	for(aIt=protein->atomsBegin(),i=0;aIt!=protein->atomsEnd();++aIt,++i)
		{
		x[i]=aIt->getPosition()[0];
		y[i]=aIt->getPosition()[1];
		z[i]=aIt->getPosition()[2];
		}
	}

void EnergyCalculator::getAtomCoordinates(float xyz[][3]) const
	{
	/* Copy atom positions: */
	MD::Protein::ConstAtomIterator aIt;
	int i;
	for(aIt=protein->atomsBegin(),i=0;aIt!=protein->atomsEnd();++aIt,++i)
		for(int j=0;j<3;++j)
			xyz[i][j]=float(aIt->getPosition()[j]);
	}

void EnergyCalculator::getAtomCoordinates(double xyz[][3]) const
	{
	/* Copy atom positions: */
	MD::Protein::ConstAtomIterator aIt;
	int i;
	for(aIt=protein->atomsBegin(),i=0;aIt!=protein->atomsEnd();++aIt,++i)
		for(int j=0;j<3;++j)
			xyz[i][j]=aIt->getPosition()[j];
	}

EnergyCalculator::EnergyCalculator(const MD::Protein* sProtein)
	:protein(sProtein),pdbFileName(new char[256])
	{
	/* Create temporary PDB file to initialize energy computation state: */
	strcpy(pdbFileName,"EFuncPDB.pdb.XXXXXX");
	mktemp(pdbFileName);
	MD::writePdbFile(*protein,pdbFileName,false);
	}

EnergyCalculator::~EnergyCalculator(void)
	{
	/* Remove temporary PDB file: */
	unlink(pdbFileName);
	delete[] pdbFileName;
	}
