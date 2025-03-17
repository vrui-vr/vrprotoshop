/***********************************************************************
ParsePdbFile - Functions to read protein structures from a protein
database file.
Copyright (c) 2001-2020 Oliver Kreylos
***********************************************************************/

#ifndef PARSEPDBFILE_INCLUDED
#define PARSEPDBFILE_INCLUDED

#include <vector>

namespace MD {

/* Forward declarations: */
class Protein;

std::vector<Protein*> parsePdbFile(const char* filename);
bool writePdbFile(const Protein& protein,const char* filename,bool writeStructure =true);
bool updatePdbFilePositions(const Protein& protein,const char* templateFilename,const char* outputFilename);

}

#endif
