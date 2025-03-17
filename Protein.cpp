/***********************************************************************
Protein - Class to represent a protein as a single chain of amino acid
residues.
Copyright (c) 2001-2021 Oliver Kreylos

This file is part of VR ProtoShop.

VR ProtoShop is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

VR ProtoShop is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with VR ProtoShop; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
***********************************************************************/

#include "Protein.h"

#include <assert.h>
#include <string.h>
#include <deque>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/AffineCombiner.h>
#include <Geometry/Matrix.h>
#include <Geometry/AffineTransformation.h>
#include <Geometry/Cylinder.h>

namespace MD {

/***********************************
Methods of class Protein::ChainAtom:
***********************************/

Protein::ChainAtom::ChainAtom(Atom::Element sType,const Position& sPosition,int sAtomIndex,const char* sPlacement,Protein::Residue* sResidue,int sBackboneIndex)
	:Atom(sType,sPosition),atomIndex(sAtomIndex),residue(sResidue),backboneIndex(sBackboneIndex),
	 pred(0),succ(0),rotationMarker(0)
	{
	strncpy(placement,sPlacement,3);
	placement[3]='\0';
	}

int Protein::ChainAtom::getBackboneDistance(void) const
	{
	/* Put atom onto a queue: */
	std::deque<std::pair<const ChainAtom*,int> > atomQueue;
	atomQueue.push_back(std::pair<const ChainAtom*,int>(this,0));
	
	/* Traverse the queue until the first item is a backbone atom: */
	while(atomQueue.front().first->backboneIndex<0)
		{
		/* Push the next atom's neighbors onto the queue: */
		const ChainAtom* a=atomQueue.front().first;
		int distance=atomQueue.front().second;
		atomQueue.pop_front();
		for(std::vector<Atom*>::const_iterator bIt=a->bonds.begin();bIt!=a->bonds.end();++bIt)
			atomQueue.push_back(std::pair<const ChainAtom*,int>(static_cast<const ChainAtom*>(*bIt),distance+1));
		}
	
	/* Return the first item's distance: */
	return atomQueue.front().second;
	}

/*****************************************
Static elements of class Protein::Residue:
*****************************************/

const char* const Protein::Residue::commonNames[23]={
	"Alanine","Arginine","Asparagine","Aspartic acid","ASP/ASN ambiguous","Cysteine",
	"Glutamine","Glutamic acid","GLU/GLN ambiguous","Glycine","Histidine","Isoleucine",
	"Leucine","Lysine","Methionine","Phenylalanine","Proline","Serine","Threonine",
	"Tryptophan","Tyrosine","Unknown","Valine"};

const char Protein::Residue::abbreviatedNames[23][4]={
	"ALA","ARG","ASN","ASP","ASX","CYS","GLN","GLU","GLX","GLY",
	"HID","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","UNK","VAL"};

const char Protein::Residue::abbreviatedLcNames[23][4]={
	"Ala","Arg","Asn","Asp","Asx","Cys","Gln","Glu","Glx","Gly",
	"His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Unk","Val"};

const char Protein::Residue::singleLetterNames[23]={
	'A','R','N','D','X','C','Q','E','Z','G',
	'H','I','L','K','M','F','P','S','T','W','Y','U','V'};

/*********************************
Methods of class Protein::Residue:
*********************************/

Protein::Residue::Residue(const char* sResidueName,int sResidueIndex,Protein::SecondaryStructure* sSecondaryStructure)
	:residueIndex(sResidueIndex),beginAtom(0),endAtom(0),secondaryStructure(sSecondaryStructure),
	 pred(0),succ(0)
	{
	/* Keep the residue name around (it might be one we don't know): */
	strncpy(residueName,sResidueName,3);
	residueName[3]='\0';
	
	/* Try to parse the residue name: */
	type=parseType(residueName);
	}

Protein::Residue::AminoAcid Protein::Residue::parseType(const char* abbreviatedName)
	{
	/* Search through the list of abbreviated names to find the given one: */
	int index;
	for(index=0;index<23;++index)
		if(strcmp(abbreviatedNames[index],abbreviatedName)==0)
			break;
	
	/* If nothing matched, punt: */
	if(index>=23)
		return UNK;
	else return AminoAcid(index);
	}

int Protein::Residue::getNumAtoms(void) const
	{
	int result=0;
	for(const ChainAtom* aPtr=beginAtom;aPtr!=endAtom;aPtr=aPtr->succ)
		++result;
	return result;
	}

Protein::Dipole Protein::Residue::getAmide(void) const
	{
	/* Assume that the residue type is neither unknown nor proline: */
	if(type!=UNK&&type!=PRO)
		{
		ChainAtom* major=backbone[0];
		ChainAtom* minor=0;
		
		/* Find hydrogen atom: */
		for(std::vector<Atom*>::const_iterator bIt=major->getBonds().begin();bIt!=major->getBonds().end();++bIt)
			if((*bIt)->getType()==Atom::H)
				{
				minor=static_cast<ChainAtom*>(*bIt);
				break;
				}
		
		return Dipole(major,minor,1);
		}
	else
		return Dipole(); // Return invalid dipole
	}

Protein::Dipole Protein::Residue::getCarboxyl(void) const
	{
	/* Assume that the residue type is not unknown: */
	if(type!=UNK)
		{
		ChainAtom* major=backbone[2];
		ChainAtom* minor=0;
		
		/* Find oxygen atom: */
		for(std::vector<Atom*>::const_iterator bIt=major->getBonds().begin();bIt!=major->getBonds().end();++bIt)
			if((*bIt)->getType()==Atom::O)
				{
				minor=static_cast<ChainAtom*>(*bIt);
				break;
				}
		
		return Dipole(major,minor,-1);
		}
	else
		return Dipole(); // Return invalid dipole
	}

Point Protein::Residue::getSideChainPosition(void) const
	{
	/* Process residue based on type: */
	if(type==UNK)
		{
		/* Just return the centroid of the entire residue: */
		Point::AffineCombiner centroidCombiner;
		for(const ChainAtom* aPtr=beginAtom;aPtr!=endAtom;aPtr=aPtr->succ)
			centroidCombiner.addPoint(aPtr->position);
		return centroidCombiner.getPoint();
		}
	else
		{
		/* Return position of carbon atom next to carbon alpha (or position of carbon alpha in case of glycine): */
		const ChainAtom* ca=backbone[1];
		Point result=ca->getPosition();
		for(std::vector<Atom*>::const_iterator bIt=ca->bonds.begin();bIt!=ca->bonds.end();++bIt)
			if((*bIt)->getType()==Atom::C&&static_cast<const ChainAtom*>(*bIt)->getBackboneIndex()<0)
				{
				result=(*bIt)->getPosition();
				break;
				}
		return result;
		}
	}

Protein::DihedralAnglePair Protein::Residue::calcDihedralAngles(void) const
	{
	DihedralAnglePair result(0,0);
	
	/* Don't calculate dihedral angles of unknown residue types: */
	if(type!=UNK)
		{
		/* Get pointers to the five atoms defining the residue's three planes: */
		const ChainAtom* atom0=pred!=0?pred->backbone[pred->backbone.size()-1]:backbone[0];
		const ChainAtom* atom1=backbone[0];
		const ChainAtom* atom2=backbone[1];
		const ChainAtom* atom3=backbone[2];
		const ChainAtom* atom4=succ!=0?succ->backbone[0]:backbone[2];

		/* Calculate the planes' normal vectors: */
		Vector d01=atom1->getPosition()-atom0->getPosition();
		Vector d12=atom2->getPosition()-atom1->getPosition();
		Vector d23=atom3->getPosition()-atom2->getPosition();
		Vector d34=atom4->getPosition()-atom3->getPosition();
		Vector n1=Geometry::cross(d01,d12);
		double n1Len=Geometry::mag(n1);
		Vector n2=Geometry::cross(d12,d23);
		double n2Len=Geometry::mag(n2);
		Vector n3=Geometry::cross(d23,d34);
		double n3Len=Geometry::mag(n3);
		
		/* Don't calculate the first dihedral angle of the first residue in a protein: */
		if(pred!=0)
			{
			/* Phi is the angle between n1 and n2: */
			result.phi=Math::acos((n1*n2)/Scalar(n1Len*n2Len));
			if(Geometry::cross(n1,n2)*d12<Scalar(0)) // Angle is positive if clockwise along main bond
				result.phi=-result.phi;
			}

		/* Don't calculate the second dihedral angle of the last residue in a protein: */
		if(succ!=0)
			{
			/* Psi is the angle between n2 and n3: */
			result.psi=Math::acos((n2*n3)/Scalar(n2Len*n3Len));
			if(Geometry::cross(n2,n3)*d23<Scalar(0)) // Angle is positive if clockwise along main bond
				result.psi=-result.psi;
			}
		}
	
	return result;
	}

/********************************************
Methods of class Protein::SecondaryStructure:
********************************************/

int Protein::SecondaryStructure::getFirstResidueIndex(void) const
	{
	int result=-1;
	for(const Residue* rPtr=residueBegin;rPtr!=0;rPtr=rPtr->pred)
		++result;
	return result;
	}

int Protein::SecondaryStructure::getNumResidues(void) const
	{
	int result=0;
	for(const Residue* rPtr=residueBegin;rPtr!=residueEnd;rPtr=rPtr->succ)
		++result;
	return result;
	}

/****************************************
Static elements of class Protein::Dipole:
****************************************/

const DistanceRange Protein::Dipole::hydrogenBondRangeMajor(2.7,3.2);
const DistanceRange Protein::Dipole::hydrogenBondRangeMinor(1.7,2.35);
const Scalar Protein::Dipole::NHdist=Scalar(1.006);
const Scalar Protein::Dipole::cosAngleMin=Scalar(0.707);

/********************************
Methods of class Protein::Dipole:
********************************/

Point Protein::Dipole::getBondSite(void) const
	{
	/* Calculate a direction vector from the major to the minor atom: */
	Point minorPos=minorAtom->getPosition();
	Vector dir=minorPos-majorAtom->getPosition();
	
	/* Offset the minor atom position by half the hydrogen bond distance: */
	double averageBondDistance=(hydrogenBondRangeMinor.getMin()+hydrogenBondRangeMinor.getMax())*0.5;
	return minorPos+dir*Scalar(averageBondDistance*0.5/Geometry::mag(dir));
	}

bool formHydrogenBond(const Protein::Dipole& amide,const Protein::Dipole& carboxyl)
	{
	/**************************************************************
	This function silently assumes that the first given dipole is
	an amide group and the second given dipole is a carboxyl group.
	**************************************************************/
	
	/* Check for distance requirements: */
	if(Protein::Dipole::hydrogenBondRangeMajor.areInRange(amide.majorAtom->getPosition(),carboxyl.minorAtom->getPosition())&&
		 Protein::Dipole::hydrogenBondRangeMinor.areInRange(amide.minorAtom->getPosition(),carboxyl.minorAtom->getPosition()))
		{
		/* Check for alignment requirements: */
		Vector A=amide.minorAtom->getPosition()-amide.majorAtom->getPosition();
		Vector B=carboxyl.minorAtom->getPosition()-amide.minorAtom->getPosition();
		Scalar cosAngle=(A*B)/Math::sqrt(Geometry::sqr(A)*Geometry::sqr(B));
		return cosAngle>=Protein::Dipole::cosAngleMin;
		}
	else
		return false;
	}

/************************************************
Static elements of class Protein::ResidueCreator:
************************************************/

const double Protein::ResidueCreator::covalentBondThreshold=1.85;

/****************************************
Methods of class Protein::ResidueCreator:
****************************************/

Protein::ResidueCreator::ResidueCreator(Protein* sProtein)
	:protein(sProtein),lastAtom(0),lastBackboneAtom(0),currentSecondaryStructure(0),currentResidue(0),
	 firstAtom(false),buildingBackbone(false),firstResidue(false)
	{
	}

Protein::ResidueCreator::~ResidueCreator(void)
	{
	}

void Protein::ResidueCreator::newSecondaryStructure(Protein::SecondaryStructure::StructureType sStructureType)
	{
	/* Create a new secondary structure: */
	SecondaryStructure* newSecondaryStructure=new SecondaryStructure(sStructureType);
	newSecondaryStructure->pred=currentSecondaryStructure;
	if(currentSecondaryStructure!=0)
		currentSecondaryStructure->succ=newSecondaryStructure;
	else
		protein->secondaryStructures=newSecondaryStructure;
	currentSecondaryStructure=newSecondaryStructure;
	
	/* Reset addition flags: */
	firstResidue=true;
	}

void Protein::ResidueCreator::newResidue(const char* newAbbreviatedName,int newResidueIndex)
	{
	/* Add a new residue to the chain of residues: */
	Residue* newResidue=new Residue(newAbbreviatedName,newResidueIndex,currentSecondaryStructure);
	newResidue->pred=currentResidue;
	if(currentResidue!=0)
		currentResidue->succ=newResidue;
	else
		protein->residues=newResidue;
	currentResidue=newResidue;
	
	/* Reset addition flags: */
	firstAtom=true;
	buildingBackbone=true;
	
	/* Start adding residues to a new secondary structure if necessary: */
	if(firstResidue)
		{
		if(currentSecondaryStructure->pred!=0)
			currentSecondaryStructure->pred->residueEnd=currentResidue;
		currentSecondaryStructure->residueBegin=currentResidue;
		firstResidue=false;
		}
	}

void Protein::ResidueCreator::addAtom(const char* elementName,int atomIndex,const Position& atomPosition,const char* placement)
	{
	/* Add the new atom to the atom list: */
	++protein->numAtoms;
	ChainAtom* newAtom=new ChainAtom(Atom::parseType(elementName),atomPosition,atomIndex,placement,currentResidue);
	newAtom->pred=lastAtom;
	if(lastAtom!=0)
		lastAtom->succ=newAtom;
	else
		protein->atoms=newAtom;
	lastAtom=newAtom;
	
	/* Add the new atom to the current residue: */
	if(firstAtom)
		{
		if(currentResidue->pred!=0)
			currentResidue->pred->endAtom=newAtom;
		currentResidue->beginAtom=newAtom;
		firstAtom=false;
		}
	
	/* Add the new atom to the backbone if appropriate: */
	if(buildingBackbone&&(newAtom->getType()==Atom::N||newAtom->getType()==Atom::C))
		{
		/* The backbone is finished once the first non-backbone carbon is encountered: */
		if(newAtom->getType()==Atom::C&&placement[0]!='\0'&&placement[0]!='A'&&placement[0]!='Y')
			buildingBackbone=false;
		else
			{
			newAtom->backboneIndex=currentResidue->backbone.size();
			if(lastBackboneAtom!=0)
				bond(*newAtom,*lastBackboneAtom);
			currentResidue->backbone.push_back(newAtom);
			lastBackboneAtom=newAtom;
			}
		}
	
	/* Bond the new atom with all appropriate existing ones: */
	for(ChainAtom* bondSearch=newAtom;bondSearch!=currentResidue->beginAtom;bondSearch=bondSearch->pred)
		{
		ChainAtom* testAtom=bondSearch->pred;
		
		/* Skip all existing hydrogen atoms - they must already be bonded: */
		if(testAtom->getType()==Atom::H)
			continue;
		
		/* Bond two atoms if their distance is less than the covalent bond threshold: */
		if(dist(*newAtom,*testAtom)<covalentBondThreshold)
			{
			/* Create backbone amide/carboxyl groups on-the-fly: */
			if(newAtom->getType()==Atom::H&&testAtom->getType()==Atom::N&&testAtom->getBackboneIndex()==0)
				{
				/* Build an amide group: */
				protein->amides.push_back(Dipole(testAtom,newAtom,1));
				}
			else if(newAtom->getType()==Atom::O&&testAtom->getType()==Atom::C&&testAtom->getBackboneIndex()==2)
				{
				/* Build a carboxyl group: */
				protein->carboxyls.push_back(Dipole(testAtom,newAtom,-1));
				}
			
			/* Bond the atoms: */
			bond(*newAtom,*testAtom);
			if(newAtom->getType()==Atom::H)
				break; // Only one bond per H atom, and to the closest "heavy" atom back in the chain
			}
		}
	}

void Protein::ResidueCreator::finishProtein(void)
	{
	/* Finish the last secondary structure: */
	if(currentSecondaryStructure!=0)
		currentSecondaryStructure->residueEnd=0;
	
	/* Associate dipoles with secondary structures: */
	SecondaryStructure* sPtr=0;
	for(std::vector<Dipole>::iterator nh=protein->amides.begin();nh!=protein->amides.end();++nh)
		{
		if(nh->getMajorAtom()->residue->secondaryStructure!=sPtr)
			{
			if(sPtr!=0)
				sPtr->amidesEnd=nh;
			sPtr=nh->getMajorAtom()->residue->secondaryStructure;
			sPtr->amidesBegin=nh;
			}
		}
	if(sPtr!=0)
		sPtr->amidesEnd=protein->amides.end();
	sPtr=0;
	for(std::vector<Dipole>::iterator co=protein->carboxyls.begin();co!=protein->carboxyls.end();++co)
		{
		if(co->getMajorAtom()->residue->secondaryStructure!=sPtr)
			{
			if(sPtr!=0)
				sPtr->carboxylsEnd=co;
			sPtr=co->getMajorAtom()->residue->secondaryStructure;
			sPtr->carboxylsBegin=co;
			}
		}
	if(sPtr!=0)
		sPtr->carboxylsEnd=protein->carboxyls.end();
	
	/* Create backbone bond iterators: */
	Residue* rPtr=protein->residues;
	std::vector<ChainAtom*>::iterator atomIt1,atomIt2;
	atomIt1=atomIt2=rPtr->backbone.begin();
	++atomIt2;
	protein->beginIterator=BackboneIterator(atomIt1,atomIt2);
	while(rPtr->succ!=0)
		rPtr=rPtr->succ;
	atomIt1=atomIt2=rPtr->backbone.end();
	--atomIt1;
	protein->endIterator=BackboneIterator(atomIt1,atomIt2);
	
	/* Create array of pointers to atoms: */
	protein->atomPointers=new ChainAtom*[protein->numAtoms];
	int index=0;
	for(ChainAtom* aPtr=protein->atoms;aPtr!=0;aPtr=aPtr->succ,++index)
		protein->atomPointers[index]=aPtr;
	}

/******************************************
Methods of class Protein::BackboneIterator:
******************************************/

int Protein::BackboneIterator::getResidueIndex(void) const
	{
	int result=-1;
	for(const Residue* rPtr=(*atom1)->residue;rPtr!=0;rPtr=rPtr->pred)
		++result;
	return result;
	}

/*******************************************
Methods of class Protein::StructureSelector:
*******************************************/

Protein::StructureSelector::StructureSelector(Protein* sProtein,SecondaryStructure* sStructure,int sStructureIndex)
	:protein(sProtein),structure(sStructure),structureIndex(sStructureIndex),
	 firstResidueIndex(-1),numResidues(0)
	{
	if(isValid())
		{
		/* Cache first residue index and number of residues: */
		firstResidueIndex=structure->getFirstResidueIndex();
		numResidues=structure->getNumResidues();
		
		/* Create iterators for beginning and end of backbone: */
		std::vector<ChainAtom*>::iterator atomIt1,atomIt2;
		
		/* Create the begin iterator out of the first residue's first two backbone atoms: */
		atomIt1=atomIt2=structure->residueBegin->backbone.begin();
		++atomIt2;
		begin=BackboneIterator(atomIt1,atomIt2);
		
		/* Find the last residue in the structure: */
		Residue* rPtr;
		for(rPtr=structure->residueBegin;rPtr->succ!=structure->residueEnd;rPtr=rPtr->succ)
			;
		
		/* Create an iterator out of the last residue's first two backbone atoms: */
		atomIt1=atomIt2=rPtr->backbone.begin();
		++atomIt2;
		end=BackboneIterator(atomIt1,atomIt2);
		
		/* Increment the end iterator until it points past the last residue's end: */
		for(size_t i=1;i<rPtr->backbone.size();++i)
			++end;
		}
	}

void Protein::StructureSelector::invalidate(void)
	{
	/* Reset the structure selector to invalid state: */
	structure=0;
	structureIndex=-1;
	firstResidueIndex=-1;
	numResidues=0;
	}

int Protein::StructureSelector::getStructureTypeIndex(void) const
	{
	int result=0;
	for(const SecondaryStructure* sPtr=protein->secondaryStructures;sPtr!=structure;sPtr=sPtr->succ)
		if(sPtr->getStructureType()==structure->getStructureType())
			++result;
	return result;
	}

/*************************************
Methods of class Protein::BondRotator:
*************************************/

void Protein::BondRotator::rotate(Scalar angle)
	{
	/* Calculate rotation transformations for both parts of the molecule: */
	Geometry::AffineTransformation<Scalar,3> rotation[2];
	for(int i=0;i<2;++i)
		{
		rotation[i]=Geometry::AffineTransformation<Scalar,3>::identity;
		rotation[i]*=Geometry::AffineTransformation<Scalar,3>::translate(getStart()-Point::origin);
		rotation[i]*=Geometry::AffineTransformation<Scalar,3>::rotate(Geometry::Rotation<Scalar,3>::rotateAxis(getAxis(),i==1?Scalar(0.5)*angle:-Scalar(0.5)*angle));
		rotation[i]*=Geometry::AffineTransformation<Scalar,3>::translate(Point::origin-getStart());
		}
	
	/* Perform a breadth-first traversal on all atoms: */
	++protein->rotationMarker;
	std::deque<MoleculeTraversal> traversalQueue;
	(*atom1)->rotationMarker=protein->rotationMarker;
	traversalQueue.push_back(MoleculeTraversal(*atom1,0));
	(*atom2)->rotationMarker=protein->rotationMarker;
	traversalQueue.push_back(MoleculeTraversal(*atom2,1));
	while(!traversalQueue.empty())
		{
		/* Pick the first element from the queue: */
		MoleculeTraversal mt=traversalQueue.front();
		traversalQueue.pop_front();
		
		/* Rotate it: */
		mt.atom->setPosition(rotation[mt.rotationIndex].transform(mt.atom->getPosition()));
		
		/* Put all its neighbours into the queue: */
		for(std::vector<Atom*>::iterator aIt=mt.atom->getBonds().begin();aIt!=mt.atom->getBonds().end();++aIt)
			{
			ChainAtom* aPtr=static_cast<ChainAtom*>(*aIt);
			if(aPtr->rotationMarker!=protein->rotationMarker)
				{
				aPtr->rotationMarker=protein->rotationMarker;
				traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex));
				}
			}
		}
	
	/* Increment version numbers of all affected secondary structures: */
	for(SecondaryStructure* sPtr=protein->secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
		++sPtr->version;
	}

Scalar Protein::BondRotator::getDihedralAngle(void) const
	{
	/* Find the atoms before atom1 and after atom2 in the backbone chain: */
	assert((*atom1)->residue==(*atom2)->residue); // Both atoms must be in the same residue (peptide bonds don't rotate)
	const Residue* res=(*atom1)->residue;
	std::vector<ChainAtom*>::const_iterator atomIt;
	
	/* Find the atom before atom1 by going backwards: */
	const ChainAtom* atom0;
	atomIt=atom1;
	if(atomIt!=res->backbone.begin())
		{
		--atomIt;
		atom0=*atomIt;
		}
	else if(res->pred!=0)
		atom0=res->pred->backbone.back();
	else
		return Scalar(0); // Ugly hack: The first phi angle is set to zero
	
	/* Find the atom behind atom2 by going forward: */
	const ChainAtom* atom3;
	atomIt=atom2;
	++atomIt;
	if(atomIt!=res->backbone.end())
		atom3=*atomIt;
	else if(res->succ!=0)
		atom3=res->succ->backbone.front();
	else
		return Scalar(0); // Ugly hack: The last psi angle is set to zero
	
	/* Now calculate the normal vectors for the planes defined by the two atom triples: */
	Vector d1=(*atom1)->getPosition()-atom0->getPosition();
	Vector d2=(*atom2)->getPosition()-(*atom1)->getPosition();
	Vector d3=atom3->getPosition()-(*atom2)->getPosition();
	Vector n1=Geometry::cross(d1,d2);
	Vector n2=Geometry::cross(d2,d3);
	
	/* The dihedral angle is the angle between the two normals: */
	Scalar angle=Math::acos((n1*n2)/(Geometry::mag(n1)*Geometry::mag(n2)));
	
	/* Determine the angle's sign by comparing the cross product of n1 and n2 with the bond axis: */
	if(Geometry::cross(n1,n2)*d2<Scalar(0))
		angle=-angle;
	
	return angle;
	}

/********************************
Static elements of class Protein:
********************************/

Scalar Protein::bondRadius=Scalar(0.25);

/************************
Methods of class Protein:
************************/

Protein::Protein(void)
	:numAtoms(0),atoms(0),residues(0),secondaryStructures(0),
	 secondaryStructureVersion(1),
	 rotationMarker(0),atomPointers(0)
	{
	}

Protein::~Protein(void)
	{
	/* Delete all secondary structures: */
	while(secondaryStructures!=0)
		{
		SecondaryStructure* succ=secondaryStructures->succ;
		delete secondaryStructures;
		secondaryStructures=succ;
		}
	
	/* Delete all residues: */
	while(residues!=0)
		{
		Residue* succ=residues->succ;
		delete residues;
		residues=succ;
		}
	
	/* Delete all atoms: */
	while(atoms!=0)
		{
		ChainAtom* succ=atoms->succ;
		delete atoms;
		atoms=succ;
		}
	
	/* Delete the atom pointers: */
	delete[] atomPointers;
	}

Position Protein::calcCentroid(void) const
	{
	/* Iterate through all atoms to calculate the protein's centroid: */
	Point::AffineCombiner centroidCombiner;
	for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
		centroidCombiner.addPoint(aPtr->getPosition());
	
	return centroidCombiner.getPoint();
	}

double Protein::calcRadius(void) const
	{
	/* Iterate through all atoms to calculate the protein's centroid: */
	Point::AffineCombiner centroidCombiner;
	for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
		centroidCombiner.addPoint(aPtr->getPosition());
	Point centroid=centroidCombiner.getPoint();
	
	/* Calculate the maximum distance of any atom from the centroid: */
	double radius2=0.0;
	for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
		{
		double d2=double(Geometry::sqrDist(Point(aPtr->getPosition()),centroid));
		if(radius2<d2)
			radius2=d2;
		}
	
	return Math::sqrt(radius2);
	}

int Protein::getNumResidues(void) const
	{
	int result=0;
	for(const Residue* rPtr=residues;rPtr!=0;rPtr=rPtr->succ)
		++result;
	return result;
	}

std::pair<int,int> Protein::getResidueIndexRange(void) const
	{
	const Residue* rPtr=residues;
	int rangeMin,rangeMax;
	rangeMin=rangeMax=rPtr->residueIndex;
	for(rPtr=rPtr->succ;rPtr!=0;rPtr=rPtr->succ)
		{
		if(rangeMin>rPtr->residueIndex)
			rangeMin=rPtr->residueIndex;
		else if(rangeMax<rPtr->residueIndex)
			rangeMax=rPtr->residueIndex;
		}
	
	return std::pair<int,int>(rangeMin,rangeMax);
	}

int Protein::getNumStructures(void) const
	{
	int result=0;
	for(const SecondaryStructure* sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
		++result;
	return result;
	}

Protein::ConstAtomIterator Protein::pickAtom(const Point& p) const
	{
	Scalar minDist2=Math::Constants<Scalar>::max;
	const ChainAtom* selected=0;
	
	/* Iterate through all atoms in the protein and check the query point against each: */
	for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
		{
		/* Check if the query point is inside the atom and closer than the previously selected atom: */
		Scalar dist2=Geometry::sqrDist(aPtr->getPosition(),p);
		if(dist2<minDist2&&dist2<=Math::sqr(aPtr->getVanDerWaalsRadius()))
			{
			selected=aPtr;
			minDist2=dist2;
			}
		}
	
	return ConstAtomIterator(selected);
	}

Protein::ConstAtomIterator Protein::pickAtom(const Ray& ray) const
	{
	Scalar lambdaMin=Math::Constants<Scalar>::max;
	const ChainAtom* selected=0;
	
	/* Iterate through all atoms in the protein and intersect the query ray with each: */
	for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
		{
		/* Intersect the query ray with the atom: */
		Scalar lambda=aPtr->intersectRay(ray);
		if(lambdaMin>lambda) // Intersection is the closest so far
			{
			selected=aPtr;
			lambdaMin=lambda;
			}
		}
	
	return ConstAtomIterator(selected);
	}

const Protein::Residue* Protein::pickResidue(int pdbResidueIndex) const
	{
	const Residue* rPtr;
	for(rPtr=residues;rPtr!=0&&rPtr->residueIndex!=pdbResidueIndex;rPtr=rPtr->succ)
		;
	return rPtr;
	}

Protein::Residue* Protein::pickResidue(int pdbResidueIndex)
	{
	Residue* rPtr;
	for(rPtr=residues;rPtr!=0&&rPtr->residueIndex!=pdbResidueIndex;rPtr=rPtr->succ)
		;
	return rPtr;
	}

Protein::Residue* Protein::pickResidue(const Point& p)
	{
	Residue* selected=0;
	Scalar distance2Min=Math::Constants<Scalar>::max;
	
	/* Iterate through all atoms in the protein and check each against the query point: */
	int i=0;
	for(SecondaryStructure* sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ,++i)
		{
		for(Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
			{
			for(ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
				{
				/* Intersect the query ray with the covalent bond sphere around the atom: */
				Scalar distance2=Geometry::sqrDist(aPtr->getPosition(),p);
				if(distance2<=Math::sqr(aPtr->getCovalentRadius()*Scalar(1.5))&&distance2<distance2Min)
					{
					distance2Min=distance2;
					selected=rPtr;
					}
				}
			}
		}
	
	return selected;
	}

Protein::Residue* Protein::pickResidue(const Ray& ray)
	{
	Scalar lambdaMin=Math::Constants<Scalar>::max;
	Residue* selected=0;
	
	/* Iterate through all backbone bonds in the protein and intersect the query ray with each: */
	for(BackboneIterator bbIt=beginIterator;bbIt!=endIterator;++bbIt)
		{
		/* Construct a cylinder around the bond and intersect the query ray with it: */
		Geometry::Cylinder<Scalar,3> bond(bbIt.getAtom1()->getPosition(),bbIt.getAtom2()->getPosition(),bondRadius);
		Scalar lambda=bond.intersectRay(ray).getParameter();
		if(lambdaMin>lambda)
			{
			lambdaMin=lambda;
			if(bbIt.isPeptideBond())
				{
				/* Check which side of the peptide bond was selected: */
				Point intersection=ray(lambda);
				if(Geometry::sqrDist(bbIt.getAtom1()->getPosition(),intersection)<=Geometry::sqrDist(bbIt.getAtom2()->getPosition(),intersection))
					selected=bbIt.getAtom1()->getResidue();
				else
					selected=bbIt.getAtom2()->getResidue();
				}
			else
				selected=bbIt.getResidue();
			}
		}
	
	return selected;
	}

int Protein::getResidueIndex(const Protein::Residue* rPtr) const
	{
	int result=0;
	const Residue* rPtr2;
	for(rPtr2=residues;rPtr2!=0&&rPtr2!=rPtr;rPtr2=rPtr2->succ,++result)
		;
	if(rPtr2==0)
		result=-1;
	return result;
	}

Point Protein::getResidueSideChainCentroid(const Protein::Residue* rPtr) const
	{
	/* Process residue based on type: */
	if(rPtr->type==Residue::UNK)
		{
		/* Just return the centroid of the entire residue: */
		Point::AffineCombiner centroidCombiner;
		for(const ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
			centroidCombiner.addPoint(aPtr->position);
		return centroidCombiner.getPoint();
		}
	
	/* Find the carbon atom next to carbon alpha: */
	const ChainAtom* ca=rPtr->backbone[1];
	const ChainAtom* cside=0;
	for(std::vector<Atom*>::const_iterator bIt=ca->bonds.begin();bIt!=ca->bonds.end();++bIt)
		if((*bIt)->getType()==Atom::C&&static_cast<const ChainAtom*>(*bIt)->getBackboneIndex()<0)
			{
			cside=static_cast<const ChainAtom*>(*bIt);
			break;
			}
	
	/* Check that a carbon was found (i.e., the residue type is not glycine): */
	if(cside==0)
		{
		/* Return the carbon alpha's position (glycine does not have a side chain): */
		return ca->position;
		}
	
	/* Compute the side chain's centroid: */
	Point::AffineCombiner centroidCombiner;
	
	/* Update the protein's traversal marker: */
	++const_cast<Protein*>(this)->rotationMarker; // Ugly; rotationMarker should be mutable
	
	/* Update the traversal markers on all backbone atoms to prevent traversing them: */
	for(std::vector<ChainAtom*>::const_iterator bbIt=rPtr->backbone.begin();bbIt!=rPtr->backbone.end();++bbIt)
		const_cast<ChainAtom*>(*bbIt)->rotationMarker=rotationMarker;
	
	/* Add the first side chain atom to the traversal queue: */
	std::deque<const ChainAtom*> queue;
	queue.push_back(cside);
	const_cast<ChainAtom*>(cside)->rotationMarker=rotationMarker;
	
	/* Traverse all side chain atoms exactly once: */
	while(!queue.empty())
		{
		/* Add the next atom from the queue to the centroid: */
		const ChainAtom* aPtr=queue.front();
		queue.pop_front();
		centroidCombiner.addPoint(aPtr->position);
		
		/* Add the next atom's neighbors to the queue: */
		for(std::vector<Atom*>::const_iterator bIt=aPtr->bonds.begin();bIt!=aPtr->bonds.end();++bIt)
			{
			/* Only add non-backbone atoms that have not been added yet: */
			const ChainAtom* nPtr=static_cast<const ChainAtom*>(*bIt);
			if(nPtr->rotationMarker!=rotationMarker)
				{
				queue.push_back(nPtr);
				const_cast<ChainAtom*>(nPtr)->rotationMarker=rotationMarker;
				}
			}
		}
	
	return centroidCombiner.getPoint();
	}

void Protein::changeResidueStructureType(Protein::Residue* rPtr,Protein::SecondaryStructure::StructureType newType)
	{
	SecondaryStructure* sPtr=rPtr->secondaryStructure;
	if(newType!=sPtr->structureType)
		{
		if(rPtr!=sPtr->residueBegin&&rPtr->succ!=sPtr->residueEnd) // Residue is in the middle of structure
			{
			/* Split current structure into three: */
			SecondaryStructure* leftStructure=new SecondaryStructure(sPtr->structureType);
			leftStructure->residueBegin=sPtr->residueBegin;
			leftStructure->residueEnd=sPtr->residueBegin=rPtr;
			SecondaryStructure* rightStructure=new SecondaryStructure(sPtr->structureType);
			rightStructure->residueEnd=sPtr->residueEnd;
			rightStructure->residueBegin=sPtr->residueEnd=rPtr->succ;
			sPtr->structureType=newType;
			leftStructure->pred=sPtr->pred;
			if(leftStructure->pred!=0)
				leftStructure->pred->succ=leftStructure;
			else
				secondaryStructures=leftStructure;
			leftStructure->succ=sPtr;
			sPtr->pred=leftStructure;
			rightStructure->succ=sPtr->succ;
			if(rightStructure->succ!=0)
				rightStructure->succ->pred=rightStructure;
			rightStructure->pred=sPtr;
			sPtr->succ=rightStructure;
			}
		else if(rPtr==sPtr->residueBegin&&rPtr->succ!=sPtr->residueEnd)
			{
			/* Change first residue in structure: */
			if(sPtr->pred==0||sPtr->pred->structureType!=newType)
				{
				/* Insert new structure before sPtr: */
				SecondaryStructure* leftStructure=new SecondaryStructure(newType);
				leftStructure->residueBegin=sPtr->residueBegin;
				leftStructure->residueEnd=sPtr->residueBegin=rPtr->succ;
				leftStructure->pred=sPtr->pred;
				if(sPtr->pred!=0)
					sPtr->pred->succ=leftStructure;
				else
					secondaryStructures=leftStructure;
				leftStructure->succ=sPtr;
				sPtr->pred=leftStructure;
				}
			else
				{
				/* Move separation between sPtr->pred and sPtr: */
				sPtr->pred->residueEnd=sPtr->residueBegin=rPtr->succ;
				}
			}
		else if(rPtr!=sPtr->residueBegin&&rPtr->succ==sPtr->residueEnd)
			{
			/* Change last residue in structure: */
			if(sPtr->succ==0||sPtr->succ->structureType!=newType)
				{
				/* Insert new structure after sPtr: */
				SecondaryStructure* rightStructure=new SecondaryStructure(newType);
				rightStructure->residueEnd=sPtr->residueEnd;
				rightStructure->residueBegin=sPtr->residueEnd=rPtr;
				rightStructure->succ=sPtr->succ;
				if(sPtr->succ!=0)
					sPtr->succ->pred=rightStructure;
				rightStructure->pred=sPtr;
				sPtr->succ=rightStructure;
				}
			else
				{
				/* Move separation between sPtr and sPtr->succ: */
				sPtr->residueEnd=sPtr->succ->residueBegin=rPtr;
				}
			}
		else
			{
			/* Change entire structure: */
			sPtr->structureType=newType;
			SecondaryStructure* left=sPtr->pred;
			if(left!=0&&left->structureType==newType)
				{
				/* Merge with structure to the left: */
				sPtr->residueBegin=left->residueBegin;
				sPtr->pred=left->pred;
				if(left->pred!=0)
					left->pred->succ=sPtr;
				else
					secondaryStructures=sPtr;
				delete left;
				}
			SecondaryStructure* right=sPtr->succ;
			if(right!=0&&right->structureType==newType)
				{
				/* Merge with structure to the right: */
				sPtr->residueEnd=right->residueEnd;
				sPtr->succ=right->succ;
				if(right->succ!=0)
					right->succ->pred=sPtr;
				delete right;
				}
			}
		
		/* Re-associate residues with secondary structures: */
		for(sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
			for(Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
				rPtr->secondaryStructure=sPtr;
		
		/* Re-associate dipoles with secondary structures: */
		std::vector<Dipole>::iterator nh=amides.begin();
		for(sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
			{
			sPtr->amidesBegin=nh;
			while(nh!=amides.end()&&nh->getMajorAtom()->residue->secondaryStructure==sPtr)
				++nh;
			sPtr->amidesEnd=nh;
			}
		std::vector<Dipole>::iterator co=carboxyls.begin();
		for(sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
			{
			sPtr->carboxylsBegin=co;
			while(co!=carboxyls.end()&&co->getMajorAtom()->residue->secondaryStructure==sPtr)
				++co;
			sPtr->carboxylsEnd=co;
			}
		}
	
	/* Update version number of secondary structure list: */
	++secondaryStructureVersion;
	}

Protein::StructureSelector Protein::pickStructure(int index)
	{
	SecondaryStructure* sPtr;
	int i;
	for(sPtr=secondaryStructures,i=0;sPtr!=0&&i!=index;sPtr=sPtr->succ,++i)
		;
	return StructureSelector(this,sPtr,i);
	}

Protein::StructureSelector Protein::pickStructure(Protein::SecondaryStructure::StructureType structureType,int index)
	{
	SecondaryStructure* sPtr=secondaryStructures;
	int i=0;
	while(true)
		{
		while(sPtr!=0&&sPtr->getStructureType()!=structureType)
			{
			sPtr=sPtr->succ;
			++i;
			}
		--index;
		if(sPtr==0||index<0)
			break;
		sPtr=sPtr->succ;
		++i;
		}
	return StructureSelector(this,sPtr,i);
	}

Protein::StructureSelector Protein::pickStructure(const Protein::Residue* residue)
	{
	int i=0;
	for(SecondaryStructure* sPtr=secondaryStructures;sPtr!=residue->getSecondaryStructure();++i,sPtr=sPtr->succ)
		;
	return StructureSelector(this,residue->getSecondaryStructure(),i);
	}

Protein::StructureSelector Protein::pickStructure(const Point& p)
	{
	SecondaryStructure* selected=0;
	int selectedIndex=-1;
	Scalar distance2Min=Math::Constants<Scalar>::max;
	
	/* Iterate through all atoms in the protein and check each against the query point: */
	int i=0;
	for(SecondaryStructure* sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ,++i)
		{
		for(Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
			{
			for(ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
				{
				/* Intersect the query ray with the covalent bond sphere around the atom: */
				Scalar distance2=Geometry::sqrDist(aPtr->getPosition(),p);
				if(distance2<=Math::sqr(aPtr->getCovalentRadius()*Scalar(1.5))&&distance2<distance2Min)
					{
					distance2Min=distance2;
					selected=sPtr;
					selectedIndex=i;
					}
				}
			}
		}
	
	return StructureSelector(this,selected,selectedIndex);
	}

Protein::StructureSelector Protein::pickStructure(const Ray& ray)
	{
	SecondaryStructure* selected=0;
	int selectedIndex=-1;
	Scalar lambdaMin=Math::Constants<Scalar>::max;
	
	/* Iterate through all atoms in the protein and intersect the query ray with each: */
	int i=0;
	for(SecondaryStructure* sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ,++i)
		{
		for(Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
			{
			for(ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
				{
				/* Intersect the query ray with the atom: */
				Scalar lambda=aPtr->intersectRay(ray);
				if(lambdaMin>lambda) // Intersection is the closest so far: */
					{
					lambdaMin=lambda;
					selected=sPtr;
					selectedIndex=i;
					}
				}
			}
		}
	
	return StructureSelector(this,selected,selectedIndex);
	}

Protein::BondRotator Protein::pickBond(const Point& p)
	{
	Scalar radius2=Math::sqr(bondRadius);
	
	/* Iterate through all backbone bonds and find the one closest to the given point: */
	std::vector<ChainAtom*>::iterator atom1,atom2;
	BondRotator::IntersectionType intersectionType=BondRotator::NONE;
	Scalar minDist2=Math::Constants<Scalar>::max;
	
	/***************************************************************************
	Due to the two-level hierarchy of storing atoms, iterating through all
	backbone bonds is slightly more complicated than one might initially expect.
	However, the fact that every third bond is a rigid peptide bond falls right
	out of the data structure, so that is a benefit at least.
	***************************************************************************/
	
	Residue* rPtr=residues;
	std::vector<ChainAtom*>::iterator atomIt=rPtr->backbone.begin();
	while(true)
		{
		/* Pick two neighbouring atoms: */
		bool peptideBond=false;
		std::vector<ChainAtom*>::iterator a1=atomIt;
		++atomIt;
		if(atomIt==rPtr->backbone.end()) // The second atom is in another residue
			{
			rPtr=rPtr->succ;
			if(rPtr==0) // Reached the end of the backbone
				break;
			atomIt=rPtr->backbone.begin();
			peptideBond=true; // This must therefore be a peptide bond
			}
		std::vector<ChainAtom*>::iterator a2=atomIt;
		if(!peptideBond) // Peptide bonds are rigid and cannot be picked!
			{
			/* Construct an upright cylinder connecting the two atoms, and test it against the point: */
			Point p1=Point((*a1)->getPosition());
			Point p2=Point((*a2)->getPosition());
			Vector d=p2-p1;
			Scalar dLen=Scalar(Geometry::mag(d));
			d/=dLen;
			Vector pp1=p-p1;
			Scalar beta=pp1*d;
			Scalar dist2=Geometry::sqr(pp1-beta*d);
			if(dist2<minDist2&&dist2<=radius2&&beta>=Scalar(0)&&beta<=dLen)
				{
				atom1=a1;
				atom2=a2;
				intersectionType=BondRotator::POINT_INSIDE;
				minDist2=dist2;
				}
			}
		}
	
	return BondRotator(this,atom1,atom2,intersectionType,p);
	}

Protein::BondRotator Protein::pickBond(const Ray& ray)
	{
	/* Iterate through all backbone bonds and find the closest intersection with the given ray: */
	std::vector<ChainAtom*>::iterator atom1,atom2;
	BondRotator::IntersectionType intersectionType=BondRotator::NONE;
	Scalar lambdaMin=Math::Constants<Scalar>::max;
	Point intersection;
	
	/***************************************************************************
	Due to the two-level hierarchy of storing atoms, iterating through all
	backbone bonds is slightly more complicated than one might initially expect.
	However, the fact that every third bond is a rigid peptide bond falls right
	out of the data structure, so that is a benefit at least.
	***************************************************************************/
	
	Residue* rPtr=residues;
	std::vector<ChainAtom*>::iterator atomIt=rPtr->backbone.begin();
	while(true)
		{
		/* Pick two neighbouring atoms: */
		bool peptideBond=false;
		std::vector<ChainAtom*>::iterator a1=atomIt;
		++atomIt;
		if(atomIt==rPtr->backbone.end()) // The second atom is in another residue
			{
			rPtr=rPtr->succ;
			if(rPtr==0) // Reached the end of the backbone
				break;
			atomIt=rPtr->backbone.begin();
			peptideBond=true; // This must therefore be a peptide bond
			}
		std::vector<ChainAtom*>::iterator a2=atomIt;
		if(!peptideBond) // Peptide bonds are rigid and cannot be picked!
			{
			/* Construct an upright cylinder connecting the two atoms, and intersect it with the ray: */
			Geometry::Cylinder<Scalar,3> bond((*a1)->getPosition(),(*a2)->getPosition(),bondRadius);
			Geometry::Cylinder<Scalar,3>::HitResult result=bond.intersectRay(ray);
			if(result.isValid()&&result.getParameter()<lambdaMin)
				{
				switch(result.getPart())
					{
					case Geometry::Cylinder<Scalar,3>::HitResult::MANTLE:
						intersectionType=BondRotator::SIDE;
						break;
					
					case Geometry::Cylinder<Scalar,3>::HitResult::BOTTOMCAP:
						intersectionType=BondRotator::BOTTOM;
						break;
					
					case Geometry::Cylinder<Scalar,3>::HitResult::TOPCAP:
						intersectionType=BondRotator::TOP;
						break;
					
					default:
						; // Just to make compiler happy
					}
				lambdaMin=result.getParameter();
				intersection=ray(lambdaMin);
				}
			}
		}
	
	return BondRotator(this,atom1,atom2,intersectionType,intersection);
	}

void Protein::getDihedralAngles(int startResidue,int numResidues,Scalar phis[],Scalar psis[]) const
	{
	/* Iterate to the first queried residue: */
	const Residue* rPtr=residues;
	for(int i=0;i<startResidue;++i)
		rPtr=rPtr->succ;
	
	/* Calculate all residues' phi and psi angles: */
	for(int i=0;i<numResidues;++i,rPtr=rPtr->succ)
		{
		DihedralAnglePair angles=rPtr->calcDihedralAngles();
		phis[i]=angles.phi;
		psis[i]=angles.psi;
		}
	}

void Protein::getDihedralAngles(int startResidue,int numResidues,Protein::DihedralAnglePair angles[]) const
	{
	/* Iterate to the first queried residue: */
	const Residue* rPtr=residues;
	for(int i=0;i<startResidue;++i)
		rPtr=rPtr->succ;
	
	/* Calculate all residues' phi and psi angles: */
	for(int i=0;i<numResidues;++i,rPtr=rPtr->succ)
		angles[i]=rPtr->calcDihedralAngles();
	}

void Protein::changeDihedralAngles(int startResidue,int numResidues,const Protein::DihedralAnglePair deltaAngles[],int direction,bool stopAtLastResidue)
	{
	typedef Geometry::AffineTransformation<Scalar,3> AtomTransform; // Type for transformations to apply along the chain
	
	/* Get pointers to the first and one-behind-last residues affected by the dihedral angle change: */
	Residue* r1=0;
	Residue* r2=0;
	
	/* Calculate cumulative transformation matrices along the backbone: */
	std::vector<AtomTransform> rotations;
	rotations.reserve(numResidues*2);
	if(direction>0)
		{
		/* Iterate to the first affected residue: */
		Residue* startResiduePtr=residues;
		for(int i=0;i<startResidue;++i)
			startResiduePtr=startResiduePtr->succ;
		
		/* Build rotation matrices: */
		AtomTransform currentRotation=AtomTransform::identity;
		Residue* rPtr=startResiduePtr;
		for(int i=0;i<numResidues;++i,rPtr=rPtr->succ)
			{
			if(rPtr->type==Residue::UNK) // Again, punt on unknown residue types
				{
				/* Create as many dummy rotations as there are backbone bonds inside this residue: */
				for(size_t j=1;j<rPtr->backbone.size();++j)
					rotations.push_back(currentRotation);
				}
			else
				{
				/* Get atom positions: */
				Point p[3];
				for(int j=0;j<3;++j)
					p[j]=rPtr->backbone[j]->getPosition();

				/* Concatenate phi rotation: */
				if(rPtr->type!=Residue::PRO) // Never rotate proline around its phi axis - it forms a rigid ring structure there
					currentRotation*=AtomTransform::rotateAround(p[0],AtomTransform::Rotation::rotateAxis(p[1]-p[0],deltaAngles[i].phi));
				rotations.push_back(currentRotation);

				/* Concatenate psi rotation: */
				currentRotation*=AtomTransform::rotateAround(p[1],AtomTransform::Rotation::rotateAxis(p[2]-p[1],deltaAngles[i].psi));
				rotations.push_back(currentRotation);
				}
			}
		
		/* Mark range of affected residues: */
		r1=startResiduePtr;
		if(stopAtLastResidue)
			r2=rPtr;
		else
			r2=0;
		
		/* Perform a breadth-first traversal on all atoms to the right of the start residue: */
		int lastRotation=rotations.size()-1;
		++rotationMarker;
		std::deque<MoleculeTraversal> traversalQueue;
		startResiduePtr->backbone[0]->rotationMarker=rotationMarker;
		startResiduePtr->backbone[1]->rotationMarker=rotationMarker;
		if(stopAtLastResidue&&rPtr!=0)
			rPtr->backbone[0]->rotationMarker=rotationMarker;
		traversalQueue.push_back(MoleculeTraversal(startResiduePtr->backbone[1],0));
		while(!traversalQueue.empty())
			{
			/* Pick the first element from the queue: */
			MoleculeTraversal mt=traversalQueue.front();
			traversalQueue.pop_front();

			/* Rotate it: */
			mt.atom->setPosition(rotations[mt.rotationIndex].transform(mt.atom->getPosition()));

			/* Put all its neighbours into the queue: */
			for(std::vector<Atom*>::iterator aIt=mt.atom->getBonds().begin();aIt!=mt.atom->getBonds().end();++aIt)
				{
				ChainAtom* aPtr=static_cast<ChainAtom*>(*aIt);
				if(aPtr->rotationMarker!=rotationMarker)
					{
					aPtr->rotationMarker=rotationMarker;
					if(mt.rotationIndex<lastRotation&&aPtr->backboneIndex>=1&&mt.atom->backboneIndex==aPtr->backboneIndex-1) // Is traversed bond a phi/psi backbone bond inside the residue sequence?
						traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex+1));
					else
						traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex));
					}
				}
			}
		}
	else
		{
		/* Iterate to the last affected residue: */
		Residue* finishResiduePtr=residues;
		for(int i=0;i<startResidue+numResidues-1;++i)
			finishResiduePtr=finishResiduePtr->succ;
		
		/* Build rotation matrices: */
		AtomTransform currentRotation=AtomTransform::identity;
		Residue* rPtr=finishResiduePtr;
		for(int i=numResidues-1;i>=0;--i,rPtr=rPtr->pred)
			{
			if(rPtr->type==Residue::UNK) // Again, punt on unknown residue types
				{
				/* Create as many dummy rotations as there are backbone bonds inside this residue: */
				for(size_t j=1;j<rPtr->backbone.size();++j)
					rotations.push_back(currentRotation);
				}
			else
				{
				/* Get atom positions: */
				Point p[3];
				for(int j=0;j<3;++j)
					p[j]=rPtr->backbone[j]->getPosition();
				
				/* Concatenate psi rotation: */
				currentRotation*=AtomTransform::rotateAround(p[2],AtomTransform::Rotation::rotateAxis(p[1]-p[2],deltaAngles[i].phi));
				rotations.push_back(currentRotation);
				
				/* Concatenate phi rotation: */
				if(rPtr->type!=Residue::PRO) // Never rotate proline around its phi axis - it forms a rigid ring structure there
					currentRotation*=AtomTransform::rotateAround(p[1],AtomTransform::Rotation::rotateAxis(p[0]-p[1],deltaAngles[i].psi));
				rotations.push_back(currentRotation);
				}
			}
		
		/* Mark range of affected residues: */
		r2=finishResiduePtr->succ;
		if(stopAtLastResidue&&rPtr!=0)
			r1=rPtr->succ;
		else
			r1=residues;
		
		/* Perform a breadth-first traversal on all atoms to the left of the finish residue: */
		int lastRotation=rotations.size()-1;
		++rotationMarker;
		std::deque<MoleculeTraversal> traversalQueue;
		finishResiduePtr->backbone[finishResiduePtr->backbone.size()-1]->rotationMarker=rotationMarker;
		finishResiduePtr->backbone[finishResiduePtr->backbone.size()-2]->rotationMarker=rotationMarker;
		if(stopAtLastResidue&&rPtr!=0)
			rPtr->backbone[rPtr->backbone.size()-1]->rotationMarker=rotationMarker;
		traversalQueue.push_back(MoleculeTraversal(finishResiduePtr->backbone[finishResiduePtr->backbone.size()-2],0));
		while(!traversalQueue.empty())
			{
			/* Pick the first element from the queue: */
			MoleculeTraversal mt=traversalQueue.front();
			traversalQueue.pop_front();

			/* Rotate it: */
			mt.atom->setPosition(rotations[mt.rotationIndex].transform(mt.atom->getPosition()));

			/* Put all its neighbours into the queue: */
			for(std::vector<Atom*>::iterator aIt=mt.atom->getBonds().begin();aIt!=mt.atom->getBonds().end();++aIt)
				{
				ChainAtom* aPtr=static_cast<ChainAtom*>(*aIt);
				if(aPtr->rotationMarker!=rotationMarker)
					{
					aPtr->rotationMarker=rotationMarker;
					if(mt.rotationIndex<lastRotation&&mt.atom->backboneIndex>=1&&aPtr->backboneIndex==mt.atom->backboneIndex-1) // Is traversed bond a phi/psi backbone bond inside the residue sequence?
						traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex+1));
					else
						traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex));
					}
				}
			}
		}
	
	/* Increment version numbers of all affected secondary structures: */
	SecondaryStructure* currentSecondaryStructure=0;
	for(Residue* rPtr=r1;rPtr!=r2;rPtr=rPtr->succ)
		if(rPtr->secondaryStructure!=currentSecondaryStructure)
			{
			currentSecondaryStructure=rPtr->secondaryStructure;
			++currentSecondaryStructure->version;
			}
	}

void Protein::setDihedralAngles(int startResidue,int numResidues,const Protein::DihedralAnglePair angles[],int direction)
	{
	/* Calculate the current dihedral angles: */
	DihedralAnglePair* deltaAngles=new DihedralAnglePair[numResidues];
	getDihedralAngles(startResidue,numResidues,deltaAngles);
	
	/* Calculate angle deltas: */
	for(int i=0;i<numResidues;++i)
		deltaAngles[i].delta(angles[i]);
	
	/* Update the dihedral angles: */
	changeDihedralAngles(startResidue,numResidues,deltaAngles,direction);
	
	/* Clean up: */
	delete[] deltaAngles;
	}

}
