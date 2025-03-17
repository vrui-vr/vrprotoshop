/***********************************************************************
UndoBuffer - Class to encapsulate undo and redo functionality for
dihedral angle updates on proteins.
Copyright (c) 2002-2020 Oliver Kreylos

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

#include "UndoBuffer.h"

#include <Geometry/AffineTransformation.h>

/*******************************************
Methods of class UndoBuffer::UndoQueueEntry:
*******************************************/

void UndoBuffer::UndoQueueEntry::set(MD::Protein* sProtein,const Geometry::OrthonormalTransformation<double,3>& sProteinTransformation,const std::list<DihedralAngleSequence> previousSequences[2])
	{
	/* Copy parameters: */
	protein=sProtein;
	proteinTransformation=sProteinTransformation;
	sequences[0].clear();
	sequences[1].clear();
	
	/* Process both lists of sequences: */
	for(int i=0;i<2;++i)
		{
		/* Calculate and store dihedral angle changes for all left/right sequences: */
		for(std::list<DihedralAngleSequence>::const_iterator pIt=previousSequences[i].begin();pIt!=previousSequences[i].end();++pIt)
			{
			/* Add sequence to sequence list: */
			DihedralAngleSequence& deltas=*sequences[i].insert(sequences[i].end(),DihedralAngleSequence());

			/* Calculate (negative) angle changes: */
			deltas.resize(pIt->getFirstResidue(),pIt->getNumResidues());
			protein->getDihedralAngles(deltas.getFirstResidue(),deltas.getNumResidues(),deltas.angles);
			for(int j=0;j<deltas.getNumResidues();++j)
				deltas.angles[j].delta(pIt->angles[j]);
			}
		}
	}

void UndoBuffer::UndoQueueEntry::undo(void)
	{
	/* Apply inverse protein transformation: */
	protein->transform(Geometry::AffineTransformation<MD::Scalar,3>(Geometry::invert(proteinTransformation)));
	
	/* Process both lists of sequences: */
	for(int i=0;i<2;++i)
		{
		/* Determine update direction: */
		int direction=i==0?-1:1;
		
		/* Process all left/right sequences: */
		for(std::list<DihedralAngleSequence>::const_iterator pIt=sequences[i].begin();pIt!=sequences[i].end();++pIt)
			{
			/* Apply negative angle changes: */
			protein->changeDihedralAngles(pIt->getFirstResidue(),pIt->getNumResidues(),pIt->angles,direction);
			}
		}
	}

void UndoBuffer::UndoQueueEntry::redo(void)
	{
	/* Process both lists of sequences: */
	for(int i=0;i<2;++i)
		{
		/* Determine update direction: */
		int direction=i==0?-1:1;
		
		/* Process all left/right sequences: */
		for(std::list<DihedralAngleSequence>::const_iterator pIt=sequences[i].begin();pIt!=sequences[i].end();++pIt)
			{
			/* Negate angle deltas: */
			for(int j=0;j<pIt->getNumResidues();++j)
				pIt->angles[j].negate();
			
			/* Apply positive angle changes: */
			protein->changeDihedralAngles(pIt->getFirstResidue(),pIt->getNumResidues(),pIt->angles,direction);
			
			/* Negate angle deltas again: */
			for(int j=0;j<pIt->getNumResidues();++j)
				pIt->angles[j].negate();
			}
		}
	
	/* Apply protein transformation: */
	protein->transform(Geometry::AffineTransformation<MD::Scalar,3>(proteinTransformation));
	}

/***************************
Methods of class UndoBuffer:
***************************/

UndoBuffer::UndoBuffer(void)
	:queueTail(queue.end()),
	 currentProtein(0),lastSavedDepth(0)
	{
	}

UndoBuffer::~UndoBuffer(void)
	{
	}

void UndoBuffer::startInteraction(MD::Protein* protein)
	{
	/* Retrieve parameters from given structure: */
	currentProtein=protein;
	}

void UndoBuffer::addResidueSequence(int startResidue,int numResidues,int direction)
	{
	/* Create a single residue sequence: */
	int i=direction==-1?0:1;
	DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
	prev.resize(startResidue,numResidues);
	
	/* Retrieve current dihedral angles: */
	currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.angles);
	}

void UndoBuffer::startInteraction(MD::Protein* protein,MD::Protein::Residue* selectedResidue,int direction)
	{
	/* Retrieve parameters from given structure: */
	currentProtein=protein;
	
	/* Create a single residue sequence: */
	int i=direction==-1?0:1;
	DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
	prev.resize(protein->getResidueIndex(selectedResidue),1);
	
	/* Retrieve current dihedral angles: */
	currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.angles);
	}

void UndoBuffer::startInteraction(MD::Protein* protein,int startResidue,int numResidues,int direction)
	{
	/* Retrieve parameters from given structure: */
	currentProtein=protein;
	
	/* Create a single residue sequence: */
	int i=direction==-1?0:1;
	DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
	prev.resize(startResidue,numResidues);
	
	/* Retrieve current dihedral angles: */
	currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.angles);
	}

void UndoBuffer::startInteraction(const MD::Protein::StructureSelector& selectedStructure,int direction)
	{
	/* Retrieve parameters from given structure: */
	currentProtein=selectedStructure.getProtein();
	
	/* Create a single residue sequence: */
	int i=direction==-1?0:1;
	DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
	prev.resize(selectedStructure.getFirstResidueIndex(),selectedStructure.getNumResidues());
	
	/* Retrieve current dihedral angles: */
	currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.angles);
	}

void UndoBuffer::startInteraction(const std::list<MD::Protein::StructureSelector>& leftStructures,const std::list<MD::Protein::StructureSelector>& rightStructures)
	{
	/* Retrieve parameters from given structure: */
	currentProtein=0;
	
	/* Create left/right residue sequences: */
	for(int i=0;i<2;++i)
		{
		const std::list<MD::Protein::StructureSelector>& str=i==0?leftStructures:rightStructures;
		for(std::list<MD::Protein::StructureSelector>::const_iterator strIt=str.begin();strIt!=str.end();++strIt)
			{
			/* Get protein from structure (silently assume all structures belong to the same protein): */
			currentProtein=strIt->getProtein();
			
			/* Create a single residue sequence: */
			DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
			prev.resize(strIt->getFirstResidueIndex(),strIt->getNumResidues());

			/* Retrieve current dihedral angles: */
			currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.angles);
			}
		}
	}

void UndoBuffer::finishInteraction(void)
	{
	if(currentProtein!=0)
		{
		/* Insert a new undo queue entry before the current tail: */
		std::list<UndoQueueEntry>::iterator newIt=queue.insert(queueTail,UndoQueueEntry());

		/* Remove any queue entries behind the current tail (therefore flushing the redo buffer): */
		queue.erase(queueTail,queue.end());
		queueTail=queue.end();

		if(lastSavedDepth<0) // The last saved state has been flushed from the undo buffer
			{
			/* Invalidate the marker: */
			lastSavedDepth=queue.size(); // One-behind-last will never reach zero
			}
		else
			++lastSavedDepth;

		/* Set the undo queue entry's data: */
		newIt->set(currentProtein,Geometry::OrthonormalTransformation<double,3>::identity,previousSequences);
		}
	
	/* Clean up: */
	currentProtein=0;
	previousSequences[0].clear();
	previousSequences[1].clear();
	}

void UndoBuffer::finishInteraction(const Geometry::OrthonormalTransformation<double,3>& proteinTransformation)
	{
	if(currentProtein!=0)
		{
		/* Insert a new undo queue entry before the current tail: */
		std::list<UndoQueueEntry>::iterator newIt=queue.insert(queueTail,UndoQueueEntry());

		/* Remove any queue entries behind the current tail (therefore flushing the redo buffer): */
		queue.erase(queueTail,queue.end());
		queueTail=queue.end();

		if(lastSavedDepth<0) // The last saved state has been flushed from the undo buffer
			{
			/* Invalidate the marker: */
			lastSavedDepth=queue.size(); // One-behind-last will never reach zero
			}
		else
			++lastSavedDepth;

		/* Set the undo queue entry's data: */
		newIt->set(currentProtein,proteinTransformation,previousSequences);
		}
	
	/* Clean up: */
	currentProtein=0;
	previousSequences[0].clear();
	previousSequences[1].clear();
	}

MD::Protein* UndoBuffer::undo(void)
	{
	/* Rewind the queue tail by one and un-apply the skipped entry's dihedral angle change: */
	MD::Protein* result=0;
	if(queueTail!=queue.begin())
		{
		--queueTail;
		result=queueTail->getProtein();
		queueTail->undo();
		--lastSavedDepth;
		}
	
	return result;
	}

MD::Protein* UndoBuffer::redo(void)
	{
	/* Advance the queue tail by one and re-apply the skipped entry's dihedral angle change: */
	MD::Protein* result=0;
	if(queueTail!=queue.end())
		{
		result=queueTail->getProtein();
		queueTail->redo();
		++queueTail;
		++lastSavedDepth;
		}
	
	return result;
	}

void UndoBuffer::deleteProtein(const MD::Protein* protein)
	{
	std::list<UndoQueueEntry>::iterator qIt=queue.begin();
	while(qIt!=queue.end())
		{
		if(qIt->getProtein()==protein)
			{
			/* Delete this entry: */
			std::list<UndoQueueEntry>::iterator del=qIt;
			++qIt;
			queue.erase(del);
			}
		else
			++qIt;
		}
	}
