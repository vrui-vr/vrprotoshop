/***********************************************************************
ProteinInteractor - Class encapsulating the manipulation state of a
protein.
Copyright (c) 2004-2021 Oliver Kreylos

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

#include "ProteinInteractor.h"

#include <stdio.h>
#include <Misc/Utility.h>
#include <Cluster/MulticastPipe.h>
#include <Math/Math.h>
#include <Geometry/PCACalculator.h>

#include "ProtoShopClient.h"

/* External global GUI functions: */
extern void updateProteinNow(void);

/***********************************************
Methods of class ProteinInteractor::LockRequest:
***********************************************/

ProteinInteractor::LockRequest::~LockRequest(void)
	{
	}

/**************************************************************
Methods of class ProteinInteractor::SelectStructureLockRequest:
**************************************************************/

void ProteinInteractor::SelectStructureLockRequest::apply(ProteinInteractor* interactor)
	{
	/* Select the new structure: */
	interactor->doSelectStructure(newSelectedStructure);
	}

/********************************************************************
Methods of class ProteinInteractor::ToggleUpdateDirectionLockRequest:
********************************************************************/

void ProteinInteractor::ToggleUpdateDirectionLockRequest::apply(ProteinInteractor* interactor)
	{
	/* Change the update direction: */
	interactor->doSetUpdateDirection(newUpdateDirection);
	}

/*********************************************************
Methods of class ProteinInteractor::ToggleCoilLockRequest:
*********************************************************/

void ProteinInteractor::ToggleCoilLockRequest::apply(ProteinInteractor* interactor)
	{
	/* Toggle the coil: */
	interactor->doToggleCoil(coil);
	}

/**********************************
Methods of class ProteinInteractor:
**********************************/

std::pair<int,int> ProteinInteractor::calcResidueRange(const MD::Protein::StructureSelector& structure,UpdateDirection direction) const
	{
	if(!structure.isValid())
		return std::pair<int,int>(-1,0);
	
	int resMin,resMax;
	resMin=structure.getFirstResidueIndex();
	resMax=resMin+structure.getNumResidues();
	
	/* Find the next coil region to the left: */
	MD::Protein::StructureSelector coil=structure;
	while(coil.isValid()&&coil.getStructureType()!=MD::Protein::SecondaryStructure::COIL)
		coil=direction==LEFT?coil.getPred():coil.getSucc();
	if(coil.isValid())
		resMin=direction==LEFT?coil.getFirstResidueIndex():coil.getFirstResidueIndex()+coil.getNumResidues();
	
	return std::pair<int,int>(resMin,resMax-resMin);
	}

std::pair<int,int> ProteinInteractor::calcResidueRange(const MD::Protein::StructureSelector& coil) const
	{
	int resMin,resMax;
	resMin=selectedStructure.getFirstResidueIndex();
	resMax=resMin+selectedStructure.getNumResidues();
	
	bool foundCoil=false;
	for(StructureList::const_iterator slIt=activeCoils.begin();slIt!=activeCoils.end();++slIt)
		{
		if((*slIt)!=coil)
			{
			if(resMin>slIt->getFirstResidueIndex())
				resMin=slIt->getFirstResidueIndex();
			else if(resMax<slIt->getFirstResidueIndex()+slIt->getNumResidues())
				resMax=slIt->getFirstResidueIndex()+slIt->getNumResidues();
			}
		else
			foundCoil=true;
		}
	if(!foundCoil)
		{
		if(resMin>coil.getFirstResidueIndex())
			resMin=coil.getFirstResidueIndex();
		else if(resMax<coil.getFirstResidueIndex()+coil.getNumResidues())
			resMax=coil.getFirstResidueIndex()+coil.getNumResidues();
		}
	
	return std::pair<int,int>(resMin,resMax-resMin);
	}

void ProteinInteractor::activateCoil(const MD::Protein::StructureSelector& coil)
	{
	#if 1
	renderer->selectStructure(coil);
	#else
	renderer->setDrawAtoms(coil,false);
	renderer->setDrawBonds(coil,true);
	renderer->setDrawBackbone(coil,false);
	renderer->setDrawBackboneRibbon(coil,false);
	renderer->setDrawCartoon(coil,false);
	#endif
	}

void ProteinInteractor::deactivateCoil(const MD::Protein::StructureSelector& coil)
	{
	#if 1
	renderer->deselectStructure(coil);
	#else
	renderer->setDrawAtoms(coil,renderer->getDrawAtoms());
	renderer->setDrawBonds(coil,false);
	renderer->setDrawBackbone(coil,renderer->getDrawBackbone());
	renderer->setDrawBackboneRibbon(coil,renderer->getDrawBackboneRibbon());
	renderer->setDrawCartoon(coil,renderer->getDrawCartoon());
	#endif
	}

void ProteinInteractor::resetCoilRegions(void)
	{
	/* Reset the active coil regions: */
	for(std::list<MD::Protein::StructureSelector>::iterator it=activeCoils.begin();it!=activeCoils.end();++it)
		deactivateCoil(*it);
	activeCoils.clear();
	
	/* Select the new active coil regions: */
	if(selectedStructure.isValid())
		{
		/* Set the default active coil region: */
		if(updateDirection==LEFT)
			{
			/* Find the next coil region to the left: */
			MD::Protein::StructureSelector coil=selectedStructure;
			while(coil.isValid()&&coil.getStructureType()!=MD::Protein::SecondaryStructure::COIL)
				coil=coil.getPred();
			if(coil.isValid())
				{
				activeCoils.push_back(coil);
				activateCoil(coil);
				}
			}
		else
			{
			/* Find the next coil region to the right: */
			MD::Protein::StructureSelector coil=selectedStructure;
			while(coil.isValid()&&coil.getStructureType()!=MD::Protein::SecondaryStructure::COIL)
				coil=coil.getSucc();
			if(coil.isValid())
				{
				activeCoils.push_back(coil);
				activateCoil(coil);
				}
			}
		}
	}

void ProteinInteractor::doSelectStructure(const MD::Protein::StructureSelector& newSelectedStructure)
	{
	/* Select the new structure: */
	selectedStructure=newSelectedStructure;
	
	/* Create a new drag box: */
	resetDragBox();
	
	/* Reset the active coil regions: */
	resetCoilRegions();
	}

void ProteinInteractor::doSetUpdateDirection(UpdateDirection newUpdateDirection)
	{
	/* Change the update direction: */
	updateDirection=newUpdateDirection;
	
	/* Reset the drag box: */
	resetDragBox();
	
	/* Reset the active coil regions: */
	resetCoilRegions();
	}

void ProteinInteractor::doToggleCoil(const MD::Protein::StructureSelector& coil)
	{
	/* Find coil region in active coil list: */
	std::list<MD::Protein::StructureSelector>::iterator it;
	for(it=activeCoils.begin();it!=activeCoils.end()&&it->getStructureIndex()<coil.getStructureIndex();++it)
		;
	if(it!=activeCoils.end()&&it->getStructureIndex()==coil.getStructureIndex())
		{
		/* Remove coil from list: */
		activeCoils.erase(it);
		deactivateCoil(coil);
		}
	else
		{
		/* Insert coil into list: */
		activeCoils.insert(it,coil);
		activateCoil(coil);
		}
	}

ProteinInteractor::ProteinInteractor(unsigned int sChainIndex,MD::Protein* sProtein,MD::ProteinRenderer* sRenderer,EnergyCalculator* sEnergyCalculator,UndoBuffer* sUndoBuffer,Cluster::MulticastPipe* sMulticastPipe)
	:chainIndex(sChainIndex),
	 protein(sProtein),
	 client(0),
	 renderer(sRenderer),
	 energyCalculator(sEnergyCalculator),
	 undoBuffer(sUndoBuffer),
	 multicastPipe(sMulticastPipe),
	 updateDirection(LEFT),
	 selectedResidue(0),
	 lockRequest(0),lockId(0),
	 ik(0),jointMap(0),
	 currentAngles(0),deltaAngles(0)
	{
	}

ProteinInteractor::~ProteinInteractor(void)
	{
	delete lockRequest;
	delete ik;
	delete[] jointMap;
	delete[] currentAngles;
	delete[] deltaAngles;
	}

void ProteinInteractor::collaborate(ProtoShopClient* newClient)
	{
	/* Set the ProtoShop client: */
	client=newClient;
	}

void ProteinInteractor::lockReply(int firstResidueIndex,int numResidues,unsigned int newLockId)
	{
	/* Check if the lock request was granted: */
	if(newLockId!=0)
		{
		/* Apply the lock request: */
		lockRequest->apply(this);
		
		/* Update the active lock ID: */
		lockId=newLockId;
		}
	
	/* Delete the lock request: */
	delete lockRequest;
	lockRequest=0;
	}

void ProteinInteractor::resetDragBox(void)
	{
	if(selectedStructure.isValid())
		{
		/* Create a drag box fitting the selected structure's backbone: */
		Geometry::PCACalculator<3> pca;
		
		/* Iterate through the selected structure's backbone bonds: */
		MD::Protein::BackboneIterator bIt=selectedStructure.beginBackbone();
		
		/* Use the second atom of the first backbone bond: */
		pca.accumulatePoint(bIt.getAtom2()->getPosition());
		
		/* Use the first atom of all subsequent backbone bonds: */
		for(++bIt;bIt!=selectedStructure.endBackbone();++bIt)
			pca.accumulatePoint(bIt.getAtom1()->getPosition());
		
		/* Calculate the backbone's centroid and eigenvectors: */
		MD::Point centroid=pca.calcCentroid();
		pca.calcCovariance();
		double eigenvalues[3];
		pca.calcEigenvalues(eigenvalues);
		MD::Vector axes[3];
		for(int i=0;i<3;++i)
			axes[i]=pca.calcEigenvector(eigenvalues[i]);
		
		/* Calculate the selected structure's extent along the three axes: */
		MD::Scalar sizes[3]={0,0,0};
		bIt=selectedStructure.beginBackbone();
		
		/* Use the second atom of the first backbone bond: */
		{
		MD::Vector d=bIt.getAtom2()->getPosition()-centroid;
		for(int i=0;i<3;++i)
			sizes[i]=Misc::max(sizes[i],Math::abs(axes[i]*d));
		}
		
		/* Use the first atom of all subsequent backbone bonds: */
		for(++bIt;bIt!=selectedStructure.endBackbone();++bIt)
			{
			MD::Vector d=bIt.getAtom1()->getPosition()-centroid;
			for(int i=0;i<3;++i)
				sizes[i]=Misc::max(sizes[i],Math::abs(axes[i]*d));
			}
		
		/* Update the drag box: */
		box.release();
		box.setCenter(centroid);
		for(int i=0;i<3;++i)
			{
			box.setAxis(i,axes[i]);
			box.setSize(i,sizes[i]*MD::Scalar(2)+MD::Scalar(2));
			}
		box.setEdgeRadius(MD::Protein::getBondRadius());
		
		/* Set the drag box's center of rotation: */
		if(selectedStructure.containsResidue(selectedResidue)&&selectedResidue->getType()!=MD::Protein::Residue::UNK)
			{
			if(selectedResidue->getType()==MD::Protein::Residue::PRO)
				{
				/* Set the drag box's center of rotation to the selected residue's hydrogen bond site: */
				MD::Protein::Dipole co=selectedResidue->getCarboxyl();
				box.setRotateCenter(co.getBondSite());
				}
			else
				{
				/* Set the drag box's center of rotation to the midpoint of the selected residue's hydrogen bond sites: */
				MD::Protein::Dipole nh=selectedResidue->getAmide();
				MD::Protein::Dipole co=selectedResidue->getCarboxyl();
				box.setRotateCenter(Geometry::mid(nh.getBondSite(),co.getBondSite()));
				}
			}
		else if(updateDirection==LEFT)
			{
			/* Set the drag box's center of rotation to the first structure atom: */
			MD::Protein::BackboneIterator bbIt=selectedStructure.beginBackbone();
			box.setRotateCenter(bbIt.getAtom1()->getPosition());
			}
		else
			{
			/* Set the drag box's center of rotation to the last structure atom: */
			MD::Protein::BackboneIterator bbIt=selectedStructure.endBackbone();
			box.setRotateCenter(bbIt.getAtom1()->getPosition());
			}
		
		/* Set drag box' picking mode: */
		box.setPickMode(DragBox::OPAQUE);
		}
	}

void ProteinInteractor::toggleUpdateDirection(void)
	{
	UpdateDirection newUpdateDirection=updateDirection==LEFT?RIGHT:LEFT;
	
	if(selectedStructure.isValid())
		{
		if(client!=0)
			{
			/* Bail out if there is an outstanding lock request: */
			if(lockRequest!=0)
				return;
			
			/* Calculate the residue range required to reverse the update direction: */
			std::pair<int,int> range=calcResidueRange(selectedStructure,newUpdateDirection);
			
			/* Ask the server to update the currently held lock to the new residue range: */
			unsigned int lockRequestId=client->updateLock(lockId,range.first,range.second);
			
			/* Create a lock request state: */
			lockRequest=new ToggleUpdateDirectionLockRequest(lockRequestId,newUpdateDirection);
			}
		else
			{
			/* Change the update direction: */
			doSetUpdateDirection(newUpdateDirection);
			}
		}
	else
		{
		/* Change the update direction: */
		updateDirection=newUpdateDirection;
		}
	}

void ProteinInteractor::unselect(void)
	{
	if(client!=0)
		{
		/* Cancel an outstanding lock request: */
		if(lockRequest!=0)
			client->cancelLockRequest(lockRequest->getRequestId());
		delete lockRequest;
		lockRequest=0;
		
		/* Release a currently-held lock: */
		client->releaseLock(lockId);
		lockId=0;
		}
	
	/* Reset the selected structure: */
	selectedStructure.invalidate();
	
	/* Reset all active coil regions: */
	resetCoilRegions();
	
	/* Reset the selected residue: */
	selectedResidue=0;
	}

void ProteinInteractor::selectStructure(const MD::Protein::StructureSelector& newSelectedStructure)
	{
	if(newSelectedStructure.getProtein()==protein)
		{
		if(client!=0)
			{
			/* Bail out if there is an outstanding lock request: */
			if(lockRequest!=0)
				return;
			
			/* Calculate the residue range required to select the new structure: */
			std::pair<int,int> range=calcResidueRange(newSelectedStructure,updateDirection);
			
			if(range.first>=0)
				{
				/* Check if there is an existing lock: */
				unsigned int lockRequestId;
				if(lockId!=0)
					{
					/* Ask the server to update the existing lock: */
					lockRequestId=client->updateLock(lockId,range.first,range.second);
					}
				else
					{
					/* Ask the server to lock the residue range: */
					lockRequestId=client->requestLock(chainIndex,range.first,range.second);
					}
				
				/* Create a lock request state: */
				lockRequest=new SelectStructureLockRequest(lockRequestId,newSelectedStructure);
				}
			else
				{
				if(lockId!=0)
					{
					/* Release the currently-held lock: */
					client->releaseLock(lockId);
					lockId=0;
					}
				
				/* Select the new (invalid) structure, which will de-select the currently selected structure: */
				doSelectStructure(newSelectedStructure);
				}
			}
		else
			{
			/* Select the new structure: */
			doSelectStructure(newSelectedStructure);
			}
		}
	}

void ProteinInteractor::setSelectedResidue(MD::Protein::Residue* newSelectedResidue)
	{
	/* Set the selected residue: */
	selectedResidue=newSelectedResidue;
	
	/* Reset the drag box to account for the new selected residue: */
	resetDragBox();
	}

void ProteinInteractor::setDihedralAngles(const MD::Protein::DihedralAnglePair angles[])
	{
	/* Bail out if the previous lock request has not been granted: */
	if(lockId==0)
		return;
	
	/* Get the affected residues' current dihedral angles: */
	int numResidues=selectedStructure.getNumResidues();
	MD::Protein::DihedralAnglePair* deltaAngles=new MD::Protein::DihedralAnglePair[numResidues];
	selectedStructure.getDihedralAngles(deltaAngles);
	
	/* Calculate and apply angle deltas to the selected structure: */
	for(int i=0;i<numResidues;++i)
		deltaAngles[i].delta(angles[i]);
	selectedStructure.changeDihedralAngles(deltaAngles,updateDirection==LEFT?1:-1);
	
	if(client!=0)
		{
		/* Send an angle update packet to the server: */
		// ...
		}
	
	/* Reset the drag box to account for the new angles: */
	resetDragBox();
	
	/* Clean up: */
	delete[] deltaAngles;
	}

void ProteinInteractor::changeDihedralAngles(const MD::Protein::DihedralAnglePair deltaAngles[])
	{
	/* Bail out if the previous lock request has not been granted: */
	if(lockId==0)
		return;
	
	/* Apply angle deltas to the selected structure: */
	selectedStructure.changeDihedralAngles(deltaAngles,updateDirection==LEFT?1:-1);
	
	if(client!=0)
		{
		/* Send an angle update packet to the server: */
		// ...
		}
	
	/* Reset the drag box to account for the new angles: */
	resetDragBox();
	}

void ProteinInteractor::toggleCoil(const MD::Protein::StructureSelector& coil)
	{
	if(selectedStructure.isValid()&&coil.isValid()&&coil.getProtein()==protein&&coil.getStructureType()==MD::Protein::SecondaryStructure::COIL)
		{
		if(client!=0)
			{
			/* Bail out if there is an outstanding lock request: */
			if(lockRequest!=0)
				return;
			
			/* Calculate the residue range required to toggle the coil region: */
			std::pair<int,int> range=calcResidueRange(coil);
			
			/* Ask the server to update the currently held lock the new residue range: */
			unsigned int lockRequestId=client->updateLock(lockId,range.first,range.second);
			
			/* Create a lock request state: */
			lockRequest=new ToggleCoilLockRequest(lockRequestId,coil);
			}
		else
			{
			/* Toggle the coil region: */
			doToggleCoil(coil);
			}
		}
	}

bool ProteinInteractor::startInteraction(void)
	{
	/* Bail out if the most recent lock request has not been granted: */
	if(client!=0&&lockId==0)
		return false;
	
	bool ok=true;
	
	if(multicastPipe==0||multicastPipe->isMaster())
		{
		/* Count number of rotatable backbone bonds in the update region: */
		numJoints=0;
		for(StructureList::iterator acIt=activeCoils.begin();acIt!=activeCoils.end();++acIt)
			for(MD::Protein::BackboneIterator bbIt=acIt->beginBackbone();bbIt!=acIt->endBackbone();++bbIt)
				if(bbIt.isRotatable())
					++numJoints;
		
		/* Create IK object: */
		ik=new IK(numJoints);
		
		/* Create translation structures from IK object to protein: */
		jointMap=new JointMapItem[numJoints];
		firstResidueIndex=activeCoils.front().getFirstResidueIndex();
		numResidues=activeCoils.back().getFirstResidueIndex()+activeCoils.back().getNumResidues()-firstResidueIndex;
		
		/* Determine the IK update direction: */
		if(activeCoils.back().getStructureIndex()<selectedStructure.getStructureIndex())
			{
			/* Unidirectional IK towards the C=O-terminus (to the right): */
			ikUpdateDirection=1;
			}
		else if(activeCoils.front().getStructureIndex()>selectedStructure.getStructureIndex())
			{
			/* Unidirectional IK towards the N-H-terminus (to the left): */
			ikUpdateDirection=-1;
			}
		else
			{
			/* Bidirectional IK: */
			ikUpdateDirection=updateDirection==LEFT?1:-1;
			}
		
		/* Initialize the joint chain: */
		if(ikUpdateDirection==1)
			{
			effectorJointIndex=numJoints;
			StructureList::iterator coilIt=activeCoils.begin();
			MD::Protein::BackboneIterator bbIt;
			int jointIndex=0;
			MD::Vector currentOffset=MD::Vector::zero;
			while(coilIt!=activeCoils.end())
				{
				if(effectorJointIndex==numJoints&&coilIt->getStructureIndex()>selectedStructure.getStructureIndex())
					effectorJointIndex=jointIndex;
				
				/* Add joints for all residues in this coil region: */
				bbIt=coilIt->beginBackbone();
				while(bbIt!=coilIt->endBackbone())
					{
					/* Add one rotation joint for each rotatable bond in the coil: */
					if(bbIt.isRotatable())
						{
						MD::Protein::BondRotator rotator(protein,bbIt);
						MD::Point rotPoint=rotator.getStart();
						MD::Point offset=rotPoint-currentOffset;
						currentOffset=rotPoint-MD::Point::origin;
						ik->setJoint(jointIndex,offset,rotator.getAxis());
						ik->setJointAngle(jointIndex,IK::Scalar(0));
						jointMap[jointIndex]=std::make_pair(bbIt.getResidueIndex()-firstResidueIndex,bbIt.getBondType());
						++jointIndex;
						}
					
					++bbIt;
					}
				
				++coilIt;
				}
			
			/* Set leaf offset: */
			MD::Point offset=bbIt.getAtom1()->getPosition()-currentOffset;
			ik->setLeafOffset(offset);
			}
		else
			{
			effectorJointIndex=numJoints;
			StructureList::iterator coilIt=activeCoils.end();
			MD::Protein::BackboneIterator bbIt;
			int jointIndex=0;
			MD::Vector currentOffset=MD::Vector::zero;
			while(coilIt!=activeCoils.begin())
				{
				--coilIt;
				
				if(effectorJointIndex==numJoints&&coilIt->getStructureIndex()<selectedStructure.getStructureIndex())
					effectorJointIndex=jointIndex;
				
				/* Add joints for all residues in this coil region: */
				bbIt=coilIt->endBackbone();
				while(bbIt!=coilIt->beginBackbone())
					{
					--bbIt;
					
					/* Add one rotation joint for each rotatable bond in the coil: */
					if(bbIt.isRotatable())
						{
						MD::Protein::BondRotator rotator(protein,bbIt);
						MD::Point rotPoint=rotator.getEnd();
						MD::Point offset=rotPoint-currentOffset;
						currentOffset=rotPoint-MD::Point::origin;
						ik->setJoint(jointIndex,offset,-rotator.getAxis());
						ik->setJointAngle(jointIndex,IK::Scalar(0));
						jointMap[jointIndex]=std::make_pair(bbIt.getResidueIndex()-firstResidueIndex,bbIt.getBondType());
						++jointIndex;
						}
					}
				}
			
			/* Set leaf offset: */
			MD::Point offset=bbIt.getAtom2()->getPosition()-currentOffset;
			ik->setLeafOffset(offset);
			}
		
		/* Initialize IK parameters: */
		IK::Scalar orientationWeights[3]={10,10,10};
		ik->setOrientationWeights(orientationWeights);
		ik->setMaxResidual(IK::Scalar(1.0e-12));
		ik->setMaxStepError(IK::Scalar(1.0e-3));
		effectorStepSize=IK::Scalar(1.0e-4);
		correctionStepSize=IK::Scalar(1.0e-4);

		/* Calculate the initial effector transformation: */
		ik->setEffectorJointIndex(effectorJointIndex);
		initialTransformation=ik->getEffectorTransformation();
		ik->setEffectorJointIndex(numJoints);
		correctionTransformation=ik->getEffectorTransformation();
		}
	
	if(multicastPipe!=0)
		{
		/* Share interaction state with the rendering cluster: */
		multicastPipe->broadcast<int>(numJoints);
		multicastPipe->broadcast<int>(firstResidueIndex);
		multicastPipe->broadcast<int>(numResidues);
		multicastPipe->broadcast<int>(ikUpdateDirection);
		multicastPipe->broadcast<int>(effectorJointIndex);
		multicastPipe->flush();
		}
	
	/* Initialize angle arrays: */
	currentAngles=new MD::Protein::DihedralAnglePair[numResidues];
	deltaAngles=new MD::Protein::DihedralAnglePair[numResidues];
	for(int i=0;i<numResidues;++i)
		{
		currentAngles[i]=MD::Protein::DihedralAnglePair(0,0);
		deltaAngles[i]=MD::Protein::DihedralAnglePair(0,0);
		}
	
	/* Open undo buffer: */
	undoBuffer->startInteraction(protein,firstResidueIndex,numResidues,ikUpdateDirection);
	
	return ok;
	}

bool ProteinInteractor::pickDragBox(const DragBox::Transformation& initialDragTransformation)
	{
	/* Try picking the interaction box: */
	if(!selectedStructure.isValid()||!box.pick(initialDragTransformation))
		return false;
	
	bool ok=startInteraction();
	if(!ok)
		box.release();
	return ok;
	}

bool ProteinInteractor::pickDragBox(const DragBox::Point& rayStart,const DragBox::Vector& rayDirection,double* maxLambda)
	{
	/* Try picking the interaction box: */
	if(!selectedStructure.isValid()||!box.pick(rayStart,rayDirection,maxLambda))
		return false;
	
	bool ok=startInteraction();
	if(!ok)
		box.release();
	return ok;
	}

bool ProteinInteractor::pickDragBox(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip)
	{
	/* Try picking the interaction box: */
	if(!selectedStructure.isValid()||!box.pick(modelView,projection,mouseClip))
		return false;
	
	bool ok=startInteraction();
	if(!ok)
		box.release();
	return ok;
	}

void ProteinInteractor::drag(const ProteinInteractor::Transformation& goalTransformation)
	{
	if(multicastPipe==0||multicastPipe->isMaster())
		{
		/* Calculate the difference transformation: */
		Transformation ikTransformation=goalTransformation*initialTransformation;
		
		/*********************************************************************
		Perform IK steps on the joint chain. First, move the dragged structure
		towards its goal position. In the case of bi-directional IK, this will
		incur (small) drift in the leaf joint in the joint chain. Therefore,
		only in the case of bi-directional IK, the leaf joint is moved back
		towards its original position in a second step.
		There are two step sizes involved here: The first is the step size for
		moving the dragged structure, the second is the step size for moving
		the leaf joint. Since these IK operations work on different scales, we
		maintain the two step sizes independently.
		*********************************************************************/
		
		/* Move the dragged structure towards its goal position: */
		ik->setEffectorJointIndex(effectorJointIndex); // Select the dragged structure for IK
		ik->setStepSize(effectorStepSize); // Set the last used IK step size
		ik->performSteps(ikTransformation,2000); // Perform the IK steps and save the residual
		effectorStepSize=ik->getStepSize(); // Save the IK step size
		
		/* Check for bi-directional IK: */
		if(effectorJointIndex<numJoints)
			{
			/* Move the chain's leaf joint back towards its original position: */
			ik->setEffectorJointIndex(numJoints); // Select the chain's leaf joint
			ik->setStepSize(correctionStepSize); // Set the last used IK step size
			ik->performSteps(correctionTransformation,2000); // Perform the IK steps
			correctionStepSize=ik->getStepSize(); // Save the IK step size
			}
		
		/* Update the coils' dihedral angles according to the joint tree: */
		for(int i=0;i<numJoints;++i)
			{
			IK::Scalar newAngle=ik->getJointAngle(i);
			int angleIndex=jointMap[i].first;
			if(jointMap[i].second==MD::Protein::BackboneIterator::PHI)
				{
				deltaAngles[angleIndex].phi+=newAngle-currentAngles[angleIndex].phi;
				currentAngles[angleIndex].phi=newAngle;
				}
			else
				{
				deltaAngles[angleIndex].psi+=newAngle-currentAngles[angleIndex].psi;
				currentAngles[angleIndex].psi=newAngle;
				}
			}
		}
	
	if(multicastPipe!=0)
		{
		/* Broadcast the dihedral angle changes: */
		multicastPipe->broadcast<MD::Protein::DihedralAnglePair>(deltaAngles,numResidues);
		multicastPipe->flush();
		}
	}

void ProteinInteractor::applyChanges(void)
	{
	if(client!=0)
		{
		/* Send a state update packet to the shared protein server: */
		client->ikUpdate(lockId,deltaAngles,ikUpdateDirection,effectorJointIndex<numJoints);
		}
	
	/* Apply dihedral angle changes to the protein: */
	protein->changeDihedralAngles(firstResidueIndex,numResidues,deltaAngles,ikUpdateDirection,effectorJointIndex<numJoints);
	
	/* Reset dihedral angle differences: */
	for(int i=0;i<numResidues;++i)
		deltaAngles[i]=MD::Protein::DihedralAnglePair(0,0);
	
	/* Update related data structures: */
	renderer->updateProtein();
	if(energyCalculator!=0)
		energyCalculator->updateProtein();
	}

void ProteinInteractor::finishInteraction(void)
	{
	/* Close undo buffer: */
	undoBuffer->finishInteraction();
	
	/* Delete IK data structures: */
	delete ik;
	ik=0;
	delete[] jointMap;
	jointMap=0;
	delete[] currentAngles;
	currentAngles=0;
	delete[] deltaAngles;
	deltaAngles=0;
	
	/* Snap drag box back to structure: */
	resetDragBox();
	}

void ProteinInteractor::glRenderAction(GLContextData& contextData) const
	{
	if(selectedStructure.isValid())
		{
		/* Render the drag box: */
		glColor4f(0.0f,0.5f,0.0f,0.33f);
		glLineWidth(3.0);
		glPointSize(5.0);
		glDisable(GL_LIGHTING);
		glDisable(GL_CULL_FACE);
		box.draw();
		glEnable(GL_CULL_FACE);
		glEnable(GL_LIGHTING);
		}
	
	if(selectedResidue!=0)
		{
		/* Highlight the selected residue: */
		renderer->highlightResidue(contextData,selectedResidue);
		}
	}

void ProteinInteractor::glTransformProtein(GLContextData& contextData) const
	{
	}

bool ProteinInteractor::readCoilRegionAngles(const char* angleFilename)
	{
	/* Open input file: */
	FILE* angleFile=fopen(angleFilename,"rt");
	if(angleFile==0)
		return false;
	
	/* Save previous protein state: */
	undoBuffer->startInteraction(protein);
	
	/* Read all coil regions from file: */
	int numAnimationSteps=100;
	std::list<CoilRegionAngles*> coilRegionAngleDeltas;
	while(true)
		{
		/* Read coil description from file: */
		int structureIndex,firstResidueIndex,numResidues;
		if(fscanf(angleFile,"%d %d %d",&structureIndex,&firstResidueIndex,&numResidues)!=3)
			break;
		
		/* Get appropriate structure selector and check consistency: */
		MD::Protein::StructureSelector structure=protein->pickStructure(structureIndex);
		if(structure.getFirstResidueIndex()==firstResidueIndex&&structure.getNumResidues()==numResidues)
			{
			/* Save current angles in undo buffer: */
			undoBuffer->addResidueSequence(firstResidueIndex,numResidues,updateDirection==LEFT?-1:1);
			
			/* Retrieve current angles from structure: */
			CoilRegionAngles currentAngles(structure);
			currentAngles.getCurrentAngles();
			
			/* Read target dihedral angles: */
			CoilRegionAngles targetAngles(structure);
			for(int i=0;i<numResidues;++i)
				{
				double phi,psi;
				if(fscanf(angleFile,"%lf %lf",&phi,&psi)==2)
					{
					targetAngles.angles[i].phi=Math::rad(MD::Scalar(phi));
					targetAngles.angles[i].psi=Math::rad(MD::Scalar(psi));
					}
				}
			
			/* Calculate and save angle deltas: */
			CoilRegionAngles* angleDeltas=new CoilRegionAngles(structure);
			for(int i=0;i<angleDeltas->getNumResidues();++i)
				{
				double dPhi=(targetAngles.angles[i].phi-currentAngles.angles[i].phi);
				if(dPhi<=-Math::Constants<double>::pi)
					dPhi+=2.0*Math::Constants<double>::pi;
				else if(dPhi>Math::Constants<double>::pi)
					dPhi-=2.0*Math::Constants<double>::pi;
				angleDeltas->angles[i].phi=MD::Scalar(dPhi/double(numAnimationSteps));
				double dPsi=(targetAngles.angles[i].psi-currentAngles.angles[i].psi);
				if(dPsi<=-Math::Constants<double>::pi)
					dPsi+=2.0*Math::Constants<double>::pi;
				else if(dPsi>Math::Constants<double>::pi)
					dPsi-=2.0*Math::Constants<double>::pi;
				angleDeltas->angles[i].psi=MD::Scalar(dPsi/double(numAnimationSteps));
				}
			coilRegionAngleDeltas.push_back(angleDeltas);
			}
		}
	
	fclose(angleFile);
	
	/* Morph protein from current coil region angles to target angles: */
	for(int i=0;i<numAnimationSteps;++i)
		{
		for(std::list<CoilRegionAngles*>::iterator craIt=coilRegionAngleDeltas.begin();craIt!=coilRegionAngleDeltas.end();++craIt)
			(*craIt)->changeDihedralAngles(updateDirection);
		updateProteinNow();
		}
	if(energyCalculator!=0)
		energyCalculator->updateProtein();
	
	/* Delete angle deltas: */
	for(std::list<CoilRegionAngles*>::iterator craIt=coilRegionAngleDeltas.begin();craIt!=coilRegionAngleDeltas.end();++craIt)
		delete *craIt;
	coilRegionAngleDeltas.clear();
	
	/* Finish interaction: */
	undoBuffer->finishInteraction();
	
	/* Snap drag box back to structure: */
	resetDragBox();
	
	return true;
	}

bool ProteinInteractor::writeCoilRegionAngles(const char* angleFilename) const
	{
	/* Open output file: */
	FILE* angleFile=fopen(angleFilename,"wt");
	if(angleFile==0)
		return false;
	
	/* Iterate through all active coils: */
	for(std::list<MD::Protein::StructureSelector>::const_iterator cIt=activeCoils.begin();cIt!=activeCoils.end();++cIt)
		{
		/* Write coil region header: */
		fprintf(angleFile,"%d %d %d\n",cIt->getStructureIndex(),cIt->getFirstResidueIndex(),cIt->getNumResidues());
		
		/* Retrieve dihedral angles: */
		MD::Protein::DihedralAnglePair* angles=new MD::Protein::DihedralAnglePair[cIt->getNumResidues()];
		cIt->getDihedralAngles(angles);
		
		/* Write dihedral angles: */
		for(int i=0;i<cIt->getNumResidues();++i)
			fprintf(angleFile,"%10.4f %10.4f\n",double(Math::deg(angles[i].phi)),double(Math::deg(angles[i].psi)));
		delete[] angles;
		}
	
	fclose(angleFile);
	return true;
	}
