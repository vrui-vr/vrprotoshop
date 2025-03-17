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

#ifndef PROTEININTERACTOR_INCLUDED
#define PROTEININTERACTOR_INCLUDED

#include <list>
#include <utility>
#include <Threads/Mutex.h>
#include <Threads/Cond.h>

#include "DragBox.h"
#include "Protein.h"
#include "ProteinRenderer.h"
#include "UndoBuffer.h"
#include "IK.h"
#include "EnergyAPI.h"

/* Forward declarations: */
namespace Cluster {
class MulticastPipe;
}
namespace Collab {
namespace Plugins {
class ProtoShopClient;
}
}

/* Namespace shortcuts: */
using Collab::Plugins::ProtoShopClient;

class ProteinInteractor
	{
	/* Embedded classes: */
	public:
	typedef IK::Transformation Transformation; // Type for transformations
	typedef std::list<MD::Protein::StructureSelector> StructureList; // Type for list of secondary structures
	typedef std::pair<int,MD::Protein::BackboneIterator::BondType> JointMapItem;
	
	enum UpdateDirection // Enumerated type for update directions
		{
		LEFT,RIGHT
		};
	
	private:
	class CoilRegionAngles // Class to hold dihedral angles for a coil region
		{
		/* Elements: */
		public:
		MD::Protein::StructureSelector structure; // Selector for the coil region
		MD::Protein::DihedralAnglePair* angles; // Pointer to C-style array of dihedral angle pairs
		
		/* Constructors and destructors: */
		CoilRegionAngles(const MD::Protein::StructureSelector& sStructure) // Creates angle arrays for secondary structure
			:structure(sStructure),
			 angles(new MD::Protein::DihedralAnglePair[structure.getNumResidues()])
			{
			};
		~CoilRegionAngles(void)
			{
			delete[] angles;
			};
		
		/* Methods: */
		int getFirstResidueIndex(void) const
			{
			return structure.getFirstResidueIndex();
			};
		int getNumResidues(void) const
			{
			return structure.getNumResidues();
			};
		void getCurrentAngles(void) // Retrieves current dihedral angles from structure
			{
			structure.getDihedralAngles(angles);
			};
		void changeDihedralAngles(UpdateDirection updateDirection)
			{
			structure.changeDihedralAngles(angles,updateDirection==LEFT?-1:1);
			};
		};
	
	class LockRequest // Base class for lock or lock update requests
		{
		/* Elements: */
		protected:
		unsigned int lockRequestId; // ID assigned to the server lock request
		
		/* Constructors and destructors: */
		LockRequest(unsigned int sLockRequestId)
			:lockRequestId(sLockRequestId)
			{
			}
		public:
		virtual ~LockRequest(void);
		
		/* Methods: */
		unsigned int getRequestId(void) const // Returns the request ID
			{
			return lockRequestId;
			}
		virtual void apply(ProteinInteractor* interactor) =0; // Applies a granted lock request to the given protein interactor
		};
	
	class SelectStructureLockRequest:public LockRequest
		{
		/* Elements: */
		private:
		MD::Protein::StructureSelector newSelectedStructure; // Selector for the structure that is to be selected
		
		/* Constructors and destructors: */
		public:
		SelectStructureLockRequest(unsigned int sLockRequestId,const MD::Protein::StructureSelector& sNewSelectedStructure)
			:LockRequest(sLockRequestId),
			 newSelectedStructure(sNewSelectedStructure)
			{
			}
		
		/* Methods from class LockRequest: */
		virtual void apply(ProteinInteractor* interactor);
		};
	
	class ToggleUpdateDirectionLockRequest:public LockRequest
		{
		/* Elements: */
		private:
		UpdateDirection newUpdateDirection; // New update direction
		
		/* Constructors and destructors: */
		public:
		ToggleUpdateDirectionLockRequest(unsigned int sLockRequestId,UpdateDirection sNewUpdateDirection)
			:LockRequest(sLockRequestId),
			 newUpdateDirection(sNewUpdateDirection)
			{
			}
		
		/* Methods from class LockRequest: */
		virtual void apply(ProteinInteractor* interactor);
		};
	
	class ToggleCoilLockRequest:public LockRequest
		{
		/* Elements: */
		private:
		MD::Protein::StructureSelector coil; // Selector for the coil region to be toggled
		
		/* Constructors and destructors: */
		public:
		ToggleCoilLockRequest(unsigned int sLockRequestId,const MD::Protein::StructureSelector& sCoil)
			:LockRequest(sLockRequestId),
			 coil(sCoil)
			{
			}
		
		/* Methods from class LockRequest: */
		virtual void apply(ProteinInteractor* interactor);
		};
	
	/* Elements: */
	private:
	unsigned int chainIndex; // Index of the polypeptide chain on which this interactor works
	MD::Protein* protein; // Pointer to the protein on which this interactor works
	ProtoShopClient* client; // Pointer to a ProtoShop protocol client, or null
	MD::ProteinRenderer* renderer; // Renderer associated with the protein worked on
	EnergyCalculator* energyCalculator; // Pointer to energy calculator associated with protein
	UndoBuffer* undoBuffer; // Pointer to the global undo buffer
	Cluster::MulticastPipe* multicastPipe; // Communications pipe for distributed protein modeling applications
	
	UpdateDirection updateDirection; // Default update direction when selecting a new structure
	
	MD::Protein::StructureSelector selectedStructure; // The selected secondary structure
	MD::Protein::Residue* selectedResidue; // A selected residue inside the protein
	
	DragBox box; // 3D interaction widget for the selected secondary structure
	StructureList activeCoils; // List of active coil regions in the protein
	
	LockRequest* lockRequest; // The state of a currently outstanding lock request, or 0
	unsigned int lockId; // ID of the lock held for the selected secondary structure and active coil regions, or 0
	
	/* Inverse kinematics state: */
	int numJoints;
	IK* ik; // Pointer to an IK object calculating dihedral angle updates for the protein
	JointMapItem* jointMap; // Array mapping from joint indices to backbone dihedral angles
	int firstResidueIndex;
	int numResidues;
	MD::Protein::DihedralAnglePair* currentAngles;
	MD::Protein::DihedralAnglePair* deltaAngles;
	int ikUpdateDirection;
	int effectorJointIndex;
	Transformation initialTransformation;
	Transformation correctionTransformation;
	IK::Scalar effectorStepSize;
	IK::Scalar correctionStepSize;
	
	/* Private methods: */
	std::pair<int,int> calcResidueRange(const MD::Protein::StructureSelector& structure,UpdateDirection direction) const; // Returns the residue range required by selecting the given structure with the given update direction
	std::pair<int,int> calcResidueRange(const MD::Protein::StructureSelector& coil) const; // Returns the residue range required by toggling the activation state of the given coil region
	void activateCoil(const MD::Protein::StructureSelector& coil); // Changes rendering state of a coil region to show it's active
	void deactivateCoil(const MD::Protein::StructureSelector& coil); // Changes rendering state of a coil region to show it's inactive
	void resetCoilRegions(void); // Resets the active coil regions to the default
	
	/* Pseudo-private methods to implement interactor state changes: */
	public:
	void doSelectStructure(const MD::Protein::StructureSelector& newSelectedtructure);
	void doSetUpdateDirection(UpdateDirection newUpdateDirection);
	void doToggleCoil(const MD::Protein::StructureSelector& coil);
	
	/* Constructors and destructors: */
	ProteinInteractor(unsigned int sChainIndex,MD::Protein* sProtein,MD::ProteinRenderer* sRenderer,EnergyCalculator* sEnergyCalculator,UndoBuffer* sUndoBuffer,Cluster::MulticastPipe* sMulticastPipe =0);
	~ProteinInteractor(void);
	
	/* Methods: */
	void collaborate(ProtoShopClient* newClient); // Enables collaborative protein editing using the given ProtoShop protocol client
	void lockReply(int firstResidueIndex,int numResidues,unsigned int newLockId); // Notifies the interactor that the ProtoShop protocol server has replied to a lock request
	void resetDragBox(void); // Resets the dragbox to the currently selected structure
	void toggleUpdateDirection(void); // Changes the default update direction
	int getUpdateDirection(void) const // Returns the current update direction for set/changeDihedralAngles
		{
		if(updateDirection==LEFT)
			return -1;
		else
			return 1;
		};
	void unselect(void); // Resets the currently selected secondary structure, all active coil regions, and the selected residue
	void selectStructure(const MD::Protein::StructureSelector& newSelectedStructure); // Selects a new secondary structure
	MD::Protein::StructureSelector getStructure(void) const // Returns the currently selected structure
		{
		return selectedStructure;
		};
	bool isValid(void) const // Returns true if a secondary structure is selected
		{
		return selectedStructure.isValid();
		};
	bool isAlphaHelix(void) const // Returns true if the currently selected structure is an alpha helix
		{
		return selectedStructure.isValid()&&selectedStructure.getStructureType()==MD::Protein::SecondaryStructure::ALPHA_HELIX;
		};
	bool isBetaStrand(void) const // Returns true if the currently selected structure is a beta strand
		{
		return selectedStructure.isValid()&&selectedStructure.getStructureType()==MD::Protein::SecondaryStructure::BETA_STRAND;
		};
	bool isCoil(void) const // Returns true if the currently selected structure is a coil region
		{
		return selectedStructure.isValid()&&selectedStructure.getStructureType()==MD::Protein::SecondaryStructure::COIL;
		};
	void setSelectedResidue(MD::Protein::Residue* newSelectedResidue); // Sets the selected residue
	MD::Protein::Residue* getSelectedResidue(void) const // Returns the selected residue
		{
		return selectedResidue;
		};
	DragBox& getBox(void) // Returns the drag box
		{
		return box;
		};
	const MD::Point& getBoxRotateCenter(void) const // Returns center of rotation of drag box
		{
		return box.getRotateCenter();
		};
	void setDihedralAngles(const MD::Protein::DihedralAnglePair angles[]); // Sets dihedral angles inside the selected structure
	void changeDihedralAngles(const MD::Protein::DihedralAnglePair deltaAngles[]); // Changes dihedral angles inside the selected structure
	void toggleCoil(const MD::Protein::StructureSelector& coil); // Toggles the "active" state of a coil region
	const std::list<MD::Protein::StructureSelector>& getLeftCoils(void) const // Returns left active coil regions
		{
		return activeCoils;
		};
	const std::list<MD::Protein::StructureSelector>& getRightCoils(void) const // Returns right active coil regions
		{
		return activeCoils;
		};
	bool startInteraction(void); // Prepares for performing IK on the active coil regions
	bool pickDragBox(const DragBox::Transformation& initialDragTransformation);
	bool pickDragBox(const DragBox::Point& rayStart,const DragBox::Vector& rayDirection,double* maxLambda =0);
	bool pickDragBox(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip);
	void dragBox(const DragBox::Point& rayStart,const DragBox::Vector& rayDirection)
		{
		box.drag(rayStart,rayDirection);
		};
	void dragBox(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip)
		{
		box.drag(modelView,projection,mouseClip);
		};
	void dragBox(const DragBox::Transformation& dragTransformation)
		{
		box.drag(dragTransformation);
		};
	const Transformation& getDragTransformation(void) const
		{
		return box.getDragTransformation();
		};
	void drag(const Transformation& goalTransformation); // Performs a step of IK on the active coil regions
	void applyChanges(void); // Applies IK changes to protein and sends state update packet to shared protein server on the given pipe if !=0
	void releaseDragBox(void)
		{
		box.release();
		};
	void finishInteraction(void); // Finishes performing IK
	void glRenderAction(GLContextData& contextData) const; // Renders the interaction object
	void glTransformProtein(GLContextData& contextData) const; // Loads a transformation matrix for visually offsetting the protein
	bool readCoilRegionAngles(const char* angleFilename); // Reads dihedral angles from a file
	bool writeCoilRegionAngles(const char* angleFilename) const; // Writes dihedral angles of active coil regions to a file
	};

#endif
