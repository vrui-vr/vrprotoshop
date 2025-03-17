/***********************************************************************
ProtoShopProtocol - Classes defining ProtoShop's client-server protocol.
Copyright (c) 2003-2021 Oliver Kreylos

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

#ifndef PROTOSHOPPROTOCOL_INCLUDED
#define PROTOSHOPPROTOCOL_INCLUDED

#include <string>
#include <Misc/SizedTypes.h>
#include <Misc/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Collaboration2/DataType.h>
#include <Collaboration2/Protocol.h>

#include "MDGeometry.h"
#include "Protein.h"

namespace Collab {

namespace Plugins {

class ProtoShopProtocol
	{
	/* Embedded classes: */
	protected:
	
	/* Protocol message IDs: */
	enum ClientMessages // Enumerated type for protocol message IDs sent by clients
		{
		LockRequest=0,
		LockUpdateRequest,
		IKUpdateRequest,
		TransformRequest,
		UnlockRequest,
		ResetTransformRequest,
		RenderFlagsUpdateRequest,
		
		NumClientMessages
		};
	
	enum ServerMessages // Enumerated type for protocol message IDs sent by servers
		{
		ConnectReply=0,
		LockReply,
		LockNotification,
		LockUpdateNotification,
		IKUpdateNotification,
		TransformNotification,
		UnlockNotification,
		ResetTransformNotification,
		RenderFlagsUpdateNotification,
		
		NumServerMessages
		};
	
	/*********************************************************************
	Data types used in protocol messages:
	*********************************************************************/
	
	struct Atom // Structure defining an atom in a polypeptide chain
		{
		/* Elements: */
		public:
		Misc::UInt8 type; // The atom's element type
		Misc::UInt16 index; // The atom's index in the polypeptide chain
		char placement[4]; // The atom's placement in its amino acid residue as a NUL-terminated string
		MD::Position position; // The atom's 3D position
		};
	
	typedef Misc::UInt16 ResidueIndex; // Type for amino acid residue indices
	
	struct Residue // Structure defining an amino acid residue in a polypeptide chain
		{
		/* Elements: */
		public:
		ResidueIndex index; // Residue's index from the input file
		char name[4]; // Residue's three-letter name from the input file as a NUL-terminated string
		Misc::Vector<Atom> atoms; // List of the residue's atoms
		};
	
	struct SecondaryStructure // Structure defining a secondary structure in a polypeptide chain
		{
		/* Elements: */
		public:
		Misc::UInt8 type; // Secondary structure type (0: None/unknown, 1: coil region, 2: alpha helix, 3: beta strand)
		Misc::Vector<Residue> residues; // List of the secondary structure's amino acid residues
		};
	
	public:
	typedef Geometry::OrthonormalTransformation<MD::Scalar,3> Transformation; // Type for 3D rigid body transformations
	
	struct RenderFlags // Structure holding per-chain rendering flags
		{
		/* Elements: */
		public:
		bool visible; // Main visibility flag
		bool drawAtoms;
		bool drawBonds;
		bool drawBackbone;
		bool drawCartoon;
		bool drawHydrogenBonds;
		bool drawHydrogenBondSites;
		bool drawHydrogenCages;
		bool drawAtomCollisions;
		};
	
	protected:
	struct Chain // Structure defining a polypeptide chain
		{
		/* Elements: */
		public:
		std::string name; // Chain's name
		Transformation transform; // The chain's transformation
		RenderFlags renderFlags; // The chain's rendering flags
		Misc::Vector<SecondaryStructure> secondaryStructures; // List of the polypeptide chain's secondary structures
		};
	
	typedef Misc::UInt8 ChainIndex; // Type for polypeptide chain indices
	
	struct Lock // Structure defining an update lock on a sequence of a amino acid residues in a polypeptide chain
		{
		/* Elements: */
		public:
		ChainIndex chain; // Index of the polypeptide chain to which the lock applies
		ResidueIndex begin; // Index of first locked residue
		ResidueIndex end; // Index one past last locked residue
		};
	
	typedef Misc::UInt8 LockID; // Type for lock IDs
	
	/*********************************************************************
	Structure declarations for shared protocol messages:
	*********************************************************************/
	
	struct IKUpdateMsg // Message structure shared by IKUpdateRequest and IKUpdateNotification
		{
		/* Elements: */
		public:
		LockID lockId; // The ID of the lock to which the update applies
		Misc::SInt8 updateDirection; // IK update direction (1: left-to-right / H-N to C=O terminus, -1: right-to-left / C=O to H-N terminus)
		bool limitUpdate; // Flag if only atoms within the affected residue sequence will be moved
		Misc::Vector<MD::Protein::DihedralAnglePair> dihedralAngles; // List of dihedral angles to apply to the locked residues
		};
	
	struct TransformMsg // Message structure shared by TransformRequest and TransformNotification
		{
		/* Elements: */
		public:
		LockID lockId; // The ID of the lock for the entire polypeptide chain
		Transformation transformation; // The chain's new transformation
		};
	
	struct UnlockMsg // Message structure shared by UnlockRequest and UnlockNotification
		{
		/* Elements: */
		public:
		LockID lockId; // The ID of the lock to be released
		};
	
	struct ResetTransformMsg // Message structure shared by ResetTransformRequest and ResetTransformNotification
		{
		/* Elements: */
		public:
		ChainIndex chain; // Index of the chain whose transformation to reset
		};
	
	struct RenderFlagsUpdateMsg // Message structure shared by RenderFlagsUpdateRequest and RenderFlagsUpdateNotification
		{
		/* Elements: */
		public:
		ChainIndex chain; // The index of the polypeptide chain whose flags are to be updated
		RenderFlags renderFlags; // The new rendering flags
		};
	
	/*********************************************************************
	Structure declarations for client protocol messages:
	*********************************************************************/
	
	struct LockRequestMsg
		{
		/* Elements: */
		public:
		LockID requestId; // Client-side ID to link lock requests and replies
		Lock lock; // The requested lock
		};
	
	struct LockUpdateRequestMsg
		{
		/* Elements: */
		public:
		LockID requestId; // Client-side ID to link lock requests and replies
		LockID lockId; // ID of the already-held lock
		Lock lock; // The requested lock
		};
	
	/*********************************************************************
	Structure declarations for server protocol messages:
	*********************************************************************/
	
	struct ConnectReplyMsg
		{
		/* Elements: */
		public:
		Misc::Vector<Chain> chains; // List of polypeptide chains in the workspace
		};
	
	struct LockReplyMsg
		{
		/* Elements: */
		public:
		LockID requestId; // Client-side ID to link lock requests and replies
		LockID lockId; // ID of the new lock if successful; 0 if unsuccessful
		};
	
	struct LockNotificationMsg
		{
		/* Elements: */
		public:
		LockID lockId; // ID of the new lock
		ClientID clientId; // ID of the client holding the new lock
		Lock lock; // The new lock
		};
	
	struct LockUpdateNotificationMsg
		{
		/* Elements: */
		public:
		LockID lockId; // ID of the updated lock
		Lock lock; // The new lock
		};
	
	/* Elements: */
	static const char* protocolName;
	static const unsigned int protocolVersion=1U<<16;
	
	/* Protocol data type declarations: */
	DataType protocolTypes; // Definitions of data types used by the NCK protocol
	DataType::TypeID scalarType; // Type for scalars
	DataType::TypeID positionType; // Type for 3D atom positions
	DataType::TypeID atomType; // Type for atom definitions
	DataType::TypeID residueIndexType; // Type for residue indices
	DataType::TypeID residueType; // Type for residues
	DataType::TypeID secondaryStructureType; // Type for secondary structures
	DataType::TypeID transformationType; // Type for 3D rigid body transformations
	DataType::TypeID renderFlagsType; // Type for rendering flags
	DataType::TypeID chainType; // Type for polypeptide chains
	DataType::TypeID chainIndexType; // Type for chain indices
	DataType::TypeID lockType; // Type for residue sequence locks
	DataType::TypeID dihedralAngleType; // Type for dihedral angles
	DataType::TypeID clientMessageTypes[NumClientMessages]; // Types for message structures sent by clients
	DataType::TypeID serverMessageTypes[NumServerMessages]; // Types for message structures sent by servers
	
	/* Constructors and destructors: */
	ProtoShopProtocol(void); // Creates protocol data types
	
	/* Methods: */
	size_t getClientMsgSize(unsigned int messageId) const // Returns the minimum size of a client protocol message
		{
		return protocolTypes.getMinSize(clientMessageTypes[messageId]);
		}
	size_t getServerMsgSize(unsigned int messageId) const // Returns the minimum size of a server protocol message
		{
		return protocolTypes.getMinSize(serverMessageTypes[messageId]);
		}
	};

}

}

#endif
