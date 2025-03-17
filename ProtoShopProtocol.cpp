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

#include "ProtoShopProtocol.h"

#include <Collaboration2/DataType.icpp>

namespace Collab {

namespace Plugins {

/******************************************
Static elements of class ProtoShopProtocol:
******************************************/

const char* ProtoShopProtocol::protocolName="ProtoShop";

/**********************************
Methods of class ProtoShopProtocol:
**********************************/

ProtoShopProtocol::ProtoShopProtocol(void)
	{
	/*********************************************************************
	Create basic protocol types:
	*********************************************************************/
	
	scalarType=DataType::getAtomicType<MD::Scalar>();
	DataType::TypeID scalar3Type=protocolTypes.createFixedArray(3,scalarType);
	positionType=scalar3Type;
	DataType::TypeID scalar4Type=protocolTypes.createFixedArray(4,scalarType);
	
	/* Struct Atom: */
	DataType::StructureElement atomElements[]=
		{
		{DataType::getAtomicType<Misc::UInt8>(),offsetof(Atom,type)},
		{DataType::getAtomicType<Misc::UInt16>(),offsetof(Atom,index)},
		{protocolTypes.createFixedArray(3,DataType::getAtomicType<char>()),offsetof(Atom,placement)},
		{positionType,offsetof(Atom,position)}
		};
	atomType=protocolTypes.createStructure(4,atomElements,sizeof(Atom));
	
	residueIndexType=DataType::getAtomicType<Misc::UInt16>();
	
	/* Struct residue: */
	DataType::StructureElement residueElements[]=
		{
		{residueIndexType,offsetof(Residue,index)},
		{protocolTypes.createFixedArray(3,DataType::getAtomicType<char>()),offsetof(Residue,name)},
		{protocolTypes.createVector(atomType),offsetof(Residue,atoms)}
		};
	residueType=protocolTypes.createStructure(3,residueElements,sizeof(Residue));
	
	/* Struct SecondaryStructure: */
	DataType::StructureElement secondaryStructureElements[]=
		{
		{DataType::getAtomicType<Misc::UInt8>(),offsetof(SecondaryStructure,type)},
		{protocolTypes.createVector(residueType),offsetof(SecondaryStructure,residues)}
		};
	secondaryStructureType=protocolTypes.createStructure(2,secondaryStructureElements,sizeof(SecondaryStructure));
	
	/* Class Transformation: */
	DataType::StructureElement transformationElements[]=
		{
		{scalar3Type,0},
		{scalar4Type,sizeof(Transformation::Vector)}
		};
	transformationType=protocolTypes.createStructure(2,transformationElements,sizeof(Transformation));
	
	/* Struct RenderFlags: */
	DataType::StructureElement renderFlagsElements[]=
		{
		{DataType::Bool,offsetof(RenderFlags,visible)},
		{DataType::Bool,offsetof(RenderFlags,drawAtoms)},
		{DataType::Bool,offsetof(RenderFlags,drawBonds)},
		{DataType::Bool,offsetof(RenderFlags,drawBackbone)},
		{DataType::Bool,offsetof(RenderFlags,drawCartoon)},
		{DataType::Bool,offsetof(RenderFlags,drawHydrogenBonds)},
		{DataType::Bool,offsetof(RenderFlags,drawHydrogenBondSites)},
		{DataType::Bool,offsetof(RenderFlags,drawHydrogenCages)},
		{DataType::Bool,offsetof(RenderFlags,drawAtomCollisions)}
		};
	renderFlagsType=protocolTypes.createStructure(9,renderFlagsElements,sizeof(RenderFlags));
	
	/* Struct Chain: */
	DataType::StructureElement chainElements[]=
		{
		{DataType::String,offsetof(Chain,name)},
		{transformationType,offsetof(Chain,transform)},
		{renderFlagsType,offsetof(Chain,renderFlags)},
		{protocolTypes.createVector(secondaryStructureType),offsetof(Chain,secondaryStructures)}
		};
	chainType=protocolTypes.createStructure(4,chainElements,sizeof(Chain));
	
	chainIndexType=DataType::getAtomicType<ChainIndex>();
	
	/* Struct Lock: */
	DataType::StructureElement lockElements[]=
		{
		{chainIndexType,offsetof(Lock,chain)},
		{residueIndexType,offsetof(Lock,begin)},
		{residueIndexType,offsetof(Lock,end)}
		};
	lockType=protocolTypes.createStructure(3,lockElements,sizeof(Lock));
	
	DataType::TypeID lockIdType=protocolTypes.getAtomicType<LockID>();
	
	/* Struct DihedralAnglePair: */
	DataType::StructureElement dihedralAnglePairElements[]=
		{
		{scalarType,offsetof(MD::Protein::DihedralAnglePair,phi)},
		{scalarType,offsetof(MD::Protein::DihedralAnglePair,psi)}
		};
	dihedralAngleType=protocolTypes.createStructure(2,dihedralAnglePairElements,sizeof(MD::Protein::DihedralAnglePair));
	
	/*********************************************************************
	Create types for shared protocol messages:
	*********************************************************************/
	
	/* IKUpdateRequest and IKUpdateNotification: */
	DataType::StructureElement ikUpdateMsgElements[]=
		{
		{lockIdType,offsetof(IKUpdateMsg,lockId)},
		{DataType::getAtomicType<Misc::SInt8>(),offsetof(IKUpdateMsg,updateDirection)},
		{DataType::Bool,offsetof(IKUpdateMsg,limitUpdate)},
		{protocolTypes.createVector(dihedralAngleType),offsetof(IKUpdateMsg,dihedralAngles)}
		};
	DataType::TypeID ikUpdateMsgType=protocolTypes.createStructure(4,ikUpdateMsgElements,sizeof(IKUpdateMsg));
	
	/* TransformRequest and TransformNotification: */
	DataType::StructureElement transformMsgElements[]=
		{
		{lockIdType,offsetof(TransformMsg,lockId)},
		{transformationType,offsetof(TransformMsg,transformation)}
		};
	DataType::TypeID transformMsgType=protocolTypes.createStructure(2,transformMsgElements,sizeof(TransformMsg));
	
	/* UnlockRequest and UnlockNotification: */
	DataType::StructureElement unlockMsgElements[]=
		{
		{lockIdType,offsetof(UnlockMsg,lockId)}
		};
	DataType::TypeID unlockMsgType=protocolTypes.createStructure(1,unlockMsgElements,sizeof(UnlockMsg));
	
	/* ResetTransformRequest and ResetTransformNotification: */
	DataType::StructureElement resetTransformMsgElements[]=
		{
		{chainIndexType,offsetof(ResetTransformMsg,chain)}
		};
	DataType::TypeID resetTransformMsgType=protocolTypes.createStructure(1,resetTransformMsgElements,sizeof(ResetTransformMsg));
	
	/* RenderFlagsUpdateRequest and RenderFlagsUpdateNotification: */
	DataType::StructureElement renderFlagsUpdateMsgElements[]=
		{
		{chainIndexType,offsetof(RenderFlagsUpdateMsg,chain)},
		{renderFlagsType,offsetof(RenderFlagsUpdateMsg,renderFlags)}
		};
	DataType::TypeID renderFlagsUpdateMsgType=protocolTypes.createStructure(2,renderFlagsUpdateMsgElements,sizeof(RenderFlagsUpdateMsg));
	
	/*********************************************************************
	Create types for client protocol messages:
	*********************************************************************/
	
	/* LockRequest: */
	DataType::StructureElement lockRequestMsgElements[]=
		{
		{lockIdType,offsetof(LockRequestMsg,requestId)},
		{lockType,offsetof(LockRequestMsg,lock)}
		};
	clientMessageTypes[LockRequest]=protocolTypes.createStructure(2,lockRequestMsgElements,sizeof(LockRequestMsg));
	
	/* LockUpdateRequest: */
	DataType::StructureElement lockUpdateRequestMsgElements[]=
		{
		{lockIdType,offsetof(LockUpdateRequestMsg,requestId)},
		{lockIdType,offsetof(LockUpdateRequestMsg,lockId)},
		{lockType,offsetof(LockUpdateRequestMsg,lock)}
		};
	clientMessageTypes[LockUpdateRequest]=protocolTypes.createStructure(3,lockUpdateRequestMsgElements,sizeof(LockUpdateRequestMsg));
	
	/* IKUpdateRequest: */
	clientMessageTypes[IKUpdateRequest]=ikUpdateMsgType;
	
	/* TransformRequest: */
	clientMessageTypes[TransformRequest]=transformMsgType;
	
	/* UnlockRequest: */
	clientMessageTypes[UnlockRequest]=unlockMsgType;
	
	/* ResetTransformRequest: */
	clientMessageTypes[ResetTransformRequest]=resetTransformMsgType;
	
	/* RenderFlagsUpdateRequest: */
	clientMessageTypes[RenderFlagsUpdateRequest]=renderFlagsUpdateMsgType;
	
	/*********************************************************************
	Create types for server protocol messages:
	*********************************************************************/
	
	/* ConnectReply: */
	DataType::StructureElement connectReplyMsgElements[]=
		{
		{protocolTypes.createVector(chainType),offsetof(ConnectReplyMsg,chains)}
		};
	serverMessageTypes[ConnectReply]=protocolTypes.createStructure(1,connectReplyMsgElements,sizeof(ConnectReplyMsg));
	
	/* LockReply: */
	DataType::StructureElement lockReplyMsgElements[]=
		{
		{lockIdType,offsetof(LockReplyMsg,requestId)},
		{lockIdType,offsetof(LockReplyMsg,lockId)}
		};
	serverMessageTypes[LockReply]=protocolTypes.createStructure(2,lockReplyMsgElements,sizeof(LockReplyMsg));
	
	/* LockNotification: */
	DataType::StructureElement lockNotificationMsgElements[]=
		{
		{lockIdType,offsetof(LockNotificationMsg,lockId)},
		{DataType::getAtomicType<ClientID>(),offsetof(LockNotificationMsg,clientId)},
		{lockType,offsetof(LockNotificationMsg,lock)}
		};
	serverMessageTypes[LockNotification]=protocolTypes.createStructure(3,lockNotificationMsgElements,sizeof(LockNotificationMsg));
	
	/* LockUpdateNotification: */
	DataType::StructureElement lockUpdateNotificationMsgElements[]=
		{
		{lockIdType,offsetof(LockUpdateNotificationMsg,lockId)},
		{lockType,offsetof(LockUpdateNotificationMsg,lock)}
		};
	serverMessageTypes[LockUpdateNotification]=protocolTypes.createStructure(2,lockUpdateNotificationMsgElements,sizeof(LockUpdateNotificationMsg));
	
	/* IKUpdateNotification: */
	serverMessageTypes[IKUpdateNotification]=ikUpdateMsgType;
	
	/* TransformNotification: */
	serverMessageTypes[TransformNotification]=transformMsgType;
	
	/* UnlockNotification: */
	serverMessageTypes[UnlockNotification]=unlockMsgType;
	
	/* ResetTransformNotification: */
	serverMessageTypes[ResetTransformNotification]=resetTransformMsgType;
	
	/* RenderFlagsUpdateNotification: */
	serverMessageTypes[RenderFlagsUpdateNotification]=renderFlagsUpdateMsgType;
	}

}

}
