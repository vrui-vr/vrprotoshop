/***********************************************************************
ProtoShopClient - Class representing a client for a ProtoShop workspace.
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

#ifndef PROTOSHOPCLIENT_INCLUDED
#define PROTOSHOPCLIENT_INCLUDED

#include <string>
#include <vector>
#include <Misc/StandardHashFunction.h>
#include <Misc/HashTable.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Collaboration2/DataType.icpp>
#include <Collaboration2/MessageBuffer.h>
#include <Collaboration2/MessageWriter.h>
#include <Collaboration2/Client.h>
#include <Collaboration2/PluginClient.h>

#include "ProtoShopProtocol.h"
#include "VRProtoShop.h"

/* Forward declarations: */
namespace Collab {
class MessageReader;
class MessageContinuation;
}
namespace MD {
class Protein;
}

namespace Collab {

namespace Plugins {

class ProtoShopClient:public PluginClient,public ProtoShopProtocol
	{
	/* Embedded classes: */
	private:
	struct ServerChain // Structure representing a polypeptide chain received from the server
		{
		/* Elements: */
		public:
		std::string name;
		Transformation transform;
		RenderFlags renderFlags;
		MD::Protein* protein;
		};
	
	struct LockRequestState // Structure defining an outstanding lock request
		{
		/* Elements: */
		public:
		unsigned int chainIndex; // Index of the polypeptide chain to which the lock applies
		unsigned int beginIndex; // Index of the first locked residue
		unsigned int endIndex; // Index one after the last locked residue
		VRProtoShop::ChainDraggingTool* chainDraggingTool; // Pointer to a chain dragging tool if this request is related to a polypeptide dragging operation
		
		/* Constructors and destructors: */
		LockRequestState(unsigned int sChainIndex,unsigned int sBeginIndex,unsigned int sEndIndex,VRProtoShop::ChainDraggingTool* sChainDraggingTool =0)
			:chainIndex(sChainIndex),beginIndex(sBeginIndex),endIndex(sEndIndex),
			 chainDraggingTool(sChainDraggingTool)
			{
			}
		};
	
	typedef Misc::HashTable<LockID,LockRequestState> LockRequestMap; // Type for hash tables mapping lock IDs to outstanding lock requests
	
	struct LockState // Structure defining the state of an active residue sequence lock
		{
		/* Elements: */
		public:
		ClientID clientId; // ID of the client holding the lock
		unsigned int chainIndex; // Index of the polypeptide chain to which the lock applies
		unsigned int beginIndex; // Index of the first locked residue
		unsigned int endIndex; // Index one after the last locked residue
		};
	
	typedef Misc::HashTable<LockID,LockState> LockMap; // Type for hash tables mapping lock IDs to lock states
	
	/* Elements: */
	VRProtoShop* application; // Pointer back to the VRProtoShop application object
	LockID lastLockRequestId; // ID of the last lock request
	LockRequestMap lockRequests; // Map of outstanding lock requests
	LockMap locks; // Map of active locks
	
	/* Private methods: */
	void queueServerMessage(unsigned int messageId,const void* messageStructure) // Writes a message structure for the given message ID into a buffer and queues to send it to the server
		{
		/* Create a message writer: */
		MessageWriter message(MessageBuffer::create(clientMessageBase+messageId,protocolTypes.calcSize(clientMessageTypes[messageId],messageStructure)));
		
		/* Write the message structure into the message: */
		protocolTypes.write(clientMessageTypes[messageId],messageStructure,message);
		
		/* Queue the message for delivery to the server: */
		client->queueServerMessage(message.getBuffer());
		}
	
	/* Message handling methods: */
	void frontendConnectReplyCallback(unsigned int messageId,MessageReader& message);
	MessageContinuation* connectReplyCallback(unsigned int messageId,MessageContinuation* continuation);
	void frontendLockReplyCallback(unsigned int messageId,MessageReader& message);
	void frontendLockNotificationCallback(unsigned int messageId,MessageReader& message);
	void frontendLockUpdateNotificationCallback(unsigned int messageId,MessageReader& message);
	void frontendIKUpdateNotificationCallback(unsigned int messageId,MessageReader& message);
	MessageContinuation* ikUpdateNotificationCallback(unsigned int messageId,MessageContinuation* continuation);
	void frontendTransformNotificationCallback(unsigned int messageId,MessageReader& message);
	void frontendUnlockNotificationCallback(unsigned int messageId,MessageReader& message);
	void frontendResetTransformNotificationCallback(unsigned int messageId,MessageReader& message);
	void frontendRenderFlagsUpdateNotificationCallback(unsigned int messageId,MessageReader& message);
	
	/* Constructors and destructors: */
	public:
	ProtoShopClient(Client* sClient,VRProtoShop* sApplication); // Creates a ProtoShop client
	virtual ~ProtoShopClient(void);
	
	/* Methods from class PluginClient: */
	virtual const char* getName(void) const;
	virtual unsigned int getVersion(void) const;
	virtual unsigned int getNumClientMessages(void) const;
	virtual unsigned int getNumServerMessages(void) const;
	virtual void setMessageBases(unsigned int newClientMessageBase,unsigned int newServerMessageBase);
	virtual void start(void);
	
	/* New methods: */
	LockID requestLock(unsigned int chainIndex,unsigned int firstResidueIndex,unsigned int numResidues); // Requests a lock for the given sequence of residues; returns the lock request ID
	LockID requestLock(unsigned int chainIndex,VRProtoShop::ChainDraggingTool* chainDraggingTool); // Requests a lock for an entire polypeptide chain about to be dragged; returns the lock request ID
	LockID updateLock(LockID lockId,unsigned int firstResidueIndex,unsigned int numResidues); // Updates an existing lock to the given new sequence of residues on the same chain; returns the lock request ID
	void cancelLockRequest(LockID lockRequestId); // Cancels a pending lock request
	void dragChain(LockID lockId,const Transformation& newTransform); // Updates the transformation of the polypeptide chain locked by the given lock ID
	void ikUpdate(LockID lockId,const MD::Protein::DihedralAnglePair deltaAngles[],int updateDirection,bool limitUpdate); // Applies an IK update to the polypeptide chain locked by the given lock ID
	void releaseLock(LockID lockId); // Releases a currently-held lock
	void resetChainTransform(unsigned int chainIndex); // Resets the transformation of the polypeptide chain with the given index
	void updateRenderFlags(void); // Called after the rendering flags for a polypeptide chain have been changed in the GUI
	};

}

}

#endif
