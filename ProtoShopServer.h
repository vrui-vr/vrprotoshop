/***********************************************************************
ProtoShopServer - Class representing a server for a ProtoShop workspace.
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

#ifndef PROTOSHOPSERVER_INCLUDED
#define PROTOSHOPSERVER_INCLUDED

#include <string>
#include <vector>
#include <Misc/SimpleSet.h>
#include <Misc/StandardHashFunction.h>
#include <Misc/HashTable.h>
#include <Collaboration2/PluginServer.h>

#include "ProtoShopProtocol.h"

/* Forward declarations: */
namespace Collab {
class MessageContinuation;
}
namespace MD {
class Protein;
}

namespace Collab {

namespace Plugins {

class ProtoShopServer:public PluginServer,public ProtoShopProtocol
	{
	/* Embedded classes: */
	private:
	struct ChainState // Structure defining the state of a polypeptide chain
		{
		/* Elements: */
		public:
		std::string name; // Chain's name
		Transformation transform; // Chain's transformation
		RenderFlags renderFlags; // Chain's rendering flags
		MD::Protein* protein; // Protein structure representing the chain
		};
	
	struct LockState // Structure defining the state of an active residue sequence lock
		{
		/* Elements: */
		public:
		ClientID clientId; // ID of the client holding the lock
		unsigned int chainIndex; // Index of the polypeptide chain to which the lock applies
		unsigned int beginIndex; // Index of the first locked residue
		unsigned int endIndex; // Index one after the last locked residue
		};
	
	class Client:public PluginServer::Client // Class representing a client participating in the ProtoShop protocol
		{
		friend class ProtoShopServer;
		
		/* Elements: */
		private:
		Misc::SimpleSet<LockID> lockIds; // List of IDs of locks held by the client
		
		/* Constructors and destructors: */
		public:
		Client(void);
		virtual ~Client(void);
		};
	
	typedef Misc::HashTable<LockID,LockState> LockMap; // Hash table from lock IDs to lock states
	
	/* Elements: */
	std::vector<ChainState> chains; // List of polypeptide chains
	LockID lastLockId; // ID assigned to the most recently created lock
	LockMap locks; // Map from lock IDs to active lock states
	
	/* Private methods: */
	void sendMessage(unsigned int clientId,bool broadcast,unsigned int messageId,const void* messageStructure); // Sends a message from the server to the given or all but the given clients
	
	/* Message handling methods: */
	MessageContinuation* lockRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation);
	MessageContinuation* lockUpdateRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation);
	MessageContinuation* ikUpdateRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation);
	MessageContinuation* transformRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation);
	MessageContinuation* unlockRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation);
	MessageContinuation* resetTransformRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation);
	MessageContinuation* renderFlagsUpdateRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation);
	
	/* Command interface methods: */
	void loadFileCommandCallback(const char* argumentBegin,const char* argumentEnd);
	
	/* Constructors and destructors: */
	public:
	ProtoShopServer(Server* sServer);
	virtual ~ProtoShopServer(void);
	
	/* Methods from class PluginServer: */
	virtual const char* getName(void) const;
	virtual unsigned int getVersion(void) const;
	virtual unsigned int getNumClientMessages(void) const;
	virtual unsigned int getNumServerMessages(void) const;
	virtual void setMessageBases(unsigned int newClientMessageBase,unsigned int newServerMessageBase);
	virtual void start(void);
	virtual void clientConnected(unsigned int clientId);
	virtual void clientDisconnected(unsigned int clientId);
	};

}

}

#endif
