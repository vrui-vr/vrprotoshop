/***********************************************************************
ProtoShopServer - Class representing a server for a ProtoShop workspace.
Copyright (c) 2003-2024 Oliver Kreylos

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

#include "ProtoShopServer.h"

#include <string.h>
#include <stdexcept>
#include <Misc/StdError.h>
#include <Collaboration2/Protocol.h>
#include <Collaboration2/DataType.icpp>
#include <Collaboration2/MessageWriter.h>
#include <Collaboration2/MessageContinuation.h>
#include <Collaboration2/NonBlockSocket.h>
#include <Collaboration2/Server.h>

#include "Protein.h"
#include "ParsePdbFile.h"

namespace Collab {

namespace Plugins {

/****************************************
Methods of class ProtoShopServer::Client:
****************************************/

ProtoShopServer::Client::Client(void)
	{
	}

ProtoShopServer::Client::~Client(void)
	{
	}

/********************************
Methods of class ProtoShopServer:
********************************/

void ProtoShopServer::sendMessage(unsigned int clientId,bool broadcast,unsigned int messageId,const void* messageStructure)
	{
	/* Create a message writer: */
	MessageWriter message(MessageBuffer::create(serverMessageBase+messageId,protocolTypes.calcSize(serverMessageTypes[messageId],messageStructure)));
	
	/* Write the message structure into the message: */
	protocolTypes.write(serverMessageTypes[messageId],messageStructure,message);
	
	/* Broadcast or send the message: */
	if(broadcast)
		broadcastMessage(clientId,message.getBuffer());
	else
		server->queueMessage(clientId,message.getBuffer());
	}

MessageContinuation* ProtoShopServer::lockRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation)
	{
	/* Access the base client state object, its TCP socket, and the ProtoShop client object: */
	Server::Client* client=server->getClient(clientId);
	NonBlockSocket& socket=client->getSocket();
	Client* psClient=client->getPlugin<Client>(pluginIndex);
	
	/* Read the request message: */
	LockRequestMsg lockRequest;
	protocolTypes.read(socket,clientMessageTypes[LockRequest],&lockRequest);
	Lock& lock=lockRequest.lock;
	
	/* Check if the lock request is valid: */
	if(lock.chain>=chains.size())
		throw std::runtime_error("ProtoShopServer::lockRequestCallback: Chain index out of bounds");
	if(lock.begin>=lock.end)
		throw std::runtime_error("ProtoShopServer::lockRequestCallback: Invalid residue sequence");
	if(lock.end>(unsigned int)(chains[lock.chain].protein->getNumResidues()))
		throw std::runtime_error("ProtoShopServer::lockRequestCallback: Residue sequence out of bounds");
	
	/* Prepare a reply message: */
	LockReplyMsg lockReply;
	lockReply.requestId=lockRequest.requestId;
	lockReply.lockId=LockID(0);
	
	/* Check if the lock can be granted: */
	bool grantLock=true;
	for(LockMap::Iterator lIt=locks.begin();lIt!=locks.end()&&grantLock;++lIt)
		{
		/* Check if the lock overlaps the requested residue sequence: */
		LockState& ls=lIt->getDest();
		grantLock=ls.chainIndex!=lock.chain||ls.beginIndex>=lock.end||lock.begin>=ls.endIndex;
		}
	
	if(grantLock)
		{
		/* Create a new lock and add it to the lock map: */
		do
			{
			++lastLockId;
			}
		while(lastLockId==0||locks.isEntry(lastLockId));
		LockState newLock;
		newLock.clientId=clientId;
		newLock.chainIndex=lock.chain;
		newLock.beginIndex=lock.begin;
		newLock.endIndex=lock.end;
		locks.setEntry(LockMap::Entry(lastLockId,newLock));
		
		/* Add the new lock to the client's lock list: */
		psClient->lockIds.add(lastLockId);
		
		/* Send the new lock's ID to the requesting client: */
		lockReply.lockId=lastLockId;
		
		/* Notify all other clients of the new lock: */
		LockNotificationMsg lockNotification;
		lockNotification.lockId=lastLockId;
		lockNotification.clientId=ClientID(clientId);
		lockNotification.lock=lock;
		sendMessage(clientId,true,LockNotification,&lockNotification);
		}
	
	/* Send the reply message to the requesting client: */
	sendMessage(clientId,false,LockReply,&lockReply);
	
	/* Done with message: */
	return 0;
	}

MessageContinuation* ProtoShopServer::lockUpdateRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation)
	{
	/* Access the base client state object and its TCP socket: */
	Server::Client* client=server->getClient(clientId);
	NonBlockSocket& socket=client->getSocket();
	
	/* Read the request message: */
	LockUpdateRequestMsg lockUpdateRequest;
	protocolTypes.read(socket,clientMessageTypes[LockUpdateRequest],&lockUpdateRequest);
	Lock& lock=lockUpdateRequest.lock;
	
	/* Check if the lock request is valid: */
	if(lock.chain>=chains.size())
		throw std::runtime_error("ProtoShopServer::lockUpdateRequestCallback: Chain index out of bounds");
	if(lock.begin>=lock.end)
		throw std::runtime_error("ProtoShopServer::lockUpdateRequestCallback: Invalid residue sequence");
	if(lock.end>(unsigned int)(chains[lock.chain].protein->getNumResidues()))
		throw std::runtime_error("ProtoShopServer::lockUpdateRequestCallback: Residue sequence out of bounds");
	
	/* Check if the client actually holds the claimed lock: */
	LockState& lockState=locks.getEntry(lockUpdateRequest.lockId).getDest();
	if(lockState.clientId==clientId)
		{
		/* Prepare a reply message: */
		LockReplyMsg lockReply;
		lockReply.requestId=lockUpdateRequest.requestId;
		lockReply.lockId=LockID(0);
		
		/* Check if the lock can be granted: */
		bool grantLock=lockState.chainIndex==lock.chain;
		for(LockMap::Iterator lIt=locks.begin();lIt!=locks.end()&&grantLock;++lIt)
			if(lIt->getSource()!=lockUpdateRequest.lockId)
				{
				/* Check if the lock overlaps the requested residue sequence: */
				LockState& ls=lIt->getDest();
				grantLock=ls.chainIndex!=lock.chain||ls.beginIndex>=lock.end||lock.begin>=ls.endIndex;
				}
		
		if(grantLock)
			{
			/* Update the existing lock: */
			lockState.beginIndex=lock.begin;
			lockState.endIndex=lock.end;
			
			/* Send the updated lock's ID to the requesting client: */
			lockReply.lockId=lockUpdateRequest.lockId;
			
			/* Notify all other clients of the updated lock: */
			LockUpdateNotificationMsg lockUpdateNotification;
			lockUpdateNotification.lockId=lockUpdateRequest.lockId;
			lockUpdateNotification.lock=lock;
			sendMessage(clientId,true,LockUpdateNotification,&lockUpdateNotification);
			}
		
		/* Send the reply message to the requesting client: */
		sendMessage(clientId,false,LockReply,&lockReply);
		}
	else
		throw std::runtime_error("ProtoShopServer::lockUpdateRequestCallback: Invalid lock ID");
	
	/* Done with message: */
	return 0;
	}

MessageContinuation* ProtoShopServer::ikUpdateRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation)
	{
	/* Access the base client state object and its TCP socket: */
	Server::Client* client=server->getClient(clientId);
	NonBlockSocket& socket=client->getSocket();
	
	/* Check if this is the start of a new message: */
	if(continuation==0)
		{
		/* Prepare to read the IK update request message: */
		continuation=protocolTypes.prepareReading(clientMessageTypes[IKUpdateRequest],new IKUpdateMsg);
		}
	
	/* Continue reading the connect reply message and check whether it's complete: */
	if(protocolTypes.continueReading(socket,continuation))
		{
		/* Retrieve the IK update request message: */
		IKUpdateMsg* msg=protocolTypes.getReadObject<IKUpdateMsg>(continuation);
		
		/* Check if the client actually holds the claimed lock: */
		LockState& lockState=locks.getEntry(msg->lockId).getDest();
		if(lockState.clientId==clientId)
			{
			/* Apply the IK update package to the affected polypeptide chain: */
			unsigned int numResidues=lockState.endIndex-lockState.beginIndex;
			if(msg->dihedralAngles.size()!=numResidues)
				throw std::runtime_error("ProtoShopServer::ikUpdateRequestCallback: Mismatching dihedral angle sequence");
			chains[lockState.chainIndex].protein->changeDihedralAngles(lockState.beginIndex,numResidues,msg->dihedralAngles.data(),msg->updateDirection,msg->limitUpdate);
			
			/* Notify all other clients of the IK update: */
			sendMessage(clientId,true,IKUpdateNotification,msg);
			}
		else
			throw std::runtime_error("ProtoShopServer::ikUpdateRequestCallback: Invalid lock ID");
		
		/* Delete the IK update request message and the continuation object: */
		delete msg;
		delete continuation;
		continuation=0;
		}
	
	return continuation;
	}

MessageContinuation* ProtoShopServer::transformRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation)
	{
	/* Access the base client state object and its TCP socket: */
	Server::Client* client=server->getClient(clientId);
	NonBlockSocket& socket=client->getSocket();
	
	/* Read the request message: */
	TransformMsg transformRequest;
	protocolTypes.read(socket,clientMessageTypes[TransformRequest],&transformRequest);
	
	/* Check if the client actually holds the claimed lock: */
	LockState& lockState=locks.getEntry(transformRequest.lockId).getDest();
	if(lockState.clientId==clientId)
		{
		/* Check if the transform request is valid: */
		if(lockState.chainIndex>=chains.size())
			throw std::runtime_error("ProtoShopServer::transformRequestCallback: Chain index out of bounds");
		
		/* Apply the transformation to the affected polypeptide chain: */
		chains[lockState.chainIndex].transform=transformRequest.transformation;
		
		/* Notify all other clients of the transformation: */
		sendMessage(clientId,true,TransformNotification,&transformRequest);
		}
	else
		throw std::runtime_error("ProtoShopServer::transformRequestCallback: Invalid lock ID");
	
	/* Done with message: */
	return 0;
	}

MessageContinuation* ProtoShopServer::unlockRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation)
	{
	/* Access the base client state object, its TCP socket, and the ProtoShop client object: */
	Server::Client* client=server->getClient(clientId);
	NonBlockSocket& socket=client->getSocket();
	Client* psClient=client->getPlugin<Client>(pluginIndex);
	
	/* Read the request message: */
	UnlockMsg unlockRequest;
	protocolTypes.read(socket,clientMessageTypes[UnlockRequest],&unlockRequest);
	
	/* Check if the client actually holds the claimed lock: */
	LockMap::Iterator lockIt=locks.findEntry(unlockRequest.lockId);
	if(!lockIt.isFinished()&&lockIt->getDest().clientId==clientId)
		{
		/* Remove the lock from the lock map: */
		locks.removeEntry(lockIt);
		
		/* Remove the lock from the client's lock list: */
		psClient->lockIds.remove(unlockRequest.lockId);
		
		/* Notify all other clients of the lock release: */
		sendMessage(clientId,true,UnlockNotification,&unlockRequest);
		}
	else
		throw std::runtime_error("ProtoShopServer::unlockRequestCallback: Invalid lock ID");
	
	/* Done with message: */
	return 0;
	}

MessageContinuation* ProtoShopServer::resetTransformRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation)
	{
	/* Access the base client state object and its TCP socket: */
	Server::Client* client=server->getClient(clientId);
	NonBlockSocket& socket=client->getSocket();
	
	/* Read the request message: */
	ResetTransformMsg resetTransformRequest;
	protocolTypes.read(socket,clientMessageTypes[ResetTransformRequest],&resetTransformRequest);
	
	/* Check if the reset transform request is valid: */
	if(resetTransformRequest.chain>=chains.size())
		throw std::runtime_error("ProtoShopServer::resetTransformRequestCallback: Chain index out of bounds");
	
	/* Check if the request can be granted: */
	bool grantLock=true;
	for(LockMap::Iterator lIt=locks.begin();lIt!=locks.end()&&grantLock;++lIt)
		{
		/* Check if the lock affects the chain to be reset: */
		grantLock=lIt->getDest().chainIndex!=resetTransformRequest.chain;
		}
	
	if(grantLock)
		{
		/* Reset the chain's transformation: */
		chains[resetTransformRequest.chain].transform=Transformation::identity;
		
		/* Notify all clients, including the requesting one, of the reset: */
		sendMessage(0,true,ResetTransformNotification,&resetTransformRequest);
		}
	
	/* Done with message: */
	return 0;
	}

MessageContinuation* ProtoShopServer::renderFlagsUpdateRequestCallback(unsigned int messageId,unsigned int clientId,MessageContinuation* continuation)
	{
	/* Access the base client state object and its TCP socket: */
	Server::Client* client=server->getClient(clientId);
	NonBlockSocket& socket=client->getSocket();
	
	/* Read the request message: */
	RenderFlagsUpdateMsg renderFlagsUpdateRequest;
	protocolTypes.read(socket,clientMessageTypes[RenderFlagsUpdateRequest],&renderFlagsUpdateRequest);
	
	/* Check if the update request is valid: */
	if(renderFlagsUpdateRequest.chain>=chains.size())
		throw std::runtime_error("ProtoShopServer::renderFlagsUpdateRequestCallback: Chain index out of bounds");
	
	/* Apply the update request: */
	chains[renderFlagsUpdateRequest.chain].renderFlags=renderFlagsUpdateRequest.renderFlags;
	
	/* Notify all other clients of the update: */
	sendMessage(clientId,true,RenderFlagsUpdateNotification,&renderFlagsUpdateRequest);
	
	/* Done with message: */
	return 0;
	}

void ProtoShopServer::loadFileCommandCallback(const char* argumentBegin,const char* argumentEnd)
	{
	/* Check if the requested file is a PDB file or a prediction file: */
	std::string filePath(argumentBegin,argumentEnd);
	const char* fnPtr=filePath.c_str();
	const char* extPtr=0;
	for(const char* aPtr=filePath.c_str();*aPtr!='\0';++aPtr)
		{
		if(*aPtr=='/')
			{
			fnPtr=aPtr+1;
			extPtr=0;
			}
		else if(*aPtr=='.')
			extPtr=aPtr+1;
		}
	if(extPtr==0)
		throw std::runtime_error("Missing file name extension");
	else if(strcmp(extPtr,"pdb")==0)
		{
		/* Load a PDB file: */
		std::vector<MD::Protein*> proteins=MD::parsePdbFile(filePath.c_str());
		
		/* Create a set of polypeptide chains: */
		chains.clear();
		chains.reserve(proteins.size());
		for(std::vector<MD::Protein*>::iterator pIt=proteins.begin();pIt!=proteins.end();++pIt)
			{
			/* Create a new chain: */
			chains.push_back(ChainState());
			ChainState& chain=chains.back();
			
			/* Generate a name for the chain: */
			chain.name=fnPtr;
			if(proteins.size()>1)
				{
				chain.name.push_back(' ');
				chain.name.push_back('A'+(pIt-proteins.begin()));
				}
			
			/* Initialize the chain transformation: */
			chain.transform=Transformation::identity;
			
			/* Initialize the chain's rendering flags: */
			chain.renderFlags.visible=true;
			chain.renderFlags.drawAtoms=false;
			chain.renderFlags.drawBonds=false;
			chain.renderFlags.drawBackbone=true;
			chain.renderFlags.drawCartoon=false;
			chain.renderFlags.drawHydrogenBonds=true;
			chain.renderFlags.drawHydrogenBondSites=true;
			chain.renderFlags.drawHydrogenCages=false;
			chain.renderFlags.drawAtomCollisions=false;
			
			/* Store the chain's protein: */
			chain.protein=*pIt;
			}
		}
	else
		throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Unrecognized file name extension .%s",extPtr);
	}

ProtoShopServer::ProtoShopServer(Server* sServer)
	:PluginServer(sServer),
	 lastLockId(0),locks(17)
	{
	/* Register pipe commands: */
	server->getCommandDispatcher().addCommandCallback("ProtoShop::loadFile",Misc::CommandDispatcher::wrapMethod<ProtoShopServer,&ProtoShopServer::loadFileCommandCallback>,this,"<file name>","Loads a PDB or prediction file");
	}

ProtoShopServer::~ProtoShopServer(void)
	{
	/* Delete all polypeptide chains: */
	for(std::vector<ChainState>::iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		delete cIt->protein;
	
	/* Unregister pipe commands: */
	server->getCommandDispatcher().removeCommandCallback("ProtoShop::loadFile");
	}

const char* ProtoShopServer::getName(void) const
	{
	return protocolName;
	}

unsigned int ProtoShopServer::getVersion(void) const
	{
	return protocolVersion;
	}

unsigned int ProtoShopServer::getNumClientMessages(void) const
	{
	return NumClientMessages;
	}

unsigned int ProtoShopServer::getNumServerMessages(void) const
	{
	return NumServerMessages;
	}

void ProtoShopServer::setMessageBases(unsigned int newClientMessageBase,unsigned int newServerMessageBase)
	{
	/* Call the base class method: */
	PluginServer::setMessageBases(newClientMessageBase,newServerMessageBase);
	
	/* Register message handlers: */
	server->setMessageHandler(clientMessageBase+LockRequest,Server::wrapMethod<ProtoShopServer,&ProtoShopServer::lockRequestCallback>,this,getClientMsgSize(LockRequest));
	server->setMessageHandler(clientMessageBase+LockUpdateRequest,Server::wrapMethod<ProtoShopServer,&ProtoShopServer::lockUpdateRequestCallback>,this,getClientMsgSize(LockUpdateRequest));
	server->setMessageHandler(clientMessageBase+IKUpdateRequest,Server::wrapMethod<ProtoShopServer,&ProtoShopServer::ikUpdateRequestCallback>,this,getClientMsgSize(IKUpdateRequest));
	server->setMessageHandler(clientMessageBase+TransformRequest,Server::wrapMethod<ProtoShopServer,&ProtoShopServer::transformRequestCallback>,this,getClientMsgSize(TransformRequest));
	server->setMessageHandler(clientMessageBase+UnlockRequest,Server::wrapMethod<ProtoShopServer,&ProtoShopServer::unlockRequestCallback>,this,getClientMsgSize(UnlockRequest));
	server->setMessageHandler(clientMessageBase+ResetTransformRequest,Server::wrapMethod<ProtoShopServer,&ProtoShopServer::resetTransformRequestCallback>,this,getClientMsgSize(ResetTransformRequest));
	server->setMessageHandler(clientMessageBase+RenderFlagsUpdateRequest,Server::wrapMethod<ProtoShopServer,&ProtoShopServer::renderFlagsUpdateRequestCallback>,this,getClientMsgSize(RenderFlagsUpdateRequest));
	}

void ProtoShopServer::start(void)
	{
	}

void ProtoShopServer::clientConnected(unsigned int clientId)
	{
	/* Add the new client to the list of clients using this protocol: */
	addClientToList(clientId);
	
	/* Associate a client structure with the new client: */
	Server::Client* client=server->getClient(clientId);
	client->setPlugin(pluginIndex,new Client);
	
	/* Send the current states of all polypeptide chains to the new client: */
	ConnectReplyMsg connectReply;
	connectReply.chains.reserve(chains.size());
	for(std::vector<ChainState>::iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		{
		/* Add a new polypeptide chain: */
		connectReply.chains.push_back(Chain());
		Chain& chain=connectReply.chains.back();
		
		/* Send chain state: */
		chain.name=cIt->name;
		chain.transform=cIt->transform;
		chain.renderFlags=cIt->renderFlags;
		
		/* Send the chain's secondary structures: */
		chain.secondaryStructures.reserve(cIt->protein->getNumStructures());
		for(MD::Protein::ConstStructureIterator sIt=cIt->protein->structuresBegin();sIt!=cIt->protein->structuresEnd();++sIt)
			{
			/* Add a new secondary structure: */
			chain.secondaryStructures.push_back(SecondaryStructure());
			SecondaryStructure& secondaryStructure=chain.secondaryStructures.back();
			
			/* Send secondary structure state: */
			secondaryStructure.type=Misc::SInt8(sIt->getStructureType());
			
			/* Send the secondary structure's residues: */
			secondaryStructure.residues.reserve(sIt.getNumResidues());
			for(MD::Protein::ConstResidueIterator rIt=sIt.residuesBegin();rIt!=sIt.residuesEnd();++rIt)
				{
				/* Add a new residue: */
				secondaryStructure.residues.push_back(Residue());
				Residue& residue=secondaryStructure.residues.back();
				
				/* Send residue state: */
				residue.index=Misc::UInt16(rIt->getPdbResidueIndex());
				memcpy(residue.name,rIt->getPdbResidueName(),4);
				
				/* Send the residue's atoms: */
				for(MD::Protein::ConstAtomIterator aIt=rIt.atomsBegin();aIt!=rIt.atomsEnd();++aIt)
					{
					/* Add a new atom: */
					residue.atoms.push_back(Atom());
					Atom& atom=residue.atoms.back();
					
					/* Send atom state: */
					atom.type=Misc::UInt8(aIt->getType());
					atom.index=Misc::UInt16(aIt->getAtomIndex());
					memcpy(atom.placement,aIt->getPlacement(),4);
					atom.position=aIt->getPosition();
					}
				}
			}
		}
	sendMessage(clientId,false,ConnectReply,&connectReply);
	
	/* Send the states of all active residue sequence locks to the new client: */
	for(LockMap::Iterator lIt=locks.begin();!lIt.isFinished();++lIt)
		{
		LockNotificationMsg lockNotification;
		lockNotification.lockId=lIt->getSource();
		lockNotification.clientId=lIt->getDest().clientId;
		lockNotification.lock.chain=lIt->getDest().chainIndex;
		lockNotification.lock.begin=lIt->getDest().beginIndex;
		lockNotification.lock.end=lIt->getDest().endIndex;
		sendMessage(clientId,false,LockNotification,&lockNotification);
		}
	}

void ProtoShopServer::clientDisconnected(unsigned int clientId)
	{
	/* Release the client's remaining active locks: */
	Client* psClient=server->getPlugin<Client>(clientId,pluginIndex);
	for(Misc::SimpleSet<LockID>::iterator clIt=psClient->lockIds.begin();clIt!=psClient->lockIds.end();++clIt)
		{
		/* Send an unlock notification to all other clients: */
		UnlockMsg unlockNotification;
		unlockNotification.lockId=*clIt;
		sendMessage(clientId,true,UnlockNotification,&unlockNotification);
		
		/* Remove the lock: */
		locks.removeEntry(*clIt);
		}
	
	/* Remove the client from the list of clients using this protocol: */
	removeClientFromList(clientId);
	}

/***********************
DSO loader entry points:
***********************/

extern "C" {

PluginServer* createObject(PluginServerLoader& objectLoader,Server* server)
	{
	return new ProtoShopServer(server);
	}

void destroyObject(PluginServer* object)
	{
	delete object;
	}

}

}

}
