/***********************************************************************
ProtoShopClient - Class representing a client for a ProtoShop workspace.
Copyright (c) 2003-2023 Oliver Kreylos

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

#include "ProtoShopClient.h"

#include <Collaboration2/MessageReader.h>
#include <Collaboration2/MessageContinuation.h>

#include "Protein.h"
#include "ProteinInteractor.h"
#include "VRProtoShop.h"

namespace Collab {

namespace Plugins {

/********************************
Methods of class ProtoShopClient:
********************************/

void ProtoShopClient::frontendConnectReplyCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message from the backend: */
	std::vector<ServerChain>* chains=message.read<std::vector<ServerChain>*>();
	
	/* Replace the current polypeptide chains: */
	for(VRProtoShop::ChainList::iterator cIt=application->chains.begin();cIt!=application->chains.end();++cIt)
		{
		delete cIt->interactor;
		delete cIt->proteinRenderer;
		delete cIt->protein;
		}
	application->chains.reserve(chains->size());
	unsigned int chainIndex=0;
	for(std::vector<ServerChain>::const_iterator cIt=chains->begin();cIt!=chains->end();++cIt,++chainIndex)
		{
		/* Create a new chain structure: */
		VRProtoShop::Chain newChain;
		newChain.name=cIt->name;
		newChain.protein=cIt->protein;
		newChain.proteinRenderer=new MD::ProteinRenderer(application->configFile.getSection("/ProteinRenderer"),newChain.protein);
		newChain.interactor=new ProteinInteractor(chainIndex,newChain.protein,newChain.proteinRenderer,0,application->undoBuffer,application->interactorPipe);
		newChain.interactor->collaborate(this);
		newChain.visible=cIt->renderFlags.visible;
		
		/* Set the new chain's initial rendering flags: */
		newChain.proteinRenderer->setDrawAtoms(cIt->renderFlags.drawAtoms);
		newChain.proteinRenderer->setDrawBonds(cIt->renderFlags.drawBonds);
		newChain.proteinRenderer->setDrawBackbone(cIt->renderFlags.drawBackbone);
		newChain.proteinRenderer->setDrawCartoon(cIt->renderFlags.drawCartoon);
		newChain.proteinRenderer->setDrawHydrogenBonds(cIt->renderFlags.drawHydrogenBonds);
		newChain.proteinRenderer->setDrawHydrogenBondSites(cIt->renderFlags.drawHydrogenBondSites);
		newChain.proteinRenderer->setDrawHydrogenCages(cIt->renderFlags.drawHydrogenCages);
		newChain.proteinRenderer->setDrawCollisions(cIt->renderFlags.drawAtomCollisions);
		
		/* Store the new chain: */
		application->chains.push_back(newChain);
		}
	
	/* Update the chain dialog: */
	application->activeChain=&application->chains.front();
	application->updateChainDialog();
	
	/* Reset the navigation transformation: */
	application->resetNavigation();
	
	/* Delete the list of chains: */
	delete chains;
	}

MessageContinuation* ProtoShopClient::connectReplyCallback(unsigned int messageId,MessageContinuation* continuation)
	{
	/* Check if this is the start of a new message: */
	if(continuation==0)
		{
		/* Prepare to read the connect reply message: */
		continuation=protocolTypes.prepareReading(serverMessageTypes[ConnectReply],new ConnectReplyMsg);
		}
	
	/* Continue reading the connect reply message and check whether it's complete: */
	if(protocolTypes.continueReading(client->getSocket(),continuation))
		{
		/* Extract the connect reply message: */
		ConnectReplyMsg* msg=protocolTypes.getReadObject<ConnectReplyMsg>(continuation);
		
		/* Parse the connect reply message: */
		std::vector<ServerChain>* chains=new std::vector<ServerChain>;
		chains->reserve(msg->chains.size());
		for(Misc::Vector<Chain>::iterator cIt=msg->chains.begin();cIt!=msg->chains.end();++cIt)
			{
			/* Parse the polypeptide chain: */
			chains->push_back(ServerChain());
			ServerChain& chain=chains->back();
			chain.name=cIt->name;
			chain.transform=cIt->transform;
			chain.renderFlags=cIt->renderFlags;
			chain.protein=new MD::Protein;
			MD::Protein::ResidueCreator residueCreator(chain.protein);
			
			/* Create the chain's secondary structures: */
			for(Misc::Vector<SecondaryStructure>::iterator ssIt=cIt->secondaryStructures.begin();ssIt!=cIt->secondaryStructures.end();++ssIt)
				{
				/* Create the secondary structure: */
				residueCreator.newSecondaryStructure(MD::Protein::SecondaryStructure::StructureType(ssIt->type));
				
				/* Create the secondary structure's residues: */
				for(Misc::Vector<Residue>::iterator rIt=ssIt->residues.begin();rIt!=ssIt->residues.end();++rIt)
					{
					/* Create the residue: */
					rIt->name[3]='\0';
					residueCreator.newResidue(rIt->name,rIt->index);
					
					/* Create the residue's atoms: */
					for(Misc::Vector<Atom>::iterator aIt=rIt->atoms.begin();aIt!=rIt->atoms.end();++aIt)
						{
						/* Add the atom to the protein: */
						aIt->placement[3]='\0';
						residueCreator.addAtom(MD::Atom::getElementName(MD::Atom::Element(aIt->type)),aIt->index,aIt->position,aIt->placement);
						}
					}
				}
			
			/* Finalize the protein: */
			residueCreator.finishProtein();
			}
		
		/* Forward the list of polypeptide chains to the frontend: */
		{
		MessageWriter frontendMessage(MessageBuffer::create(serverMessageBase+ConnectReply,sizeof(std::vector<ServerChain>*)));
		frontendMessage.write(chains);
		client->queueFrontendMessage(frontendMessage.getBuffer());
		}
		
		/* Delete the connect reply message and the continuation object: */
		delete msg;
		delete continuation;
		continuation=0;
		}
	
	return continuation;
	}

void ProtoShopClient::frontendLockReplyCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	LockReplyMsg lockReply;
	protocolTypes.read(message,serverMessageTypes[LockReply],&lockReply);
	
	/* Find the outstanding lock: */
	LockRequestMap::Iterator lrIt=lockRequests.findEntry(lockReply.requestId);
	if(!lrIt.isFinished())
		{
		/* Check if this is a chain dragging lock: */
		LockRequestState& lr=lrIt->getDest();
		if(lr.chainDraggingTool!=0)
			{
			/* Notify the affected chain dragging tool: */
			lr.chainDraggingTool->lockReply(lockReply.lockId);
			}
		else
			{
			/* Notify the affected protein interactor: */
			application->chains[lr.chainIndex].interactor->lockReply(lr.beginIndex,lr.endIndex-lr.beginIndex,lockReply.lockId);
			}
		
		/* Check if the lock request was granted: */
		if(lockReply.lockId!=LockID(0))
			{
			/* Create a lock state for the local client: */
			LockState ls;
			ls.clientId=0;
			ls.chainIndex=lr.chainIndex;
			ls.beginIndex=lr.beginIndex;
			ls.endIndex=lr.endIndex;
			locks.setEntry(LockMap::Entry(lockReply.lockId,ls));
			}
		
		/* Remove the outstanding lock request: */
		lockRequests.removeEntry(lrIt);
		}
	else
		{
		/* Assume the lock request was cancelled, and release a granted lock immediately: */
		if(lockReply.lockId!=LockID(0))
			{
			UnlockMsg unlockRequest;
			unlockRequest.lockId=lockReply.lockId;
			queueServerMessage(UnlockRequest,&unlockRequest);
			}
		}
	}

void ProtoShopClient::frontendLockNotificationCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	LockNotificationMsg lockNotification;
	protocolTypes.read(message,serverMessageTypes[LockNotification],&lockNotification);
	
	/* Enter the new lock into the lock map: */
	LockState ls;
	ls.clientId=lockNotification.clientId;
	ls.chainIndex=lockNotification.lock.chain;
	ls.beginIndex=lockNotification.lock.begin;
	ls.endIndex=lockNotification.lock.end;
	locks.setEntry(LockMap::Entry(lockNotification.lockId,ls));
	
	/* Update the affected chain's visual state: */
	application->chains[ls.chainIndex].proteinRenderer->lockResidueRange(ls.beginIndex,ls.endIndex-ls.beginIndex,true);
	}

void ProtoShopClient::frontendLockUpdateNotificationCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	LockUpdateNotificationMsg lockUpdateNotification;
	protocolTypes.read(message,serverMessageTypes[LockUpdateNotification],&lockUpdateNotification);
	
	/* Find the existing lock in the lock map: */
	LockState& ls=locks.getEntry(lockUpdateNotification.lockId).getDest();
	
	/* Update the affected chain's visual state and the lock state: */
	application->chains[ls.chainIndex].proteinRenderer->lockResidueRange(ls.beginIndex,ls.endIndex-ls.beginIndex,false);
	ls.beginIndex=lockUpdateNotification.lock.begin;
	ls.endIndex=lockUpdateNotification.lock.end;
	application->chains[ls.chainIndex].proteinRenderer->lockResidueRange(ls.beginIndex,ls.endIndex-ls.beginIndex,true);
	}

void ProtoShopClient::frontendIKUpdateNotificationCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	IKUpdateMsg ikUpdateNotification;
	protocolTypes.read(message,serverMessageTypes[IKUpdateNotification],&ikUpdateNotification);
	
	/* Find the existing lock in the lock map: */
	LockState& ls=locks.getEntry(ikUpdateNotification.lockId).getDest();
	
	/* Update the affected chain: */
	VRProtoShop::Chain& chain=application->chains[ls.chainIndex];
	{
	Threads::Mutex::Lock proteinLock(chain.proteinMutex);
	chain.protein->changeDihedralAngles(ls.beginIndex,ls.endIndex-ls.beginIndex,ikUpdateNotification.dihedralAngles.data(),ikUpdateNotification.updateDirection,ikUpdateNotification.limitUpdate);
	chain.proteinRenderer->updateProtein();
	}
	}

MessageContinuation* ProtoShopClient::ikUpdateNotificationCallback(unsigned int messageId,MessageContinuation* continuation)
	{
	/* Check if this is the start of a new message: */
	if(continuation==0)
		{
		/* Prepare to read the IK update notification message: */
		continuation=protocolTypes.prepareReading(serverMessageTypes[IKUpdateNotification],new IKUpdateMsg);
		}
	
	/* Continue reading the IK update notification message and check whether it's complete: */
	if(protocolTypes.continueReading(client->getSocket(),continuation))
		{
		/* Extract the IK update notification message: */
		IKUpdateMsg* msg=protocolTypes.getReadObject<IKUpdateMsg>(continuation);
		
		/* Forward the IK update notification message to the front end: */
		{
		MessageWriter frontendMessage(MessageBuffer::create(serverMessageBase+IKUpdateNotification,protocolTypes.calcSize(serverMessageTypes[IKUpdateNotification],msg)));
		protocolTypes.write(serverMessageTypes[IKUpdateNotification],msg,frontendMessage);
		client->queueFrontendMessage(frontendMessage.getBuffer());
		}
		
		/* Delete the IK update notification message and the continuation object: */
		delete msg;
		delete continuation;
		continuation=0;
		}
	
	return continuation;
	}

void ProtoShopClient::frontendTransformNotificationCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	TransformMsg transformNotification;
	protocolTypes.read(message,serverMessageTypes[TransformNotification],&transformNotification);
	
	/* Find the existing lock in the lock map: */
	LockState& ls=locks.getEntry(transformNotification.lockId).getDest();
	
	/* Update the affected chain's transformation: */
	application->chains[ls.chainIndex].transformation=transformNotification.transformation;
	}

void ProtoShopClient::frontendUnlockNotificationCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	UnlockMsg unlockNotification;
	protocolTypes.read(message,serverMessageTypes[UnlockNotification],&unlockNotification);
	
	/* Find the existing lock in the lock map: */
	LockMap::Iterator lIt=locks.findEntry(unlockNotification.lockId);
	
	/* Update the affected chain's visual state and release the lock: */
	LockState& ls=lIt->getDest();
	application->chains[ls.chainIndex].proteinRenderer->lockResidueRange(ls.beginIndex,ls.endIndex-ls.beginIndex,false);
	locks.removeEntry(lIt);
	}

void ProtoShopClient::frontendResetTransformNotificationCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	ResetTransformMsg resetTransformNotification;
	protocolTypes.read(message,serverMessageTypes[ResetTransformNotification],&resetTransformNotification);
	
	/* Reset the affected chain's transformation: */
	application->chains[resetTransformNotification.chain].transformation=Transformation::identity;
	}

void ProtoShopClient::frontendRenderFlagsUpdateNotificationCallback(unsigned int messageId,MessageReader& message)
	{
	/* Read the message: */
	RenderFlagsUpdateMsg renderFlagsUpdateNotification;
	protocolTypes.read(message,serverMessageTypes[RenderFlagsUpdateNotification],&renderFlagsUpdateNotification);
	
	/* Update the affected chain's rendering flags: */
	VRProtoShop::Chain& chain=application->chains[renderFlagsUpdateNotification.chain];
	chain.visible=renderFlagsUpdateNotification.renderFlags.visible;
	chain.proteinRenderer->setDrawAtoms(renderFlagsUpdateNotification.renderFlags.drawAtoms);
	chain.proteinRenderer->setDrawBonds(renderFlagsUpdateNotification.renderFlags.drawBonds);
	chain.proteinRenderer->setDrawBackbone(renderFlagsUpdateNotification.renderFlags.drawBackbone);
	chain.proteinRenderer->setDrawCartoon(renderFlagsUpdateNotification.renderFlags.drawCartoon);
	chain.proteinRenderer->setDrawHydrogenBonds(renderFlagsUpdateNotification.renderFlags.drawHydrogenBonds);
	chain.proteinRenderer->setDrawHydrogenBondSites(renderFlagsUpdateNotification.renderFlags.drawHydrogenBondSites);
	chain.proteinRenderer->setDrawHydrogenCages(renderFlagsUpdateNotification.renderFlags.drawHydrogenCages);
	chain.proteinRenderer->setDrawCollisions(renderFlagsUpdateNotification.renderFlags.drawAtomCollisions);
	
	/* Update the GUI: */
	if(application->activeChain==&chain)
		application->updateChainDialogToggles();
	}

ProtoShopClient::ProtoShopClient(Client* sClient,VRProtoShop* sApplication)
	:PluginClient(sClient),
	 application(sApplication),
	 lastLockRequestId(0),lockRequests(5),
	 locks(17)
	{
	}

ProtoShopClient::~ProtoShopClient(void)
	{
	}

const char* ProtoShopClient::getName(void) const
	{
	return protocolName;
	}

unsigned int ProtoShopClient::getVersion(void) const
	{
	return protocolVersion;
	}

unsigned int ProtoShopClient::getNumClientMessages(void) const
	{
	return NumClientMessages;
	}

unsigned int ProtoShopClient::getNumServerMessages(void) const
	{
	return NumServerMessages;
	}

void ProtoShopClient::setMessageBases(unsigned int newClientMessageBase,unsigned int newServerMessageBase)
	{
	/* Call the base class method: */
	PluginClient::setMessageBases(newClientMessageBase,newServerMessageBase);
	
	/* Register message handlers: */
	client->setTCPMessageHandler(serverMessageBase+ConnectReply,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::connectReplyCallback>,this,getServerMsgSize(ConnectReply));
	client->setFrontendMessageHandler(serverMessageBase+ConnectReply,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendConnectReplyCallback>,this);
	client->setMessageForwarder(serverMessageBase+LockReply,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendLockReplyCallback>,this,getServerMsgSize(LockReply));
	client->setMessageForwarder(serverMessageBase+LockNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendLockNotificationCallback>,this,getServerMsgSize(LockNotification));
	client->setMessageForwarder(serverMessageBase+LockUpdateNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendLockUpdateNotificationCallback>,this,getServerMsgSize(LockUpdateNotification));
	client->setTCPMessageHandler(serverMessageBase+IKUpdateNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::ikUpdateNotificationCallback>,this,getServerMsgSize(IKUpdateNotification));
	client->setFrontendMessageHandler(serverMessageBase+IKUpdateNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendIKUpdateNotificationCallback>,this);
	client->setMessageForwarder(serverMessageBase+TransformNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendTransformNotificationCallback>,this,getServerMsgSize(TransformNotification));
	client->setMessageForwarder(serverMessageBase+UnlockNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendUnlockNotificationCallback>,this,getServerMsgSize(UnlockNotification));
	client->setMessageForwarder(serverMessageBase+ResetTransformNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendResetTransformNotificationCallback>,this,getServerMsgSize(ResetTransformNotification));
	client->setMessageForwarder(serverMessageBase+RenderFlagsUpdateNotification,Client::wrapMethod<ProtoShopClient,&ProtoShopClient::frontendRenderFlagsUpdateNotificationCallback>,this,getServerMsgSize(RenderFlagsUpdateNotification));
	}

void ProtoShopClient::start(void)
	{
	}

ProtoShopProtocol::LockID ProtoShopClient::requestLock(unsigned int chainIndex,unsigned int firstResidueIndex,unsigned int numResidues)
	{
	/* Enter a new lock request into the outstanding lock map: */
	do
		{
		++lastLockRequestId;
		}
	while(lastLockRequestId==0||lockRequests.isEntry(lastLockRequestId));
	lockRequests.setEntry(LockRequestMap::Entry(lastLockRequestId,LockRequestState(chainIndex,firstResidueIndex,firstResidueIndex+numResidues)));
	
	/* Send a lock request message to the server: */
	LockRequestMsg lockRequest;
	lockRequest.requestId=lastLockRequestId;
	lockRequest.lock.chain=chainIndex;
	lockRequest.lock.begin=firstResidueIndex;
	lockRequest.lock.end=firstResidueIndex+numResidues;
	queueServerMessage(LockRequest,&lockRequest);
	
	return lastLockRequestId;
	}

ProtoShopProtocol::LockID ProtoShopClient::requestLock(unsigned int chainIndex,VRProtoShop::ChainDraggingTool* chainDraggingTool)
	{
	/* Enter a new lock request into the outstanding lock map: */
	do
		{
		++lastLockRequestId;
		}
	while(lastLockRequestId==0||lockRequests.isEntry(lastLockRequestId));
	unsigned int numResidues=application->chains[chainIndex].protein->getNumResidues();
	lockRequests.setEntry(LockRequestMap::Entry(lastLockRequestId,LockRequestState(chainIndex,0,numResidues,chainDraggingTool)));
	
	/* Send a lock request message to the server: */
	LockRequestMsg lockRequest;
	lockRequest.requestId=lastLockRequestId;
	lockRequest.lock.chain=chainIndex;
	lockRequest.lock.begin=0;
	lockRequest.lock.end=numResidues;
	queueServerMessage(LockRequest,&lockRequest);
	
	return lastLockRequestId;
	}

ProtoShopProtocol::LockID ProtoShopClient::updateLock(ProtoShopProtocol::LockID lockId,unsigned int firstResidueIndex,unsigned int numResidues)
	{
	/* Access the existing lock: */
	LockState& ls=locks.getEntry(lockId).getDest();
	
	/* Enter a new lock request into the outstanding lock map: */
	do
		{
		++lastLockRequestId;
		}
	while(lastLockRequestId==0||lockRequests.isEntry(lastLockRequestId));
	lockRequests.setEntry(LockRequestMap::Entry(lastLockRequestId,LockRequestState(ls.chainIndex,firstResidueIndex,firstResidueIndex+numResidues)));
	
	/* Send a lock update request message to the server: */
	LockUpdateRequestMsg lockUpdateRequest;
	lockUpdateRequest.requestId=lastLockRequestId;
	lockUpdateRequest.lockId=lockId;
	lockUpdateRequest.lock.chain=ls.chainIndex;
	lockUpdateRequest.lock.begin=firstResidueIndex;
	lockUpdateRequest.lock.end=firstResidueIndex+numResidues;
	queueServerMessage(LockUpdateRequest,&lockUpdateRequest);
	
	return lastLockRequestId;
	}

void ProtoShopClient::cancelLockRequest(ProtoShopProtocol::LockID lockRequestId)
	{
	/* Remove the lock request from the map; ignore if it's not in there: */
	lockRequests.removeEntry(lockRequestId);
	}

void ProtoShopClient::dragChain(ProtoShopProtocol::LockID lockId,const ProtoShopProtocol::Transformation& newTransform)
	{
	/* Check if the lock is valid: */
	if(lockId!=0)
		{
		/* Send a transform request message to the server: */
		TransformMsg transformRequest;
		transformRequest.lockId=lockId;
		transformRequest.transformation=newTransform;
		queueServerMessage(TransformRequest,&transformRequest);
		}
	}

void ProtoShopClient::ikUpdate(ProtoShopProtocol::LockID lockId,const MD::Protein::DihedralAnglePair deltaAngles[],int updateDirection,bool limitUpdate)
	{
	/* Access the lock: */
	LockState& ls=locks.getEntry(lockId).getDest();
	
	/* Send an IK update request to the server: */
	IKUpdateMsg ikUpdateRequest;
	ikUpdateRequest.lockId=lockId;
	ikUpdateRequest.updateDirection=updateDirection;
	ikUpdateRequest.limitUpdate=limitUpdate;
	unsigned int numResidues=ls.endIndex-ls.beginIndex;
	ikUpdateRequest.dihedralAngles.reserve(numResidues);
	for(unsigned int i=0;i<numResidues;++i)
		ikUpdateRequest.dihedralAngles.push_back(deltaAngles[i]);
	queueServerMessage(IKUpdateRequest,&ikUpdateRequest);
	}

void ProtoShopClient::releaseLock(ProtoShopProtocol::LockID lockId)
	{
	/* Check if the lock is valid: */
	if(lockId!=0)
		{
		/* Remove the lock from the lock map: */
		locks.removeEntry(lockId);
		
		/* Send an unlock request message to the server: */
		UnlockMsg unlockRequest;
		unlockRequest.lockId=lockId;
		queueServerMessage(UnlockRequest,&unlockRequest);
		}
	}

void ProtoShopClient::resetChainTransform(unsigned int chainIndex)
	{
	/* Send a reset transform request to the server: */
	ResetTransformMsg resetTransformRequest;
	resetTransformRequest.chain=chainIndex;
	queueServerMessage(ResetTransformRequest,&resetTransformRequest);
	}

void ProtoShopClient::updateRenderFlags(void)
	{
	/* Send a render flags update request message to the server: */
	RenderFlagsUpdateMsg renderFlagsUpdateRequest;
	renderFlagsUpdateRequest.chain=application->chainListBox->getSelectedItem();
	renderFlagsUpdateRequest.renderFlags.visible=application->activeChain->visible;
	MD::ProteinRenderer* renderer=application->activeChain->proteinRenderer;
	renderFlagsUpdateRequest.renderFlags.drawAtoms=renderer->getDrawAtoms();
	renderFlagsUpdateRequest.renderFlags.drawBonds=renderer->getDrawBonds();
	renderFlagsUpdateRequest.renderFlags.drawBackbone=renderer->getDrawBackbone();
	renderFlagsUpdateRequest.renderFlags.drawCartoon=renderer->getDrawCartoon();
	renderFlagsUpdateRequest.renderFlags.drawHydrogenBonds=renderer->getDrawHydrogenBonds();
	renderFlagsUpdateRequest.renderFlags.drawHydrogenBondSites=renderer->getDrawHydrogenBondSites();
	renderFlagsUpdateRequest.renderFlags.drawHydrogenCages=renderer->getDrawHydrogenCages();
	renderFlagsUpdateRequest.renderFlags.drawAtomCollisions=renderer->getDrawCollisions();
	queueServerMessage(RenderFlagsUpdateRequest,&renderFlagsUpdateRequest);
	}

}

}
