/***********************************************************************
VRProtoShop - Virtual reality version of ProtoShop (will hopefully be
fully functional at some point).
Copyright (c) 2002-2025 Oliver Kreylos

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

#include "VRProtoShop.h"

#include <string.h>
#include <iostream>
#include <stdexcept>
#include <Misc/StdError.h>
#include <Misc/Timer.h>
#include <Cluster/MulticastPipe.h>
#include <Math/Math.h>
#include <Geometry/Vector.h>
#include <Geometry/Box.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/OrthogonalTransformation.h>
#include <GL/gl.h>
#include <GL/GLColorMap.h>
#include <GL/GLContextData.h>
#include <GL/GLLineIlluminator.h>
#include <GL/GLTransformationWrappers.h>
#include <GLMotif/StyleSheet.h>
#include <GLMotif/WidgetManager.h>
#include <GLMotif/PopupMenu.h>
#include <GLMotif/PopupWindow.h>
#include <GLMotif/Margin.h>
#include <GLMotif/RowColumn.h>
#include <GLMotif/Separator.h>
#include <GLMotif/Button.h>
#include <GLMotif/CascadeButton.h>
#include <GLMotif/TextField.h>
#include <Vrui/Vrui.h>
#include <Vrui/ClusterSupport.h>
#include <Collaboration2/Client.h>

#include "DragBox.h"
#include "Protein.h"
#include "ProtoShopClient.h"
#include "ProteinRenderer.h"
#include "UndoBuffer.h"
#include "ParsePdbFile.h"
#include "CreateProtein.h"
#include "Config.h"

// DEBUGGING
#include <Realtime/Time.h>

/**********************************
Methods of class VRProtoShop::Tool:
**********************************/

Vrui::ToolFactory* VRProtoShop::Tool::initClass(void)
	{
	/* Create a factory object for the base tool class: */
	ToolFactory* factory=new ToolFactory("VRProtoShopTool","VR ProtoShop",0,*Vrui::getToolManager());
	
	/* Register the custom tool class with the Vrui tool manager: */
	Vrui::getToolManager()->addAbstractClass(factory,Vrui::ToolManager::defaultToolFactoryDestructor);
	
	return factory;
	}

/*******************************************************
Static elements of class VRProtoShop::ChainDraggingTool:
*******************************************************/

VRProtoShop::ChainDraggingToolFactory* VRProtoShop::ChainDraggingTool::factory=0;

/***********************************************
Methods of class VRProtoShop::ChainDraggingTool:
***********************************************/

void VRProtoShop::ChainDraggingTool::initClass(Vrui::ToolFactory* baseClass)
	{
	Vrui::ToolManager& tm=*Vrui::getToolManager();
	
	/* Create a factory object for the tool class: */
	ChainDraggingToolFactory* factory=new ChainDraggingToolFactory("ChainDraggingTool","Drag Chains",baseClass,tm);
	
	/* Set the custom tool class' input layout: */
	factory->setNumButtons(1,false);
	factory->setButtonFunction(0,"Drag Chain");
	
	/* Register the custom tool class with the Vrui tool manager: */
	tm.addClass(factory,Vrui::ToolManager::defaultToolFactoryDestructor);
	}

VRProtoShop::ChainDraggingTool::ChainDraggingTool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment)
	:Tool(factory,inputAssignment),
	 rayBased(getButtonDevice(0)->isRayDevice()),
	 dragging(false),dragChain(0),lockRequestId(0),lockId(0)
	{
	}

VRProtoShop::ChainDraggingTool::~ChainDraggingTool(void)
	{
	/* Cancel any outstanding lock request: */
	if(lockRequestId!=0)
		application->psClient->cancelLockRequest(lockRequestId);
	}

const Vrui::ToolFactory* VRProtoShop::ChainDraggingTool::getFactory(void) const
	{
	return factory;
	}

void VRProtoShop::ChainDraggingTool::buttonCallback(int,Vrui::InputDevice::ButtonCallbackData* cbData)
	{
	if(cbData->newButtonState) // Button was just pressed
		{
		/* Bail out if still waiting for a lock reply from a previous dragging operation: */
		if(lockRequestId!=0)
			return;
		
		/* Try picking all visible polypeptide chains: */
		Vrui::NavTransform invNav=Vrui::getInverseNavigationTransformation();
		unsigned int chainIndex=0;
		if(rayBased)
			{
			/* Get the device ray in navigation coordinates: */
			Vrui::Ray deviceRay=getButtonDeviceRay(0);
			deviceRay.transform(invNav);
			deviceRay.normalizeDirection();
			
			/* Pick all chains: */
			for(ChainList::iterator cIt=application->chains.begin();cIt!=application->chains.end();++cIt,++chainIndex)
				{
				/* Transform the device ray to the chain's coordinate space: */
				Vrui::Ray chainDeviceRay=deviceRay;
				chainDeviceRay.transform(Geometry::invert(cIt->transformation));
				
				/* Pick the chain: */
				if(cIt->protein->pickAtom(chainDeviceRay)!=cIt->protein->atomsEnd())
					{
					/* Drag this chain: */
					dragChain=&*cIt;
					
					/* Stop looking: */
					break;
					}
				}
			}
		else
			{
			/* Get the device position in navigation coordinates: */
			Vrui::Point devicePos=invNav.transform(getButtonDevicePosition(0));
			
			/* Pick all chains: */
			for(ChainList::iterator cIt=application->chains.begin();cIt!=application->chains.end();++cIt,++chainIndex)
				{
				/* Transform the device position to the chain's coordinate space: */
				Vrui::Point chainDevicePos=cIt->transformation.inverseTransform(devicePos);
				
				/* Pick the chain: */
				if(cIt->protein->pickAtom(chainDevicePos)!=cIt->protein->atomsEnd())
					{
					/* Drag this chain: */
					dragChain=&*cIt;
					
					/* Stop looking: */
					break;
					}
				}
			}
		
		/* Check if the pick was successful: */
		if(dragChain!=0)
			{
			/* Start dragging the picked chain: */
			dragging=true;
			initial=dragChain->transformation;
			drag=initial;
			Vrui::NavTransform invDev0(getButtonDeviceTransformation(0));
			invDev0.leftMultiply(invNav);
			IK::Transformation invDev(invDev0.getTranslation(),invDev0.getRotation());
			invDev.doInvert();
			drag.leftMultiply(invDev);
			
			if(application->psClient!=0)
				{
				/* Request a lock for the entire picked chain: */
				lockRequestId=application->psClient->requestLock(chainIndex,this);
				}
			}
		}
	else // Button was just released
		{
		/* Stop dragging: */
		dragging=false;
		
		if(application->psClient!=0)
			{
			/* Check if there is a valid lock: */
			if(lockId!=0)
				{
				/* Release the drag chain: */
				dragChain=0;
				
				/* Release the lock: */
				application->psClient->releaseLock(lockId);
				lockId=0;
				}
			}
		else
			{
			/* Release the drag chain: */
			dragChain=0;
			}
		}
	}

void VRProtoShop::ChainDraggingTool::frame(void)
	{
	if(dragging)
		{
		/* Apply the new device transformation to the dragged chain: */
		dragChain->transformation=drag;
		Vrui::NavTransform dev(getButtonDeviceTransformation(0));
		dev.leftMultiply(Vrui::getInverseNavigationTransformation());
		dragChain->transformation.leftMultiply(IK::Transformation(dev.getTranslation(),dev.getRotation()));
		dragChain->transformation.renormalize();
		
		if(application->psClient!=0)
			{
			/* Send a dragging update to the server: */
			application->psClient->dragChain(lockId,dragChain->transformation);
			}
		}
	}

void VRProtoShop::ChainDraggingTool::lockReply(unsigned int newLockId)
	{
	/* Release the lock request: */
	lockRequestId=0;
	lockId=newLockId;
	
	/* Check if the lock was granted: */
	if(lockId!=0)
		{
		/* Send the chain's current transformation to the server: */
		application->psClient->dragChain(lockId,dragChain->transformation);
		
		/* Check if the dragging operation already finished: */
		if(!dragging)
			{
			/* Release the drag chain: */
			dragChain=0;
			
			/* Release the lock immediately: */
			application->psClient->releaseLock(lockId);
			lockId=0;
			}
		}
	else
		{
		/* Reset the dragged polypeptide chain to its original state: */
		dragChain->transformation=initial;
		
		/* Stop dragging: */
		dragging=false;
		dragChain=0;
		}
	}

/***********************************************************
Static elements of class VRProtoShop::StructureDraggingTool:
***********************************************************/

VRProtoShop::StructureDraggingToolFactory* VRProtoShop::StructureDraggingTool::factory=0;

/***************************************************
Methods of class VRProtoShop::StructureDraggingTool:
***************************************************/

void VRProtoShop::StructureDraggingTool::initClass(Vrui::ToolFactory* baseClass)
	{
	Vrui::ToolManager& tm=*Vrui::getToolManager();
	
	/* Create a factory object for the tool class: */
	StructureDraggingToolFactory* factory=new StructureDraggingToolFactory("StructureDraggingTool","Drag Secondary Structure",baseClass,tm);
	
	/* Set the custom tool class' input layout: */
	factory->setNumButtons(1,false);
	factory->setButtonFunction(0,"Drag Structure");
	
	/* Register the custom tool class with the Vrui tool manager: */
	tm.addClass(factory,Vrui::ToolManager::defaultToolFactoryDestructor);
	}

VRProtoShop::StructureDraggingTool::StructureDraggingTool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment)
	:Tool(factory,inputAssignment),
	 rayBased(getButtonDevice(0)->isRayDevice()),
	 dragChain(0)
	{
	}

const Vrui::ToolFactory* VRProtoShop::StructureDraggingTool::getFactory(void) const
	{
	return factory;
	}

void VRProtoShop::StructureDraggingTool::buttonCallback(int,Vrui::InputDevice::ButtonCallbackData* cbData)
	{
	if(cbData->newButtonState) // Button was just pressed
		{
		/* Try picking the structure dragging box: */
		if(rayBased)
			{
			/* Get the device ray in the active chain's coordinates: */
			Vrui::Ray ray=getButtonDeviceRay(0);
			ray.transform(Vrui::getInverseNavigationTransformation());
			ray.transform(Geometry::invert(application->activeChain->transformation));
			ray.normalizeDirection();
			DragBox::Point start(ray.getOrigin());
			DragBox::Vector direction(ray.getDirection());
			
			/* Pick the active polypeptide chain's drag box: */
			{
			Threads::Mutex::Lock proteinLock(application->activeChain->proteinMutex);
			if(application->activeChain->interactor->pickDragBox(start,direction))
				{
				/* Start dragging the active chain: */
				dragChain=application->activeChain;
				}
			}
			}
		else
			{
			/* Get the device position and orientation in the active chain's coordinates: */
			Vrui::NavTransform chainTrans=Vrui::NavTransform(getButtonDeviceTransformation(0));
			chainTrans.leftMultiply(Vrui::getInverseNavigationTransformation());
			chainTrans.leftMultiply(application->activeChain->transformation);
			DragBox::Transformation transform(chainTrans.getTranslation(),chainTrans.getRotation());
			
			/* Pick the active polypeptide chain's drag box: */
			{
			Threads::Mutex::Lock proteinLock(application->activeChain->proteinMutex);
			if(application->activeChain->interactor->pickDragBox(transform))
				{
				/* Start dragging the active chain: */
				dragChain=application->activeChain;
				}
			}
			}
		}
	else // Button was just released
		{
		if(dragChain!=0)
			{
			/* Check if the IK thread is already done with the final IK request: */
			{
			Threads::Mutex::Lock proteinLock(dragChain->proteinMutex);
			if(application->ikResult==application->ikRequest)
				{
				/* Finish the IK sequence: */
				application->finishIKSequence(dragChain);
				}
			
			/* Mark the most recent IK request as the last one: */
			application->ikLastRequest=application->ikRequest;
			}
			
			/* Release the drag chain: */
			dragChain=0;
			}
		}
	}

void VRProtoShop::StructureDraggingTool::frame(void)
	{
	if(dragChain!=0)
		{
		if(rayBased)
			{
			/* Get the device ray in the active chain's coordinates: */
			Vrui::Ray ray=getButtonDeviceRay(0);
			ray.transform(Vrui::getInverseNavigationTransformation());
			ray.transform(Geometry::invert(dragChain->transformation));
			ray.normalizeDirection();
			
			/* Drag the dragging box: */
			DragBox::Point start(ray.getOrigin());
			DragBox::Vector direction(ray.getDirection());
			{
			Threads::Mutex::Lock proteinLock(dragChain->proteinMutex);
			dragChain->interactor->dragBox(start,direction);
			}
			}
		else
			{
			/* Get the device position and orientation in navigation coordinates: */
			Vrui::NavTransform chainTrans=Vrui::NavTransform(getButtonDeviceTransformation(0));
			chainTrans.leftMultiply(Vrui::getInverseNavigationTransformation());
			chainTrans.leftMultiply(dragChain->transformation);
			DragBox::Transformation transform(chainTrans.getTranslation(),chainTrans.getRotation());
			
			/* Drag the dragging box: */
			{
			Threads::Mutex::Lock proteinLock(dragChain->proteinMutex);
			dragChain->interactor->dragBox(transform);
			}
			}
		
		/* Post an update request to the IK thread: */
		{
		Threads::MutexCond::Lock ikRequestLock(application->ikRequestCond);
		++application->ikRequest;
		application->ikChain=dragChain;
		application->ikGoalTransformation=dragChain->interactor->getDragTransformation();
		application->ikRequestCond.signal();
		}
		}
	}

/****************************
Methods of class VRProtoShop:
****************************/

void VRProtoShop::finishIKSequence(Chain* ikChain)
	{
	/* Reset the chain's interactor: */
	ikChain->interactor->finishInteraction();
	
	/* Check if the active polypeptide chain has changed: */
	if(activeChain!=ikChain)
		{
		/* Unselect the manipulated chain: */
		ikChain->interactor->unselect();
		}
	
	/* Update the states of the undo/redo buttons: */
	undoButton->setEnabled(undoBuffer->canUndo());
	redoButton->setEnabled(undoBuffer->canRedo());
	}

void VRProtoShop::undoCallback(Misc::CallbackData* cbData)
	{
	/* Forward to the undo buffer's method: */
	MD::Protein* protein=undoBuffer->undo();
	
	/* Update derived protein state: */
	for(ChainList::iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		if(cIt->protein==protein)
			{
			/* Invalidate the protein's graphics cache: */
			cIt->proteinRenderer->updateProtein();
			
			/* Reset the drag box: */
			cIt->interactor->resetDragBox();
			
			/* Stop looking: */
			break;
			}
	
	/* Update the states of the undo/redo buttons: */
	undoButton->setEnabled(undoBuffer->canUndo());
	redoButton->setEnabled(undoBuffer->canRedo());
	}

void VRProtoShop::redoCallback(Misc::CallbackData* cbData)
	{
	/* Forward to the undo buffer's method: */
	MD::Protein* protein=undoBuffer->redo();
	
	/* Update derived protein state: */
	for(ChainList::iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		if(cIt->protein==protein)
			{
			/* Invalidate the protein's graphics cache: */
			cIt->proteinRenderer->updateProtein();
			
			/* Reset the drag box: */
			cIt->interactor->resetDragBox();
			
			/* Stop looking: */
			break;
			}
	
	/* Update the states of the undo/redo buttons: */
	undoButton->setEnabled(undoBuffer->canUndo());
	redoButton->setEnabled(undoBuffer->canRedo());
	}

void VRProtoShop::showChainDialogCallback(Misc::CallbackData* cbData)
	{
	/* Pop up the polypeptide chain dialog: */
	Vrui::popupPrimaryWidget(chainDialog);
	}

void VRProtoShop::showResidueDialogCallback(Misc::CallbackData* cbData)
	{
	/* Pop up the amino acid residue dialog: */
	Vrui::popupPrimaryWidget(residueDialog);
	}

void VRProtoShop::showInteractionDialogCallback(Misc::CallbackData* cbData)
	{
	/* Pop up the interaction dialog: */
	Vrui::popupPrimaryWidget(interactionDialog);
	}

GLMotif::PopupMenu* VRProtoShop::createMainMenu(void)
	{
	GLMotif::PopupMenu* mainMenu=new GLMotif::PopupMenu("MainMenu",Vrui::getWidgetManager());
	mainMenu->setTitle("VR ProtoShop");
	
	undoButton=new GLMotif::Button("UndoButton",mainMenu,"Undo");
	undoButton->setEnabled(undoBuffer->canUndo());
	undoButton->getSelectCallbacks().add(this,&VRProtoShop::undoCallback);
	
	redoButton=new GLMotif::Button("RedoButton",mainMenu,"Redo");
	redoButton->setEnabled(undoBuffer->canRedo());
	redoButton->getSelectCallbacks().add(this,&VRProtoShop::redoCallback);
	
	GLMotif::Button* showChainDialogButton=new GLMotif::Button("ShowChainDialogButton",mainMenu,"Show Chain Dialog");
	showChainDialogButton->getSelectCallbacks().add(this,&VRProtoShop::showChainDialogCallback);
	
	GLMotif::Button* showResidueDialogButton=new GLMotif::Button("ShowResidueDialogButton",mainMenu,"Show Residue Dialog");
	showResidueDialogButton->getSelectCallbacks().add(this,&VRProtoShop::showResidueDialogCallback);
	
	GLMotif::Button* showInteractionDialogButton=new GLMotif::Button("ShowInteractionDialogButton",mainMenu,"Show Interaction Dialog");
	showInteractionDialogButton->getSelectCallbacks().add(this,&VRProtoShop::showInteractionDialogCallback);
	
	mainMenu->manageMenu();
	
	return mainMenu;
	}

void VRProtoShop::updateChainDialogToggles(void)
	{
	chainVisibleToggle->setToggle(activeChain!=0&&activeChain->visible);
	chainDrawAtomsToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawAtoms());
	chainDrawBondsToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawBonds());
	chainDrawBackboneToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawBackboneRibbon());
	chainDrawCartoonToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawCartoon());
	chainDrawHydrogenBondsToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawHydrogenBonds());
	chainDrawHydrogenBondSitesToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawHydrogenBondSites());
	chainDrawHydrogenCagesToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawHydrogenCages());
	chainDrawAtomCollisionsToggle->setToggle(activeChain!=0&&activeChain->proteinRenderer->getDrawCollisions());
	}

void VRProtoShop::updateChainDialog(void)
	{
	/* Replace the chain list box's contents: */
	chainListBox->clear();
	int activeChainIndex=-1;
	for(size_t i=0;i<chains.size();++i)
		{
		chainListBox->addItem(chains[i].name.c_str());
		if(activeChain==&chains[i])
			activeChainIndex=int(i);
		}
	if(activeChainIndex>=0)
		chainListBox->selectItem(activeChainIndex);
	
	/* Enable or disable all chain widgets: */
	chainVisibleToggle->setEnabled(activeChainIndex>=0);
	chainDrawAtomsToggle->setEnabled(activeChainIndex>=0);
	chainDrawBondsToggle->setEnabled(activeChainIndex>=0);
	chainDrawBackboneToggle->setEnabled(activeChainIndex>=0);
	chainDrawCartoonToggle->setEnabled(activeChainIndex>=0);
	chainDrawHydrogenBondsToggle->setEnabled(activeChainIndex>=0);
	chainDrawHydrogenBondSitesToggle->setEnabled(activeChainIndex>=0);
	chainDrawHydrogenCagesToggle->setEnabled(activeChainIndex>=0);
	chainDrawAtomCollisionsToggle->setEnabled(activeChainIndex>=0);
	resetTransformButton->setEnabled(activeChainIndex>=0);
	
	/* Update the chain rendering toggles: */
	updateChainDialogToggles();
	}

void VRProtoShop::chainListBoxValueChangedCallback(GLMotif::ListBox::ValueChangedCallbackData* cbData)
	{
	/* Reset the currently active chain if there is no on-going IK operation: */
	{
	Threads::Mutex::Lock proteinLock(activeChain->proteinMutex);
	if(ikResult==ikLastRequest)
		activeChain->interactor->unselect();
	}
	
	/* Update the active chain: */
	activeChain=&chains[cbData->newSelectedItem];
	
	/* Update the chain rendering toggles: */
	updateChainDialogToggles();
	}

void VRProtoShop::chainVisibleToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Hide or show the active chain: */
	chains[chainListBox->getSelectedItem()].visible=cbData->set;
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

namespace {

/****************
Helper functions:
****************/

void checkRenderFlags(MD::ProteinRenderer* renderer)
	{
	/* Enable basic backbone drawing if no other rendering flags are enabled: */
	bool enableBackbone=!(renderer->getDrawAtoms()||renderer->getDrawBonds()||renderer->getDrawBackboneRibbon()||renderer->getDrawCartoon());
	renderer->setDrawBackbone(enableBackbone);
	}

}

void VRProtoShop::chainDrawAtomsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawAtoms(cbData->set);
	checkRenderFlags(chains[chainListBox->getSelectedItem()].proteinRenderer);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::chainDrawBondsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawBonds(cbData->set);
	checkRenderFlags(chains[chainListBox->getSelectedItem()].proteinRenderer);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::chainDrawBackboneToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawBackboneRibbon(cbData->set);
	checkRenderFlags(chains[chainListBox->getSelectedItem()].proteinRenderer);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::chainDrawCartoonToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawCartoon(cbData->set);
	checkRenderFlags(chains[chainListBox->getSelectedItem()].proteinRenderer);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::chainDrawHydrogenBondsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawHydrogenBonds(cbData->set);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::chainDrawHydrogenBondSitesToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawHydrogenBondSites(cbData->set);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::chainDrawHydrogenCagesToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawHydrogenCages(cbData->set);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::chainDrawAtomCollisionsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Update the selected chain's renderer: */
	chains[chainListBox->getSelectedItem()].proteinRenderer->setDrawCollisions(cbData->set);
	if(psClient!=0)
		psClient->updateRenderFlags();
	}

void VRProtoShop::resetTransformCallback(Misc::CallbackData* cbData)
	{
	if(psClient!=0)
		{
		/* Ask the server to reset the active chain's transformation; will be done when server's reply arrives: */
		psClient->resetChainTransform(chainListBox->getSelectedItem());
		}
	else
		{
		/* Reset the active chain's transformation: */
		chains[chainListBox->getSelectedItem()].transformation=IK::Transformation::identity;
		}
	}

GLMotif::PopupWindow* VRProtoShop::createChainDialog(void)
	{
	GLMotif::PopupWindow* chainPopup=new GLMotif::PopupWindow("ChainPopup",Vrui::getWidgetManager(),"Polypeptide Chains");
	chainPopup->setCloseButton(true);
	chainPopup->setResizableFlags(true,true);
	chainPopup->popDownOnClose();
	
	GLMotif::RowColumn* chainPanel=new GLMotif::RowColumn("ChainPanel",chainPopup,false);
	chainPanel->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	chainPanel->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	chainPanel->setNumMinorWidgets(1);
	
	chainListBox=new GLMotif::ListBox("ChainListBox",chainPanel,GLMotif::ListBox::ALWAYS_ONE,20,10,false);
	chainListBox->getValueChangedCallbacks().add(this,&VRProtoShop::chainListBoxValueChangedCallback);
	
	chainListBox->manageChild();
	
	GLMotif::Margin* renderTogglesMargin=new GLMotif::Margin("RenderTogglesMargin",chainPanel,false);
	renderTogglesMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::HFILL,GLMotif::Alignment::VCENTER));
	
	GLMotif::RowColumn* renderTogglesBox=new GLMotif::RowColumn("RenderTogglesBox",renderTogglesMargin,false);
	renderTogglesBox->setOrientation(GLMotif::RowColumn::VERTICAL);
	renderTogglesBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	renderTogglesBox->setNumMinorWidgets(1);
	
	chainVisibleToggle=new GLMotif::ToggleButton("ChainVisibleToggle",renderTogglesBox,"Draw Chain");
	chainVisibleToggle->setBorderWidth(0.0f);
	chainVisibleToggle->setHAlignment(GLFont::Left);
	chainVisibleToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainVisibleToggleCallback);
	
	new GLMotif::Separator("Sep1",renderTogglesBox,GLMotif::Separator::HORIZONTAL,0.0f,GLMotif::Separator::LOWERED);
	
	chainDrawAtomsToggle=new GLMotif::ToggleButton("DrawAtomsToggle",renderTogglesBox,"Draw atoms");
	chainDrawAtomsToggle->setBorderWidth(0.0f);
	chainDrawAtomsToggle->setHAlignment(GLFont::Left);
	chainDrawAtomsToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawAtomsToggleCallback);
	
	chainDrawBondsToggle=new GLMotif::ToggleButton("DrawBondsToggle",renderTogglesBox,"Draw sidechain bonds");
	chainDrawBondsToggle->setBorderWidth(0.0f);
	chainDrawBondsToggle->setHAlignment(GLFont::Left);
	chainDrawBondsToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawBondsToggleCallback);
	
	chainDrawBackboneToggle=new GLMotif::ToggleButton("DrawBackboneToggle",renderTogglesBox,"Draw backbone ribbon");
	chainDrawBackboneToggle->setBorderWidth(0.0f);
	chainDrawBackboneToggle->setHAlignment(GLFont::Left);
	chainDrawBackboneToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawBackboneToggleCallback);
	
	chainDrawCartoonToggle=new GLMotif::ToggleButton("DrawCartoonToggle",renderTogglesBox,"Draw structure cartoons");
	chainDrawCartoonToggle->setBorderWidth(0.0f);
	chainDrawCartoonToggle->setHAlignment(GLFont::Left);
	chainDrawCartoonToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawCartoonToggleCallback);
	
	chainDrawHydrogenBondsToggle=new GLMotif::ToggleButton("DrawHydrogenBondsToggle",renderTogglesBox,"Draw hydrogen bonds");
	chainDrawHydrogenBondsToggle->setBorderWidth(0.0f);
	chainDrawHydrogenBondsToggle->setHAlignment(GLFont::Left);
	chainDrawHydrogenBondsToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawHydrogenBondsToggleCallback);
	
	chainDrawHydrogenBondSitesToggle=new GLMotif::ToggleButton("DrawHydrogenBondSitesToggle",renderTogglesBox,"Draw hydrogen bond sites");
	chainDrawHydrogenBondSitesToggle->setBorderWidth(0.0f);
	chainDrawHydrogenBondSitesToggle->setHAlignment(GLFont::Left);
	chainDrawHydrogenBondSitesToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawHydrogenBondSitesToggleCallback);
	
	chainDrawHydrogenCagesToggle=new GLMotif::ToggleButton("DrawHydrogenCagesToggle",renderTogglesBox,"Draw hydrogen cages");
	chainDrawHydrogenCagesToggle->setBorderWidth(0.0f);
	chainDrawHydrogenCagesToggle->setHAlignment(GLFont::Left);
	chainDrawHydrogenCagesToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawHydrogenCagesToggleCallback);
	
	chainDrawAtomCollisionsToggle=new GLMotif::ToggleButton("DrawAtomCollisionsToggle",renderTogglesBox,"Visualize atom collisions");
	chainDrawAtomCollisionsToggle->setBorderWidth(0.0f);
	chainDrawAtomCollisionsToggle->setHAlignment(GLFont::Left);
	chainDrawAtomCollisionsToggle->getValueChangedCallbacks().add(this,&VRProtoShop::chainDrawAtomCollisionsToggleCallback);
	
	new GLMotif::Separator("Sep2",renderTogglesBox,GLMotif::Separator::HORIZONTAL,0.0f,GLMotif::Separator::LOWERED);
	
	resetTransformButton=new GLMotif::Button("ResetTransformButton",renderTogglesBox,"Reset Transform");
	resetTransformButton->getSelectCallbacks().add(this,&VRProtoShop::resetTransformCallback);
	
	renderTogglesBox->manageChild();
	
	renderTogglesMargin->manageChild();
	
	chainPanel->setColumnWeight(0,1.0f);
	chainPanel->manageChild();
	
	/* Update the chain dialog to match initial state: */
	updateChainDialog();
	
	return chainPopup;
	}

void VRProtoShop::updateResidueFields(MD::Protein::Residue* residue)
	{
	if(residue!=0)
		{
		/* Update the residue name: */
		residueTypeTextField->setEnabled(true);
		residueTypeTextField->setString(residue->getPdbResidueName());
		
		/* Update the residue structure type: */
		structureTypeBox->setEnabled(true);
		switch(residue->getSecondaryStructure()->getStructureType())
			{
			case MD::Protein::SecondaryStructure::NONE:
				structureTypeBox->setSelectedToggle(-1);
				break;
			
			case MD::Protein::SecondaryStructure::COIL:
				structureTypeBox->setSelectedToggle(0);
				break;
			
			case MD::Protein::SecondaryStructure::ALPHA_HELIX:
				structureTypeBox->setSelectedToggle(1);
				break;
			
			case MD::Protein::SecondaryStructure::BETA_STRAND:
				structureTypeBox->setSelectedToggle(2);
				break;
			}
		}
	else
		{
		/* Disable residue UI elements: */
		residueTypeTextField->setEnabled(false);
		residueTypeTextField->setString("");
		structureTypeBox->setEnabled(false);
		structureTypeBox->setSelectedToggle(-1);
		}
	}

void VRProtoShop::updateResidueDialog(void)
	{
	/* Check if there is an active chain and a selected residue: */
	if(activeChain!=0)
		{
		/* Enable the residue selection slider: */
		residueIndexSlider->setEnabled(true);
		
		/* Update the range of the chain's residue indices: */
		std::pair<int,int> residueIndexRange=activeChain->protein->getResidueIndexRange();
		residueIndexSlider->setValueRange(residueIndexRange.first-1,residueIndexRange.second,1.0);
		
		/* Update the residue index slider: */
		MD::Protein::Residue* residue=activeChain->interactor->getSelectedResidue();
		if(residue!=0)
			residueIndexSlider->setValue(residue->getPdbResidueIndex());
		else
			{
			/* Set the slider to the "none selected" index: */
			residueIndexSlider->setValue(residueIndexRange.first-1);
			}
		
		/* Update the residue fields: */
		updateResidueFields(residue);
		}
	else
		{
		/* Disable all UI elements: */
		residueIndexSlider->setEnabled(false);
		residueTypeTextField->setEnabled(false);
		residueTypeTextField->setString("");
		structureTypeBox->setEnabled(false);
		structureTypeBox->setSelectedToggle(-1);
		}
	}

void VRProtoShop::residueIndexSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Select the selected residue: */
	int index(Math::floor(cbData->value+0.5));
	MD::Protein::Residue* residue=index>=activeChain->protein->getResidueIndexRange().first?activeChain->protein->pickResidue(index):0;
	activeChain->interactor->setSelectedResidue(residue);
	
	/* Update the residue fields: */
	updateResidueFields(residue);
	}

void VRProtoShop::structureTypeBoxCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData)
	{
	if(cbData->newSelectedToggle!=0)
		{
		switch(cbData->radioBox->getToggleIndex(cbData->newSelectedToggle))
			{
			case 0:
				activeChain->protein->changeResidueStructureType(activeChain->interactor->getSelectedResidue(),MD::Protein::SecondaryStructure::COIL);
				break;
			
			case 1:
				activeChain->protein->changeResidueStructureType(activeChain->interactor->getSelectedResidue(),MD::Protein::SecondaryStructure::ALPHA_HELIX);
				break;
			
			case 2:
				activeChain->protein->changeResidueStructureType(activeChain->interactor->getSelectedResidue(),MD::Protein::SecondaryStructure::BETA_STRAND);
				break;
			}
		}
	else
		activeChain->protein->changeResidueStructureType(activeChain->interactor->getSelectedResidue(),MD::Protein::SecondaryStructure::NONE);
	
	activeChain->proteinRenderer->updateProtein();
	}

GLMotif::PopupWindow* VRProtoShop::createResidueDialog(void)
	{
	const GLMotif::StyleSheet& ss=*Vrui::getUiStyleSheet();
	
	GLMotif::PopupWindow* residuePopup=new GLMotif::PopupWindow("ResiduePopup",Vrui::getWidgetManager(),"Amino Acid Residues");
	residuePopup->setCloseButton(true);
	residuePopup->setResizableFlags(true,true);
	residuePopup->popDownOnClose();
	
	GLMotif::RowColumn* residuePanel=new GLMotif::RowColumn("ResiduePanel",residuePopup,false);
	residuePanel->setOrientation(GLMotif::RowColumn::VERTICAL);
	residuePanel->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	residuePanel->setNumMinorWidgets(2);
	
	new GLMotif::Label("ResidueIndexLabel",residuePanel,"Residue Index");
	
	residueIndexSlider=new GLMotif::TextFieldSlider("ResidueIndexSlider",residuePanel,6,ss.fontHeight*10.0f);
	residueIndexSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	residueIndexSlider->setValueType(GLMotif::TextFieldSlider::INT);
	residueIndexSlider->getValueChangedCallbacks().add(this,&VRProtoShop::residueIndexSliderCallback);
	
	new GLMotif::Label("ResidueTypeLabel",residuePanel,"Residue Type");
	
	residueTypeTextField=new GLMotif::TextField("ResidueTypeTextField",residuePanel,6);
	residueTypeTextField->setValueType(GLMotif::TextField::ALPHA);
	residueTypeTextField->setHAlignment(GLFont::HAlignment::Left);
	
	new GLMotif::Label("StructureTypeLabel",residuePanel,"Structure Type");
	
	structureTypeBox=new GLMotif::RadioBox("StructureTypeBox",residuePanel,false);
	structureTypeBox->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	structureTypeBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	structureTypeBox->setNumMinorWidgets(1);
	structureTypeBox->setSelectionMode(GLMotif::RadioBox::ATMOST_ONE);
	
	structureTypeBox->addToggle("Coil");
	structureTypeBox->addToggle("Alpha Helix");
	structureTypeBox->addToggle("Beta Strand");
	structureTypeBox->setSelectedToggle(-1);
	
	structureTypeBox->getValueChangedCallbacks().add(this,&VRProtoShop::structureTypeBoxCallback);
	
	structureTypeBox->manageChild();
	
	residuePanel->manageChild();
	
	updateResidueDialog();
	
	return residuePopup;
	}

void VRProtoShop::ikUpdateDirectionValueChangedCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData)
	{
	}

void VRProtoShop::resetActiveCoilRegionsCallback(Misc::CallbackData*)
	{
	}

void VRProtoShop::autoBondingValueChangedCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData)
	{
	}

GLMotif::PopupWindow* VRProtoShop::createInteractionDialog(void)
	{
	GLMotif::PopupWindow* interactionPopup=new GLMotif::PopupWindow("InteractionPopup",Vrui::getWidgetManager(),"Interaction Dialog");
	interactionPopup->setCloseButton(true);
	interactionPopup->setResizableFlags(false,false);
	interactionPopup->popDownOnClose();
	
	GLMotif::RowColumn* interactionPanel=new GLMotif::RowColumn("InteractionPanel",interactionPopup,false);
	
	GLMotif::RadioBox* ikUpdateDirectionBox=new GLMotif::RadioBox("IkUpdateDirectionBox",interactionPanel,false);
	ikUpdateDirectionBox->setSelectionMode(GLMotif::RadioBox::ALWAYS_ONE);
	
	ikUpdateDirectionBox->addToggle("Left (C to N)");
	ikUpdateDirectionBox->addToggle("Right (N to C)");
	
	ikUpdateDirectionBox->manageChild();
	switch(updateDirection)
		{
		case UPDATE_LEFT:
			ikUpdateDirectionBox->setSelectedToggle(0);
			break;
		
		case UPDATE_RIGHT:
			ikUpdateDirectionBox->setSelectedToggle(0);
			break;
		}
	ikUpdateDirectionBox->getValueChangedCallbacks().add(this,&VRProtoShop::ikUpdateDirectionValueChangedCallback);
	
	GLMotif::Button* resetActiveCoilRegionsButton=new GLMotif::Button("ResetActiveCoilRegionsButton",interactionPanel,"Reset Active Coil Regions");
	resetActiveCoilRegionsButton->getSelectCallbacks().add(this,&VRProtoShop::resetActiveCoilRegionsCallback);
	
	autoBondingBox=new GLMotif::RadioBox("AutoBondingBox",interactionPanel,false);
	autoBondingBox->setSelectionMode(GLMotif::RadioBox::ATMOST_ONE);
	
	autoBondingBox->addToggle("Parallel");
	autoBondingBox->addToggle("Anti-parallel");
	
	autoBondingBox->manageChild();
	autoBondingBox->setSelectedToggle(0);
	autoBondingBox->getValueChangedCallbacks().add(this,&VRProtoShop::autoBondingValueChangedCallback);
	
	interactionPanel->manageChild();
	
	return interactionPopup;
	}

void* VRProtoShop::ikThreadMethod(void)
	{
	/* Calculate IK requests until shut down: */
	while(true)
		{
		/* Wait while there is no pending IK request: */
		int request;
		Chain* chain;
		ProteinInteractor::Transformation goalTransformation;
		{
		Threads::MutexCond::Lock ikRequestLock(ikRequestCond);
		while(ikKeepRunning&&ikRequest==ikResult)
			ikRequestCond.wait(ikRequestLock);
		request=ikRequest;
		chain=ikChain;
		goalTransformation=ikGoalTransformation;
		}
		
		/* Bail out if shut down: */
		if(!ikKeepRunning)
			break;
		
		/* Perform IK steps: */
		
		// DEBUGGING
		Realtime::TimePointMonotonic now;
		
		chain->interactor->drag(goalTransformation);
		
		// DEBUGGING
		std::cout<<"IK time "<<double(now.setAndDiff())*1000.0<<" ms"<<std::endl;
		
		/* Apply calculated dihedral angle changes to the active polypeptide chain: */
		{
		Threads::Mutex::Lock proteinLock(chain->proteinMutex);
		
		// DEBUGGING
		now.set();
		
		chain->interactor->applyChanges();
		
		// DEBUGGING
		std::cout<<"Apply time "<<double(now.setAndDiff())*1000.0<<" ms"<<std::endl;
		
		/* Check if the current dragging operation is done: */
		if(ikLastRequest==request)
			{
			/* Finish the IK interaction: */
			finishIKSequence(chain);
			}
		
		/* Mark the request as fulfilled: */
		ikResult=request;
		}
		}
	
	return 0;
	}

VRProtoShop::VRProtoShop(int& argc,char**& argv)
	:Vrui::Application(argc,argv),
	 clusterPipe(Vrui::openPipe()),
	 configFile(VRPROTOSHOP_CONFIG_ETCDIR "/" VRPROTOSHOP_CONFIG_CONFIGFILENAME),
	 psClient(0),
	 drawDensity(false),densityRenderer(0),densityPalette(0),
	 undoBuffer(new UndoBuffer),
	 activeChain(0),
	 autoBondingMode(AUTOBOND_NONE),
	 updateDirection(UPDATE_RIGHT),
	 interactorPipe(Vrui::openPipe()),
	 ikRequest(0),ikLastRequest(0),ikResult(0),ikKeepRunning(false),
	 mainMenu(0),chainDialog(0),residueDialog(0),interactionDialog(0)
	{
	/* Check if there is a collaboration client: */
	Collab::Client* client=Collab::Client::getTheClient();
	if(client!=0)
		{
		/* Register a ProtoShop client: */
		psClient=new ProtoShopClient(client,this);
		client->addPluginProtocol(psClient);
		}
	else
		{
		/* Scan the command line: */
		std::vector<char*> inputFileNames;
		const char* densityFileName=0;
		const char* paletteFileName=0;
		for(int i=1;i<argc;++i)
			{
			if(argv[i][0]=='-')
				{
				if(strcasecmp(argv[i]+1,"vol")==0||strcasecmp(argv[i]+1,"fvol")==0)
					{
					/* Treat the next argument as a density file name: */
					++i;
					densityFileName=argv[i];
					}
				else if(strcasecmp(argv[i]+1,"pal")==0)
					{
					/* Treat the next argument as a palette file name: */
					++i;
					paletteFileName=argv[i];
					}
				}
			else
				{
				/* Treat this argument as a protein file name: */
				inputFileNames.push_back(argv[i]);
				}
			}
		
		/* Initialize protein creator: */
		MD::initializeProteinCreation(configFile.getSection("/ProteinCreation"));
		
		/* Process all input file names: */
		if(inputFileNames.empty())
			throw std::runtime_error("No prediction or PDB file name(s) provided");
		
		for(std::vector<char*>::iterator ifnIt=inputFileNames.begin();ifnIt!=inputFileNames.end();++ifnIt)
			{
			/* Query base name and extension of input file: */
			const char* baseNamePtr=*ifnIt;
			const char* extPtr=0;
			const char* cPtr;
			for(cPtr=*ifnIt;*cPtr!='\0';++cPtr)
				{
				if(*cPtr=='/')
					baseNamePtr=cPtr+1;
				else if(*cPtr=='.')
					extPtr=cPtr;
				}
			if(extPtr==0)
				extPtr=cPtr;
			
			/* Load or create a protein model: */
			std::vector<MD::Protein*> fileChains;
			if(strcasecmp(extPtr,".pdb")==0)
				{
				/* Load a protein model: */
				fileChains=MD::parsePdbFile(*ifnIt);
				}
			else if(strcasecmp(extPtr,".pred")==0)
				{
				/* Build a protein model: */
				MD::ReadStandards(VRPROTOSHOP_CONFIG_STANDARDSDIR);
				fileChains.push_back(MD::ReadPredictionFile(*ifnIt));
				}
			else
				throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Unrecognized input file name extension %s",extPtr);
			
			unsigned int chainIndex=0;
			for(std::vector<MD::Protein*>::iterator cIt=fileChains.begin();cIt!=fileChains.end();++cIt,++chainIndex)
				{
				/* Create a new chain structure: */
				Chain newChain;
				newChain.name=std::string(baseNamePtr,extPtr);
				if(fileChains.size()>1)
					{
					newChain.name.push_back(' ');
					newChain.name.push_back('A'+(cIt-fileChains.begin()));
					}
				newChain.protein=*cIt;
				newChain.proteinRenderer=new MD::ProteinRenderer(configFile.getSection("/ProteinRenderer"),newChain.protein);
				newChain.interactor=new ProteinInteractor(chainIndex,newChain.protein,newChain.proteinRenderer,0,undoBuffer,interactorPipe);
				
				#if 0
				/* Disable bond cylinder rendering for all structures, but enable it for the entire protein: */
				for(int i=0;i<newChain.protein->getNumStructures();++i)
					newChain.proteinRenderer->setDrawBonds(newChain.protein->pickStructure(i),false);
				newChain.proteinRenderer->setDrawBonds(true);
				#endif
				
				chains.push_back(newChain);
				}
			}
		
		if(densityFileName!=0)
			{
			/* Load palette renderer from file: */
			densityRenderer=new PaletteRenderer(densityFileName);
			
			/* Initialize renderer settings: */
			// densityRenderer->setUseNPOTDTextures(true);
			densityRenderer->setVoxelAlignment(VolumeRenderer::VERTEX_CENTERED);
			// densityRenderer->setRenderingMode(VolumeRenderer::AXIS_ALIGNED);
			densityRenderer->setRenderingMode(VolumeRenderer::VIEW_PERPENDICULAR);
			densityRenderer->setInterpolationMode(VolumeRenderer::LINEAR);
			densityRenderer->setTextureFunction(VolumeRenderer::REPLACE);
			densityRenderer->setSliceFactor(VolumeRenderer::Scalar(1.4142));
			densityRenderer->setAutosaveGLState(true);
			densityRenderer->setTextureCaching(true);
			
			/* Create a color map: */
			if(paletteFileName!=0)
				densityPalette=new GLColorMap(paletteFileName,0.0f,1.0f);
			else
				densityPalette=new GLColorMap(GLColorMap::RAINBOW|GLColorMap::RAMP_ALPHA,1.0f,1.0f,0.0f,1.0f);
			densityPalette->changeTransparency(1.0);
			densityPalette->premultiplyAlpha();
			densityRenderer->setSharePalette(false);
			densityRenderer->setColorMap(densityPalette);
			}
		
		/* Activate the first chain: */
		activeChain=&chains.front();
		}
	
	/* Create the main menu: */
	mainMenu=createMainMenu();
	Vrui::setMainMenu(mainMenu);
	
	/* Create the polypeptide chain dialog: */
	chainDialog=createChainDialog();
	
	/* Create the amino acid residue dialog: */
	residueDialog=createResidueDialog();
	
	/* Create the interaction dialog: */
	interactionDialog=createInteractionDialog();
	
	/* Create all tool classes: */
	Vrui::ToolFactory* baseClass=Tool::initClass();
	addEventTool("Select Secondary Structure",baseClass,0);
	addEventTool("Select Residue",baseClass,1);
	addEventTool("Toggle Coil Region",baseClass,2);
	ChainDraggingTool::initClass(baseClass);
	StructureDraggingTool::initClass(baseClass);
	
	/* Start the IK calculation thread: */
	ikKeepRunning=true;
	ikThread.start(this,&VRProtoShop::ikThreadMethod);
	}

VRProtoShop::~VRProtoShop(void)
	{
	delete mainMenu;
	delete chainDialog;
	delete residueDialog;
	delete interactionDialog;
	
	if(!ikThread.isJoined())
		{
		/* Shut down the IK calculation thread: */
		ikKeepRunning=false;
		ikRequestCond.signal();
		ikThread.join();
		}
	
	/* Delete polypeptide chain state: */
	for(ChainList::iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		{
		delete cIt->interactor;
		delete cIt->proteinRenderer;
		delete cIt->protein;
		}
	delete interactorPipe;
	delete undoBuffer;
	delete densityPalette;
	delete densityRenderer;
	
	/* Shut down cluster communication: */
	delete clusterPipe;
	}

void VRProtoShop::frame(void)
	{
	/* Set the viewing direction for illuminated lines and density rendering: */
	viewDirection=VolumeRenderer::Vector(Vrui::getViewDirection());
	GLLineIlluminator::Vector lineViewDirection(viewDirection.getComponents());
	
	/* Set the model scales of all protein renderers: */
	for(ChainList::iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		{
		cIt->proteinRenderer->setModelScale(Vrui::getNavigationTransformation().getScaling());
		cIt->proteinRenderer->setViewDirection(lineViewDirection);
		cIt->proteinRenderer->setLightDirection(lineViewDirection);
		}
	}

void VRProtoShop::display(GLContextData& contextData) const
	{
	/* Draw all visible protein chains: */
	for(ChainList::const_iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		if(cIt->visible)
			{
			glPushMatrix();
			glMultMatrix(cIt->transformation);
			
			/* Draw the chain: */
			{
			Threads::Mutex::Lock proteinLock(cIt->proteinMutex);
			cIt->proteinRenderer->glRenderAction(contextData);
			}
			
			glPopMatrix();
			}
	
	if(densityRenderer!=0&&drawDensity)
		{
		/* Render the density distribution: */
		densityRenderer->renderBlock(contextData,viewDirection);
		}
	
	/* Draw all structure selection boxes: */
	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	for(ChainList::const_iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		if(cIt->visible)
			{
			glPushMatrix();
			glMultMatrix(cIt->transformation);
			
			cIt->interactor->glRenderAction(contextData);
			
			glPopMatrix();
			}
	glPopAttrib();
	}

void VRProtoShop::resetNavigation(void)
	{
	/* Calculate the overall model's position and size: */
	Geometry::Box<Vrui::Scalar,3> bbox=Geometry::Box<Vrui::Scalar,3>::empty;
	for(ChainList::iterator cIt=chains.begin();cIt!=chains.end();++cIt)
		{
		Vrui::Point modelCenter=Vrui::Point(cIt->protein->calcCentroid());
		Vrui::Scalar modelSize=Vrui::Scalar(cIt->protein->calcRadius());
		bbox.addBox(Geometry::Box<Vrui::Scalar,3>(modelCenter-Vrui::Vector(modelSize),modelCenter+Vrui::Vector(modelSize)));
		}
	
	/* Set the navigation transformation: */
	Vrui::setNavigationTransformation(Geometry::mid(bbox.min,bbox.max),Geometry::dist(bbox.min,bbox.max));
	}

void VRProtoShop::eventCallback(Vrui::Application::EventID eventId,Vrui::InputDevice::ButtonCallbackData* cbData)
	{
	/* Bail out if this isn't a button press event, or there isn't an active polypeptide chain: */
	if(!cbData->newButtonState||activeChain==0)
		return;
	
	/* Get the initiating device's position and ray in the active chain's coordinates: */
	Vrui::Point pos;
	Vrui::Ray ray;
	bool rayBased=cbData->inputDevice->isRayDevice();
	if(rayBased)
		{
		/* Convert the device's ray to navigation coordinates: */
		ray=cbData->inputDevice->getRay();
		ray.transform(Vrui::getInverseNavigationTransformation());
		ray.transform(Geometry::invert(activeChain->transformation));
		ray.normalizeDirection();
		}
	else
		{
		/* Convert the device's position to navigation coordinates: */
		pos=Vrui::getInverseNavigationTransformation().transform(cbData->inputDevice->getPosition());
		pos=activeChain->transformation.inverseTransform(pos);
		}
	
	{
	Threads::Mutex::Lock proteinLock(activeChain->proteinMutex);
	switch(eventId)
		{
		case 0: // Select secondary structure
			if(rayBased)
				activeChain->interactor->selectStructure(activeChain->protein->pickStructure(ray));
			else
				activeChain->interactor->selectStructure(activeChain->protein->pickStructure(pos));
			break;
		
		case 1: // Select residue
			if(rayBased)
				activeChain->interactor->setSelectedResidue(activeChain->protein->pickResidue(ray));
			else
				activeChain->interactor->setSelectedResidue(activeChain->protein->pickResidue(pos));
			updateResidueDialog();
			break;
		
		case 2: // Toggle coil region
			if(rayBased)
				activeChain->interactor->toggleCoil(activeChain->protein->pickStructure(ray));
			else
				activeChain->interactor->toggleCoil(activeChain->protein->pickStructure(pos));
			break;
		}
	}
	}

/****************
Global functions:
****************/

void updateProteinNow(void)
	{
	Vrui::requestUpdate();
	}

VRUI_APPLICATION_RUN(VRProtoShop)
