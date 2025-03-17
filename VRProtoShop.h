/***********************************************************************
VRProtoShop - Virtual reality version of ProtoShop (will hopefully be
fully functional at some point).
Copyright (c) 2002-2024 Oliver Kreylos

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

#ifndef VRPROTOSHOP_INCLUDED
#define VRPROTOSHOP_INCLUDED

#include <string>
#include <vector>
#include <Misc/ConfigurationFile.h>
#include <Threads/Thread.h>
#include <Threads/Mutex.h>
#include <Threads/MutexCond.h>
#include <Geometry/Ray.h>
#include <GL/gl.h>
#include <GLMotif/RadioBox.h>
#include <GLMotif/ToggleButton.h>
#include <GLMotif/ListBox.h>
#include <GLMotif/TextFieldSlider.h>
#include <Vrui/Tool.h>
#include <Vrui/GenericToolFactory.h>
#include <Vrui/GenericAbstractToolFactory.h>
#include <Vrui/ToolManager.h>
#include <Vrui/Application.h>
#include <PaletteRenderer.h>

#include "IK.h"
#include "ProteinInteractor.h"

/* Forward declarations: */
namespace Cluster {
class MulticastPipe;
}
class GLColorMap;
namespace GLMotif {
class PopupMenu;
class PopupWindow;
class Button;
class TextField;
}
namespace MD {
class Protein;
class ProteinRenderer;
}
namespace Collab {
namespace Plugins {
class ProtoShopClient;
}
}
class UndoBuffer;

/* Namespace shortcuts: */
using Collab::Plugins::ProtoShopClient;

struct VRProtoShop:public Vrui::Application
	{
	friend class ProtoShopClient;
	
	/* Embedded classes: */
	private:
	enum AutoBondingMode
		{
		AUTOBOND_NONE,AUTOBOND_PARALLEL,AUTOBOND_ANTIPARALLEL
		};
	
	enum UpdateDirection
		{
		UPDATE_LEFT,UPDATE_RIGHT
		};
	
	struct Chain // Structure describing the interaction and rendering state of a single polypeptide chain
		{
		/* Elements: */
		public:
		std::string name; // Chain's name (like protein name + chain identifier)
		IK::Transformation transformation; // Transformation applied to the entire chain
		mutable Threads::Mutex proteinMutex; // Mutex protecting the chain during interaction
		MD::Protein* protein; // The chain's 3D representation
		MD::ProteinRenderer* proteinRenderer; // Renderer attached to the chain
		ProteinInteractor* interactor; // Interaction object working on this chain
		bool visible; // Flag whether the chain is currently rendered
		
		/* Constructors and destructors: */
		Chain(void) // Creates "empty" chain
			:transformation(IK::Transformation::identity),
			 protein(0),proteinRenderer(0),interactor(0),
			 visible(true)
			{
			}
		Chain(const Chain& source) // Copy constructor
			:name(source.name),transformation(source.transformation),
			 protein(source.protein),proteinRenderer(source.proteinRenderer),interactor(source.interactor),
			 visible(source.visible)
			{
			}
		};
	
	typedef std::vector<Chain> ChainList; // Type for lists of polypeptide chains
	
	class Tool; // Forward declaration
	typedef Vrui::GenericAbstractToolFactory<Tool> ToolFactory; // Factory class for the abstract tool base class
	
	class Tool:public Vrui::Tool,public Vrui::Application::Tool<VRProtoShop> // Base class for application-specific tools
		{
		friend class Vrui::GenericAbstractToolFactory<Tool>;
		
		/* Constructors and destructors: */
		public:
		static Vrui::ToolFactory* initClass(void);
		Tool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment)
			:Vrui::Tool(factory,inputAssignment)
			{
			}
		};
	
	class ChainDraggingTool; // Forward declaration
	typedef Vrui::GenericToolFactory<ChainDraggingTool> ChainDraggingToolFactory; // Factory class for the chain dragging tool class
	
	class ChainDraggingTool:public Tool // Tool class to drag polypeptide chains
		{
		friend class Vrui::GenericToolFactory<ChainDraggingTool>;
		
		/* Elements: */
		private:
		static ChainDraggingToolFactory* factory; // Pointer to the tool class' factory
		bool rayBased; // Flag whether the dragging tool works using a ray-based or 6-DOF interface
		bool dragging; // Flag if the tool is currently dragging a chain; can be false when dragChain!=0 while waiting for lock reply
		Chain* dragChain; // Polypeptide chain with which the tool is interacting
		IK::Transformation initial; // The chain's original transformation before dragging started
		IK::Transformation drag; // The chain dragging transformation
		unsigned int lockRequestId; // ID of the dragging lock request
		unsigned int lockId; // The ID of the granted dragging lock
		
		/* Constructors and destructors: */
		public:
		static void initClass(Vrui::ToolFactory* baseClass);
		ChainDraggingTool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment);
		virtual ~ChainDraggingTool(void);
		
		/* Methods from Vrui::Tool: */
		virtual const Vrui::ToolFactory* getFactory(void) const;
		virtual void buttonCallback(int buttonSlotIndex,Vrui::InputDevice::ButtonCallbackData* cbData);
		virtual void frame(void);
		
		/* New methods: */
		void lockReply(unsigned int newLockId); // Method to notify the tool when a reply to its dragging lock request arrives; lock was not granted if newLockId==0
		};
	
	class StructureDraggingTool; // Forward declaration
	typedef Vrui::GenericToolFactory<StructureDraggingTool> StructureDraggingToolFactory; // Factory class for the structure dragging tool class
	
	class StructureDraggingTool:public Tool // Tool class to drag secondary structures
		{
		friend class Vrui::GenericToolFactory<StructureDraggingTool>;
		
		/* Elements: */
		private:
		static StructureDraggingToolFactory* factory; // Pointer to the tool class' factory
		bool rayBased; // Flag whether the dragging tool works using a ray-based or 6-DOF interface
		Chain* dragChain; // Polypeptide chain with which the tool is interacting
		
		/* Constructors and destructors: */
		public:
		static void initClass(Vrui::ToolFactory* baseClass);
		StructureDraggingTool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment);
		
		/* Methods from Vrui::Tool: */
		virtual const Vrui::ToolFactory* getFactory(void) const;
		virtual void buttonCallback(int buttonSlotIndex,Vrui::InputDevice::ButtonCallbackData* cbData);
		virtual void frame(void);
		};
	
	/* Elements: */
	private:
	Cluster::MulticastPipe* clusterPipe; // A multicast pipe to synchronize a distributed rendering cluster
	Misc::ConfigurationFile configFile; // ProtoShop configuration file
	
	/* Protein sharing state: */
	ProtoShopClient* psClient; // Client object for a shared ProtoShop protocol
	
	/* Local protein representation: */
	ChainList chains; // List of polypeptide chains being visualized and manipulated
	
	/* Global rendering stuff: */
	bool drawDensity; // Flag whether to render the density distribution
	PaletteRenderer* densityRenderer; // A volume renderer to render density distributions
	GLColorMap* densityPalette; // Color map to render density distributions
	VolumeRenderer::Vector viewDirection; // Current view direction
	
	/* Interaction state: */
	UndoBuffer* undoBuffer; // An undo buffer shared by all polypeptide chains
	Chain* activeChain; // The polypeptide chain being manipulated
	AutoBondingMode autoBondingMode; // Current auto-bonding mode
	UpdateDirection updateDirection; // Current update direction
	Cluster::MulticastPipe* interactorPipe; // Multicast pipe to synchronize IK updates across a distributed rendering cluster
	Threads::Thread ikThread; // Thread running IK calculations
	Chain* ikChain; // Polypeptide chain on which to perform IK
	ProteinInteractor::Transformation ikGoalTransformation; // Goal transformation for the most recent IK request
	Threads::MutexCond ikRequestCond; // Condition variable to wake up the IK thread
	volatile int ikRequest; // Counter for IK requests, updated by main thread
	int ikLastRequest; // Last IK request issued by a dragging tool
	volatile int ikResult; // Counter for IK results, updated by IK thread
	volatile bool ikKeepRunning; // Flag to shut down the IK thread
	
	/* Vrui stuff: */
	GLMotif::Button* undoButton;
	GLMotif::Button* redoButton;
	GLMotif::PopupMenu* mainMenu; // The program's main menu
	GLMotif::PopupWindow* chainDialog; // The polypeptide chain selection dialog
	GLMotif::ListBox* chainListBox;
	GLMotif::ToggleButton* chainVisibleToggle;
	GLMotif::ToggleButton* chainDrawAtomsToggle;
	GLMotif::ToggleButton* chainDrawBondsToggle;
	GLMotif::ToggleButton* chainDrawBackboneToggle;
	GLMotif::ToggleButton* chainDrawCartoonToggle;
	GLMotif::ToggleButton* chainDrawHydrogenBondsToggle;
	GLMotif::ToggleButton* chainDrawHydrogenBondSitesToggle;
	GLMotif::ToggleButton* chainDrawHydrogenCagesToggle;
	GLMotif::ToggleButton* chainDrawAtomCollisionsToggle;
	GLMotif::Button* resetTransformButton;
	GLMotif::PopupWindow* residueDialog; // The residue dialog
	GLMotif::TextFieldSlider* residueIndexSlider;
	GLMotif::TextField* residueTypeTextField;
	GLMotif::RadioBox* structureTypeBox;
	GLMotif::PopupWindow* interactionDialog; // The interaction dialog
	GLMotif::RadioBox* autoBondingBox;
	
	/* Private methods: */
	void finishIKSequence(Chain* ikChain); // Finishes an IK operation on the given polypeptide chain
	void undoCallback(Misc::CallbackData* cbData);
	void redoCallback(Misc::CallbackData* cbData);
	void showChainDialogCallback(Misc::CallbackData* cbData);
	void showResidueDialogCallback(Misc::CallbackData* cbData);
	void showInteractionDialogCallback(Misc::CallbackData* cbData);
	GLMotif::PopupMenu* createMainMenu(void);
	void updateChainDialogToggles(void);
	void updateChainDialog(void);
	void chainListBoxValueChangedCallback(GLMotif::ListBox::ValueChangedCallbackData* cbData);
	void chainVisibleToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawAtomsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawBondsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawBackboneToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawCartoonToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawHydrogenBondsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawHydrogenBondSitesToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawHydrogenCagesToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void chainDrawAtomCollisionsToggleCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void resetTransformCallback(Misc::CallbackData* cbData);
	GLMotif::PopupWindow* createChainDialog(void);
	
	void updateResidueFields(MD::Protein::Residue* residue);
	void updateResidueDialog(void);
	void residueIndexSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData);
	void structureTypeBoxCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData);
	GLMotif::PopupWindow* createResidueDialog(void);
	
	void autoBondingValueChangedCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData);
	void resetActiveCoilRegionsCallback(Misc::CallbackData* cbData);
	void ikUpdateDirectionValueChangedCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData);
	GLMotif::PopupWindow* createInteractionDialog(void);
	void* ikThreadMethod(void);
	
	/* Constructors and destructors: */
	public:
	VRProtoShop(int& argc,char**& argv);
	virtual ~VRProtoShop(void);
	
	/* Methods from Vrui::Application: */
	virtual void frame(void);
	virtual void display(GLContextData& contextData) const;
	virtual void eventCallback(EventID eventId,Vrui::InputDevice::ButtonCallbackData* cbData);
	virtual void resetNavigation(void);
	};

#endif
