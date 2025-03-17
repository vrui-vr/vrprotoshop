/***********************************************************************
ProteinRenderer - Class to visualize protein structures.
Copyright (c) 2002-2021 Oliver Kreylos

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

#ifndef PROTEINRENDERER_INCLUDED
#define PROTEINRENDERER_INCLUDED

#define USE_BONDCYLINDER_VERTEXARRAY 0

#include <utility>
#include <Misc/ConfigurationFile.h>
#include <Geometry/Box.h>
#include <GL/gl.h>
#include <GL/GLVertex.h>
#include <GL/GLMaterial.h>
#include <GL/GLColorMap.h>
#include <GL/GLObject.h>
#include <GL/GLLineIlluminator.h>
#include <GL/GLSphereRenderer.h>
#include <GL/GLCylinderRenderer.h>

#include "Protein.h"
#include "CatmullRomSpline.h"

namespace MD {

class ProteinRenderer:public GLObject
	{
	/* Embedded classes: */
	public:
	typedef GLMaterial::Color Color; // Data type for colors
	
	private:
	struct CollisionCell;
	
	struct AtomListItem // Structure associating atoms with collision detection cells
		{
		/* Elements: */
		public:
		const Protein::ChainAtom* atom; // Pointer to atom
		CollisionCell* cell; // Pointer to the cell containing the atom
		AtomListItem* cellSucc; // Pointer to next atom in the same cell
		};
	
	struct CollisionCell // Structure for collision detection cells
		{
		/* Elements: */
		public:
		AtomListItem* atoms; // Pointer to list of atoms inside this cell
		
		/* Constructors and destructors: */
		CollisionCell(void) // Creates empty cell
			:atoms(0)
			{
			};
		};
	
	struct StructureFlags // Structure holding rendering flags for each secondary structure
		{
		/* Elements: */
		public:
		const Protein::SecondaryStructure* structure; // The secondary structure
		unsigned int settingsVersion; // Version number of per-structure settings
		bool drawAtoms;
		bool drawBonds;
		GLfloat bondWidth; // Cosmetic bond (line) width
		Color bondColor;
		bool drawBackbone;
		GLfloat backboneWidth;
		Color backboneColor;
		bool drawBackboneRibbon;
		Scalar backboneRibbonWidth;
		int backboneRibbonNumSamples;
		Scalar backboneRibbonParameterMin,backboneRibbonParameterMax;
		bool drawCartoon;
		int numCartoonSplineSamples;
		Point* cartoonP;
		Vector* cartoonX;
		Vector* cartoonY;
		Vector* cartoonZ;
		bool drawHydrogenBondSites;
		GLfloat hydrogenBondSiteDiameter;
		GLfloat hydrogenBondSiteWidth;
		Color amideColor,carboxylColor;
		bool drawHydrogenCages;
		GLfloat hydrogenCageWidth;
		Color hydrogenCageColor;
		bool hydrogenCageLarge;
		
		/* Constructors and destructors: */
		StructureFlags(const Misc::ConfigurationFileSection& configFileSection,const Protein::SecondaryStructure* sStructure,int ri0,int ri1); // Creates default rendering flags
		~StructureFlags(void);
		};
	
	struct DataItem:public GLObject::DataItem
		{
		/* Elements: */
		public:
		GLuint hydrogenCageSmallDisplayListId; // Display list for small hydrogen cages
		GLuint hydrogenCageLargeDisplayListId; // Display list for large hydrogen cages
		unsigned int globalSettingsVersion; // Version number of global rendering settings
		unsigned int secondaryStructureListVersion; // Version number of secondary structure list
		int numSecondaryStructures; // Number of secondary structures in list
		GLuint secondaryStructureListBaseId; // Base of display lists for secondary structures
		unsigned int* secondaryStructureSettingsVersions; // Version numbers of per-structure rendering settings
		unsigned int* secondaryStructureVersions; // Version numbers of secondary structures
		
		/* Constructors and destructors: */
		DataItem(void);
		virtual ~DataItem(void);
		};
	
	/* Elements: */
	Misc::ConfigurationFileSection configFileSection; // Configuration file section containing rendering parameters
	const Protein* protein; // Pointer to the visualized protein
	AtomListItem* atomListItems; // C-style array of atom list items; stays valid throughout lifetime of renderer
	double collisionCellSize[3];
	int numCollisionCells[3],collisionCellsWidth[3];
	CollisionCell* collisionCells;
	CollisionCell* collisionCellBase;
	static const double maxAtomRadius; // Maximum radius of any atom sphere
	Geometry::Box<Scalar,3> boundingBox; // Protein's current bounding box
	std::vector<StructureFlags> structureFlags; // Rendering flags for each secondary structure
	unsigned int globalSettingsVersion; // Version number for global rendering settings
	bool drawAtoms; // Global flag for atom rendering
	GLMaterial atomMaterial; // Material for atom rendering
	GLSphereRenderer atomRenderer; // Renderer for atoms
	int atomTesselation; // Number of subdivisions for atom spheres
	static const Color elementColors[118]; // Standard colors to render atoms
	bool mapAtomValues; // Flag to enable value mapping for atoms
	GLColorMap atomColorMap; // Color map for mapping atom values
	float* atomValues; // Array of values that can be used to colormap atoms
	bool drawBonds; // Global flag for bond rendering
	GLLineIlluminator bondIlluminator;
	GLMaterial bondMaterial;
	unsigned int bondMaterialVersion;
	GLSphereRenderer bondAtomRenderer; // Renderer for atoms in bond rendering
	GLCylinderRenderer bondRenderer; // Renderer for bonds
	int numBondVertices; // Number of vertices for fancy stick model rendering
	GLfloat bondRadius; // Cylinder radius for fancy stick model rendering
	Color alphaHelixColor; // Color for alpha helices
	Color betaStrandColor; // Color for beta strands
	Color coilColor; // Color for coils
	Color highlightColor; // Color for selected structures
	bool drawBackbone; // Global flag for backbone rendering
	bool drawBackboneRibbon; // Global flag for backbone ribbon rendering
	GLMaterial backboneRibbonMaterial;
	bool backboneRibbonUseAllAtoms; // Flag if backbone spline uses all backbone atoms or only one per residue as control points
	int backboneRibbonDegree; // B-spline degree of backbone ribbon
	int backboneRibbonSampleDensity; // Sampling density for rendering the backbone ribbon
	CatmullRomSpline<Scalar,3> backboneRibbonSpline1,backboneRibbonSpline2; // Major and minor spline curves used to render the backbone
	bool drawCartoon; // Global flag for cartoon rendering
	GLMaterial cartoonMaterial;
	int cartoonSampleDensity; // Sampling density for rendering cartoon splines
	Scalar alphaHelixWidth,alphaHelixThickness; // Width and thickness of alpha helix strip
	Scalar betaStrandWidth,betaStrandThickness; // Width and thickness of beta strand arrow
	Scalar betaStrandHeadWidth; // Scale factor for width of beta strand arrow heads
	int numCoilVertices; // Number of vertices for coil tubes
	Scalar coilRadius; // Radius of coil tubes
	bool drawHydrogenBonds; // Global flag for hydrogen bond rendering
	GLfloat hydrogenBondWidth;
	Color hydrogenBondColor;
	bool drawHydrogenBondSites; // Global flag for hydrogen bond site rendering
	bool drawHydrogenCages; // Global flag for hydrogen cage rendering
	bool drawCollisions; // Global flag for collision sphere rendering
	GLMaterial collisionSphereMaterial; // Material for collision spheres
	GLSphereRenderer collisionSphereRenderer; // Renderer for collision spheres
	int collisionSphereTesselation; // Number of subdivisions for collision spheres
	GLfloat modelScale; // Current scaling factor of the modelview transformation
	
	/* Private methods: */
	void createBackboneRibbonSpline(void);
	void createCartoonSpline(StructureFlags& sf);
	void updateCartoonSpline(StructureFlags& sf);
	void glDrawBonds(const StructureFlags& sf,DataItem& dataItem) const;
	void glDrawBackbone(const StructureFlags& sf,DataItem& dataItem) const;
	void glDrawBackboneRibbon(const StructureFlags& sf,DataItem& dataItem) const;
	void glDrawCartoon(const StructureFlags& sf,DataItem& dataItem) const;
	void glDrawHydrogenBondSites(const StructureFlags& sf,DataItem& dataItem) const;
	void glDrawHydrogenCages(const StructureFlags& sf,DataItem& dataItem) const;
	void glDrawHydrogenBonds(GLContextData& contextData) const;
	void glDrawCollisions(GLContextData& contextData) const;
	
	/* Constructors and destructors: */
	public:
	ProteinRenderer(const Misc::ConfigurationFileSection& sConfigFileSection,const Protein* sProtein); // Creates a renderer for the given protein
	virtual ~ProteinRenderer(void);
	
	/* Methods: */
	virtual void initContext(GLContextData& contextData) const;
	void updateStructureFlags(void);  // Updates structure rendering flags after change to protein's secondary structure sequence
	void updateProtein(void); // Updates internal representation of the protein structure after changes
	void setModelScale(GLfloat newModelScale); // Updates the scale factor of the modelview matrix for subsequent rendering
	template <class TransformationParam>
	std::pair<double,double> calcDepthRange(const TransformationParam& modelView) const // Returns depth value range for given modelview transformation
		{
		double frontDist,backDist;
		frontDist=backDist=-modelView.transform(boundingBox.getVertex(0))[2];
		for(int i=1;i<8;++i)
			{
			double dist=-modelView.transform(boundingBox.getVertex(i))[2];
			if(frontDist>dist)
				frontDist=dist;
			else if(backDist<dist)
				backDist=dist;
			}
		return std::pair<double,double>(frontDist-maxAtomRadius,backDist+maxAtomRadius);
		}
	void setViewDirection(const GLLineIlluminator::Vector& viewDirection); // Sets the view direction for rendering
	void setLightDirection(const GLLineIlluminator::Vector& lightDirection); // Sets the light direction for rendering
	void glRenderAction(GLContextData& contextData) const; // Renders the protein
	void highlightResidue(GLContextData& contextData,const Protein::Residue* rPtr) const; // Highlights a single residue
	int glPick(GLContextData& contextData,int which) const; // Uses OpenGL picking to return the index of the secondary structure (which==0) or residue (which==1) hit by the query point; -1 if none found
	Color getResidueBackboneColor(const Protein::Residue* rPtr) const; // Returns the color a residue's ba
	void lockResidueRange(int firstResidueIndex,int numResidues,bool locked); // Marks a range of residues (must coincide with structure boundaries) as locked or unlocked
	
	/* Parameter access methods: */
	bool getDrawAtoms(void) const
		{
		return drawAtoms;
		}
	void setDrawAtoms(bool newDrawAtoms);
	bool getMapAtomValues(void) const
		{
		return mapAtomValues;
		}
	void setMapAtomValues(bool newMapAtomValues);
	void setMapAtomValueRange(float newMinValue,float newMaxValue) // Sets mapping range for atom values
		{
		atomColorMap.setScalarRange(newMinValue,newMaxValue);
		}
	void setAtomValue(int atomIndex,float value) // Sets colormapping source value for an atom
		{
		atomValues[atomIndex]=value;
		}
	bool getDrawBonds(void) const
		{
		return drawBonds;
		}
	void setDrawBonds(bool newDrawBonds);
	bool getDrawBackbone(void) const
		{
		return drawBackbone;
		}
	void setDrawBackbone(bool newDrawBackbone);
	bool getDrawBackboneRibbon(void) const
		{
		return drawBackboneRibbon;
		}
	void setDrawBackboneRibbon(bool newDrawBackboneRibbon);
	bool getDrawCartoon(void) const
		{
		return drawCartoon;
		}
	void setDrawCartoon(bool newDrawCartoon);
	bool getDrawHydrogenBonds(void) const
		{
		return drawHydrogenBonds;
		}
	void setDrawHydrogenBonds(bool newDrawHydrogenBonds);
	bool getDrawHydrogenBondSites(void) const
		{
		return drawHydrogenBondSites;
		}
	void setDrawHydrogenBondSites(bool newDrawHydrogenBondSites);
	bool getDrawHydrogenCages(void) const
		{
		return drawHydrogenCages;
		}
	void setDrawHydrogenCages(bool newDrawHydrogenCages);
	bool getDrawCollisions(void) const
		{
		return drawCollisions;
		}
	void setDrawCollisions(bool newDrawCollisions);
	
	/* Per-structure parameter access methods: */
	bool getDrawAtoms(const Protein::StructureSelector& selector) const;
	void setDrawAtoms(const Protein::StructureSelector& selector,bool newDrawAtoms);
	bool getDrawBonds(const Protein::StructureSelector& selector) const;
	void setDrawBonds(const Protein::StructureSelector& selector,bool newDrawBonds);
	bool getDrawBackbone(const Protein::StructureSelector& selector) const;
	void setDrawBackbone(const Protein::StructureSelector& selector,bool newDrawBackbone);
	bool getDrawBackboneRibbon(const Protein::StructureSelector& selector) const;
	void setDrawBackboneRibbon(const Protein::StructureSelector& selector,bool newDrawBackboneRibbon);
	bool getDrawCartoon(const Protein::StructureSelector& selector) const;
	void setDrawCartoon(const Protein::StructureSelector& selector,bool newDrawCartoon);
	void resetAllBackboneColor(void);
	void resetBackboneColor(const Protein::StructureSelector& selector);
	void setBackboneColor(const Protein::StructureSelector& selector,const Color& newBackboneColor);
	void selectStructure(const Protein::StructureSelector& selector)
		{
		setBackboneColor(selector,highlightColor);
		}
	void deselectStructure(const Protein::StructureSelector& selector)
		{
		resetBackboneColor(selector);
		}
	bool getDrawHydrogenBondSites(const Protein::StructureSelector& selector) const;
	void setAllDrawHydrogenBondSites(bool newDrawHydrogenBondSites);
	void setDrawHydrogenBondSites(const Protein::StructureSelector& selector,bool newDrawHydrogenBondSites);
	bool getDrawHydrogenCages(const Protein::StructureSelector& selector) const;
	void setDrawHydrogenCages(const Protein::StructureSelector& selector,bool newDrawHydrogenCages);
	bool getDrawLargeHydrogenCages(const Protein::StructureSelector& selector) const;
	void setDrawLargeHydrogenCages(const Protein::StructureSelector& selector,bool newDrawLargeHydrogenCages);
	};

}

#endif
