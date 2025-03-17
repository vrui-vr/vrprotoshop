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

#include "ProteinRenderer.h"

#include <Misc/StandardValueCoders.h>
#include <Math/Constants.h>
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLValueCoders.h>
#include <GL/GLContextData.h>
#include <GL/GLModels.h>
#include <GL/GLGeometryWrappers.h>

#define USE_GL_PICKING 0

namespace MD {

/************************************************
Methods of class ProteinRenderer::StructureFlags:
************************************************/

ProteinRenderer::StructureFlags::StructureFlags(const Misc::ConfigurationFileSection& configFileSection,const Protein::SecondaryStructure* sStructure,int ri0,int ri1)
	:structure(sStructure),
	 settingsVersion(1),
	 drawAtoms(configFileSection.retrieveValue<bool>("./drawAtoms",true)),
	 drawBonds(configFileSection.retrieveValue<bool>("./drawBonds",true)),
	 bondWidth(configFileSection.retrieveValue<float>("./bondWidth",2.0f)),
	 bondColor(configFileSection.retrieveValue<Color>("./bondColor",Color(1.0,0.0,0.0))),
	 drawBackbone(configFileSection.retrieveValue<bool>("./drawBackbone",true)),
	 backboneWidth(configFileSection.retrieveValue<float>("./backboneWidth",4.0f)),
	 backboneColor(0.5f,0.5f,0.5f),
	 drawBackboneRibbon(configFileSection.retrieveValue<bool>("./drawBackboneRibbon",true)),
	 backboneRibbonWidth(configFileSection.retrieveValue<Scalar>("./backboneRibbonWidth",Scalar(2))),
	 backboneRibbonNumSamples((ri1-ri0)*5+1),backboneRibbonParameterMin(ri0),backboneRibbonParameterMax(ri1),
	 drawCartoon(configFileSection.retrieveValue<bool>("./drawCartoon",true)),
	 cartoonP(0),cartoonX(0),cartoonY(0),cartoonZ(0),
	 drawHydrogenBondSites(configFileSection.retrieveValue<bool>("./drawHydrogenBondSites",structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND)),
	 hydrogenBondSiteDiameter(configFileSection.retrieveValue<float>("./hydrogenBondSiteDiameter",5.0f)),
	 hydrogenBondSiteWidth(configFileSection.retrieveValue<float>("./hydrogenBondSiteWidth",1.0f)),
	 amideColor(configFileSection.retrieveValue<Color>("./amideColor",Color(0.5,0.5,1.0))),
	 carboxylColor(configFileSection.retrieveValue<Color>("./carboxylColor",Color(1.0,1.0,0.0))),
	 drawHydrogenCages(configFileSection.retrieveValue<bool>("./drawHydrogenCages",structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND)),
	 hydrogenCageWidth(configFileSection.retrieveValue<float>("./hydrogenCageWidth",1.5f)),
	 hydrogenCageColor(configFileSection.retrieveValue<Color>("./hydrogenCageColor",Color(0.6,0.5,0.0))),
	 hydrogenCageLarge(configFileSection.retrieveValue<bool>("./hydrogenCageLarge",structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND))
	{
	}

ProteinRenderer::StructureFlags::~StructureFlags(void)
	{
	delete[] cartoonP;
	delete[] cartoonX;
	delete[] cartoonY;
	delete[] cartoonZ;
	}

/******************************************
Methods of class ProteinRenderer::DataItem:
******************************************/

ProteinRenderer::DataItem::DataItem(void)
	:hydrogenCageSmallDisplayListId(glGenLists(2)),
	 hydrogenCageLargeDisplayListId(hydrogenCageSmallDisplayListId+1),
	 globalSettingsVersion(0),
	 secondaryStructureListVersion(0),
	 numSecondaryStructures(0),
	 secondaryStructureListBaseId(0),
	 secondaryStructureSettingsVersions(0),
	 secondaryStructureVersions(0)
	{
	}

ProteinRenderer::DataItem::~DataItem(void)
	{
	glDeleteLists(hydrogenCageSmallDisplayListId,2);
	if(numSecondaryStructures>0)
		{
		glDeleteLists(secondaryStructureListBaseId,numSecondaryStructures);
		delete[] secondaryStructureSettingsVersions;
		delete[] secondaryStructureVersions;
		}
	}

/****************************************
Static elements of class ProteinRenderer:
****************************************/

const double ProteinRenderer::maxAtomRadius=1.8;
const ProteinRenderer::Color ProteinRenderer::elementColors[118]=
	{
	Color(0.8f,0.8f,0.8f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.2f,0.2f,0.2f),Color(0.5f,0.5f,1.0f),Color(1.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,1.0f,0.0f),Color(1.0f,1.0f,0.0f),
	Color(0.0f,1.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(1.0f,0.5f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
	Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f)
	};

/********************************
Methods of class ProteinRenderer:
********************************/

void ProteinRenderer::createBackboneRibbonSpline(void)
	{
	/* Create the major and minor backbone ribbon splines: */
	size_t numResidues=protein->getNumResidues();
	backboneRibbonSpline1.setNumControlPoints(numResidues+1);
	CatmullRomSpline<Scalar,3>& s1=backboneRibbonSpline1;
	backboneRibbonSpline2.setNumControlPoints(numResidues+1);
	CatmullRomSpline<Scalar,3>& s2=backboneRibbonSpline2;
	
	Scalar tanScale=Scalar(1)/Scalar(3);
	
	/* Calculate major backbone ribbon spline control point positions: */
	const Protein::Residue* rPtr0=protein->residues;
	s1[0].position=rPtr0->backbone.front()->getPosition();
	size_t i=1;
	const Protein::Residue* rPtr1;
	for(rPtr1=rPtr0->succ;rPtr1!=0;rPtr0=rPtr1,rPtr1=rPtr1->succ,++i)
		s1[i].position=Geometry::mid(rPtr0->backbone.back()->getPosition(),rPtr1->backbone.front()->getPosition());
	s1[i].position=rPtr0->backbone.back()->getPosition();
	
	/* Calculate major backbone ribbon spline control point tangents: */
	if(numResidues==1)
		{
		s1[1].tangent=s1[0].tangent=(s1[1].position-s1[0].position)*tanScale;
		}
	else
		{
		Vector d1,d2;
		d1=s1[1].position-s1[0].position;
		d2=s1[2].position-s1[0].position;
		s1[0].tangent=d1*(Scalar(2)*tanScale)-d2*(Scalar(0.5)*tanScale);
		for(size_t i=1;i<numResidues;++i)
			s1[i].tangent=(s1[i+1].position-s1[i-1].position)*(Scalar(0.5)*tanScale);
		d1=s1[numResidues].position-s1[numResidues-1].position;
		d2=s1[numResidues].position-s1[numResidues-2].position;
		s1[numResidues].tangent=d1*(Scalar(2)*tanScale)-d2*(Scalar(0.5)*tanScale);
		}
	
	/* Calculate minor backbone ribbon spline control point positions: */
	rPtr0=protein->residues;
	Vector d0=rPtr0->backbone[1]->getPosition()-rPtr0->backbone[0]->getPosition();
	d0.orthogonalize(s1[0].tangent);
	s2[0].position=s1[0].position+d0/d0.mag();
	i=1;
	for(rPtr1=rPtr0->succ;rPtr1!=0;rPtr0=rPtr1,rPtr1=rPtr1->succ,++i)
		{
		Vector d1=rPtr1->backbone.front()->getPosition()-rPtr0->backbone.back()->getPosition();
		d1.orthogonalize(s1[i].tangent);
		if(d0*d1<Scalar(0))
			d1=-d1;
		s2[i].position=s1[i].position+d1/d1.mag();
		d0=d1;
		}
	size_t n=rPtr0->backbone.size();
	Vector d1=rPtr0->backbone[n-1]->getPosition()-rPtr0->backbone[n-2]->getPosition();
	d1.orthogonalize(s1[i].tangent);
	if(d0*d1<Scalar(0))
		d1=-d1;
	s2[i].position=s1[i].position+d1/d1.mag();
	
	/* Calculate minor backbone ribbon spline control point tangents: */
	if(numResidues==1)
		{
		s2[1].tangent=s2[0].tangent=(s2[1].position-s2[0].position)*tanScale;
		}
	else
		{
		Vector d1,d2;
		d1=s2[1].position-s2[0].position;
		d2=s2[2].position-s2[0].position;
		s2[0].tangent=d1*(Scalar(2)*tanScale)-d2*(Scalar(0.5)*tanScale);
		for(size_t i=1;i<numResidues;++i)
			s2[i].tangent=(s2[i+1].position-s2[i-1].position)*(Scalar(0.5)*tanScale);
		d1=s2[numResidues].position-s2[numResidues-1].position;
		d2=s2[numResidues].position-s2[numResidues-2].position;
		s2[numResidues].tangent=d1*(Scalar(2)*tanScale)-d2*(Scalar(0.5)*tanScale);
		}
	}

void ProteinRenderer::createCartoonSpline(ProteinRenderer::StructureFlags& sf)
	{
	/* Create evaluation arrays for the cartoon rendering splines: */
	sf.numCartoonSplineSamples=sf.structure->getNumResidues()*cartoonSampleDensity;
	delete sf.cartoonP;
	sf.cartoonP=new Point[sf.numCartoonSplineSamples+1];
	delete sf.cartoonX;
	sf.cartoonX=new Vector[sf.numCartoonSplineSamples+1];
	delete sf.cartoonY;
	sf.cartoonY=new Vector[sf.numCartoonSplineSamples+1];
	delete sf.cartoonZ;
	sf.cartoonZ=new Vector[sf.numCartoonSplineSamples+1];
	}

void ProteinRenderer::updateCartoonSpline(ProteinRenderer::StructureFlags& sf)
	{
	if(sf.drawCartoon)
		{
		/* Evaluate cartoon splines based on secondary structure type: */
		switch(sf.structure->getStructureType())
			{
			case Protein::SecondaryStructure::COIL:
				{
				/* Evaluate the minor spline once: */
				Vector x=backboneRibbonSpline2(sf.backboneRibbonParameterMin)-backboneRibbonSpline1(sf.backboneRibbonParameterMin);
				for(int i=0;i<=sf.numCartoonSplineSamples;++i)
					{
					/* Evaluate the major spline: */
					Scalar u=Scalar(i)*(sf.backboneRibbonParameterMax-sf.backboneRibbonParameterMin)/Scalar(sf.numCartoonSplineSamples)+sf.backboneRibbonParameterMin;
					Point c=backboneRibbonSpline1(u);
					Vector t=backboneRibbonSpline1.diff(u);
					
					/* Construct a coordinate frame around the major spline point: */
					sf.cartoonP[i]=c;
					x.orthogonalize(t);
					sf.cartoonX[i]=x;
					sf.cartoonX[i].normalize();
					sf.cartoonY[i]=t;
					sf.cartoonY[i].normalize();
					sf.cartoonZ[i]=sf.cartoonX[i]^sf.cartoonY[i];
					sf.cartoonZ[i].normalize();
					}
				break;
				}
			
			case Protein::SecondaryStructure::ALPHA_HELIX:
				{
				for(int i=0;i<=sf.numCartoonSplineSamples;++i)
					{
					/* Evaluate the major and minor splines: */
					Scalar u=Scalar(i)*(sf.backboneRibbonParameterMax-sf.backboneRibbonParameterMin)/Scalar(sf.numCartoonSplineSamples)+sf.backboneRibbonParameterMin;
					Point c=backboneRibbonSpline1(u);
					Vector t=backboneRibbonSpline1.diff(u);
					Point s=backboneRibbonSpline2(u);
					
					/* Construct a coordinate frame around the major spline point: */
					sf.cartoonP[i]=c;
					sf.cartoonX[i]=s-c;
					sf.cartoonX[i].normalize();
					sf.cartoonY[i]=t;
					sf.cartoonY[i].normalize();
					sf.cartoonZ[i]=sf.cartoonX[i]^sf.cartoonY[i];
					sf.cartoonZ[i].normalize();
					}
				
				break;
				}
			
			case Protein::SecondaryStructure::BETA_STRAND:
				{
				Scalar alpha=Math::rad(Scalar(-22.5));
				Scalar ca=Math::cos(alpha);
				Scalar sa=Math::sin(alpha);
				for(int i=0;i<=sf.numCartoonSplineSamples;++i)
					{
					/* Evaluate the major and minor splines: */
					Scalar u=Scalar(i)*(sf.backboneRibbonParameterMax-sf.backboneRibbonParameterMin)/Scalar(sf.numCartoonSplineSamples)+sf.backboneRibbonParameterMin;
					Point c=backboneRibbonSpline1(u);
					Vector t=backboneRibbonSpline1.diff(u);
					Point s=backboneRibbonSpline2(u);
					
					/* Construct a rotated coordinate frame around the major spline point: */
					sf.cartoonP[i]=c;
					Vector x=s-c;
					x.normalize();
					Vector z=x^t;
					z.normalize();
					sf.cartoonX[i]=x*ca+z*sa;
					sf.cartoonY[i]=t;
					sf.cartoonY[i].normalize();
					sf.cartoonZ[i]=z*ca-x*sa;
					}
				
				break;
				}
			
			default:
				; // Just to make compiler happy
			}
		}
	}

void ProteinRenderer::glDrawBonds(const ProteinRenderer::StructureFlags& sf,ProteinRenderer::DataItem& dataItem) const
	{
	/* Set rendering parameters for this structure: */
	glLineWidth(sf.bondWidth);
	// glColor(sf.bondColor);
	glColor(1.0f,1.0f,1.0f);

	/* Iterate through all residues in this secondary structure: */
	for(const Protein::Residue* rPtr=sf.structure->residueBegin;rPtr!=sf.structure->residueEnd;rPtr=rPtr->succ)
		{
		#if USE_GL_PICKING
		/* Push residue (PDB) index onto the selection name stack: */
		glPushName(rPtr->residueIndex);
		#endif
		
		/* Iterate through all atoms in this residue: */
		glBegin(GL_LINES);
		for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
			{
			/* Iterate through all bonds for this atom: */
			for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
				{
				if(static_cast<const Atom*>(aPtr)<*bIt)
					{
					/* Calculate bond axis for illumination: */
					Vector bondAxis=(*bIt)->getPosition()-aPtr->getPosition();
					glTexCoord(bondAxis.normalize());
					glVertex(aPtr->getPosition());
					glVertex((*bIt)->getPosition());
					}
				}
			}
		glEnd();
		
		#if USE_GL_PICKING
		glPopName();
		#endif
		}
	}

void ProteinRenderer::glDrawBackbone(const ProteinRenderer::StructureFlags& sf,ProteinRenderer::DataItem& dataItem) const
	{
	/* Set rendering parameters for this structure: */
	glLineWidth(sf.backboneWidth);
	glColor(sf.backboneColor);
	
	/* Render residues in this structure independently: */
	for(const Protein::Residue* rPtr=sf.structure->residueBegin;rPtr!=sf.structure->residueEnd;rPtr=rPtr->succ)
		{
		#if USE_GL_PICKING
		/* Push residue (PDB) index onto the selection name stack: */
		glPushName(rPtr->residueIndex);
		#endif
		
		glBegin(GL_LINE_STRIP);
		if(rPtr->pred!=0)
			{
			/* Draw midpoint between last atom in previous residue and first atom in this one: */
			glVertex(Geometry::mid(rPtr->pred->backbone[rPtr->pred->backbone.size()-1]->getPosition(),rPtr->backbone[0]->getPosition()));
			}
		for(std::vector<Protein::ChainAtom*>::const_iterator bbIt=rPtr->backbone.begin();bbIt!=rPtr->backbone.end();++bbIt)
			glVertex((*bbIt)->getPosition());
		if(rPtr->succ!=0)
			{
			/* Draw midpoint between last atom in this residue and first atom in the next one: */
			glVertex(Geometry::mid(rPtr->backbone[rPtr->backbone.size()-1]->getPosition(),rPtr->succ->backbone[0]->getPosition()));
			}
		glEnd();
		
		#if USE_GL_PICKING
		glPopName();
		#endif
		}
	
	/* Reset OpenGL state: */
	}

void ProteinRenderer::glDrawBackboneRibbon(const ProteinRenderer::StructureFlags& sf,ProteinRenderer::DataItem& dataItem) const
	{
	/* Set ribbon rendering parameters: */
	glDisable(GL_CULL_FACE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
	glMaterial(GLMaterialEnums::FRONT_AND_BACK,backboneRibbonMaterial);
	
	/* Set rendering parameters for this structure: */
	glColor(sf.backboneColor);
	glBegin(GL_QUAD_STRIP);
	
	/* Render all samples: */
	Scalar sideFactor=sf.backboneRibbonWidth*Scalar(0.5);
	for(int i=0;i<sf.backboneRibbonNumSamples;++i)
		{
		/* Evaluate the spline curve: */
		Scalar u=Scalar(i)*(sf.backboneRibbonParameterMax-sf.backboneRibbonParameterMin)/Scalar(sf.backboneRibbonNumSamples-1)+sf.backboneRibbonParameterMin;
		Point c=backboneRibbonSpline1(u);
		Vector d=backboneRibbonSpline1.diff(u);
		Vector s=backboneRibbonSpline2(u)-c;
		s*=sideFactor/s.mag();
		
		/* Draw a quad on the strip: */
		glNormal(s^d);
		glVertex(c-s);
		glVertex(c+s);
		}
	
	glEnd();
	
	/* Reset OpenGL state: */
	glEnable(GL_CULL_FACE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
	glDisable(GL_COLOR_MATERIAL);
	}

void ProteinRenderer::glDrawCartoon(const ProteinRenderer::StructureFlags& sf,ProteinRenderer::DataItem& dataItem) const
	{
	/* Set cartoon rendering parameters: */
	glMaterial(GLMaterialEnums::FRONT,cartoonMaterial);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT,GL_DIFFUSE);
	
	switch(sf.structure->getStructureType())
		{
		case Protein::SecondaryStructure::COIL:
			/* Render coil region: */
			glColor(sf.backboneColor);
			for(int i=0;i<numCoilVertices;++i)
				{
				Scalar angle1=Scalar(i)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
				Scalar xFac1=Math::cos(angle1);
				Scalar zFac1=Math::sin(angle1);
				Scalar angle2=Scalar(i+1)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
				Scalar xFac2=Math::cos(angle2);
				Scalar zFac2=Math::sin(angle2);
				glBegin(GL_QUAD_STRIP);
				for(int j=0;j<=sf.numCartoonSplineSamples;++j)
					{
					Vector v2=sf.cartoonX[j]*xFac2+sf.cartoonZ[j]*zFac2;
					glNormal(v2);
					glVertex(sf.cartoonP[j]+v2*coilRadius);
					Vector v1=sf.cartoonX[j]*xFac1+sf.cartoonZ[j]*zFac1;
					glNormal(v1);
					glVertex(sf.cartoonP[j]+v1*coilRadius);
					}
				glEnd();
				}
			break;
		
		case Protein::SecondaryStructure::ALPHA_HELIX:
			{
			/* Render alpha helix: */
			glColor(sf.backboneColor);
			glBegin(GL_QUADS);
			Vector x0=sf.cartoonX[0]*alphaHelixWidth;
			Vector z0=sf.cartoonZ[0]*alphaHelixThickness;
			glNormal(-sf.cartoonY[0]);
			glVertex(sf.cartoonP[0]-z0-x0);
			glVertex(sf.cartoonP[0]-z0+x0);
			glVertex(sf.cartoonP[0]+z0+x0);
			glVertex(sf.cartoonP[0]+z0-x0);
			Vector xn=sf.cartoonX[sf.numCartoonSplineSamples]*alphaHelixWidth;
			Vector zn=sf.cartoonZ[sf.numCartoonSplineSamples]*alphaHelixThickness;
			glNormal(sf.cartoonY[sf.numCartoonSplineSamples]);
			glVertex(sf.cartoonP[sf.numCartoonSplineSamples]-zn-xn);
			glVertex(sf.cartoonP[sf.numCartoonSplineSamples]+zn-xn);
			glVertex(sf.cartoonP[sf.numCartoonSplineSamples]+zn+xn);
			glVertex(sf.cartoonP[sf.numCartoonSplineSamples]-zn+xn);
			glEnd();
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=sf.numCartoonSplineSamples;++i)
				{
				Vector x=sf.cartoonX[i]*alphaHelixWidth;
				Vector z=sf.cartoonZ[i]*alphaHelixThickness;
				glNormal(sf.cartoonZ[i]);
				glVertex(sf.cartoonP[i]+z-x);
				glVertex(sf.cartoonP[i]+z+x);
				}
			glEnd();
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=sf.numCartoonSplineSamples;++i)
				{
				Vector x=sf.cartoonX[i]*alphaHelixWidth;
				Vector z=sf.cartoonZ[i]*alphaHelixThickness;
				glNormal(-sf.cartoonZ[i]);
				glVertex(sf.cartoonP[i]-z+x);
				glVertex(sf.cartoonP[i]-z-x);
				}
			glEnd();
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=sf.numCartoonSplineSamples;++i)
				{
				Vector x=sf.cartoonX[i]*alphaHelixWidth;
				Vector z=sf.cartoonZ[i]*alphaHelixThickness;
				glNormal(sf.cartoonX[i]);
				glVertex(sf.cartoonP[i]+z+x);
				glVertex(sf.cartoonP[i]-z+x);
				}
			glEnd();
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=sf.numCartoonSplineSamples;++i)
				{
				Vector x=sf.cartoonX[i]*alphaHelixWidth;
				Vector z=sf.cartoonZ[i]*alphaHelixThickness;
				glNormal(-sf.cartoonX[i]);
				glVertex(sf.cartoonP[i]-z-x);
				glVertex(sf.cartoonP[i]+z-x);
				}
			glEnd();
			break;
			}
		
		case Protein::SecondaryStructure::BETA_STRAND:
			{
			/* Render beta strand: */
			glColor(sf.backboneColor);
			glBegin(GL_QUADS);
			Vector x0=sf.cartoonX[0]*betaStrandWidth;
			Vector z0=sf.cartoonZ[0]*betaStrandThickness;
			glNormal(-sf.cartoonY[0]);
			glVertex(sf.cartoonP[0]-z0-x0);
			glVertex(sf.cartoonP[0]-z0+x0);
			glVertex(sf.cartoonP[0]+z0+x0);
			glVertex(sf.cartoonP[0]+z0-x0);
			glEnd();
			int tailEnd=sf.numCartoonSplineSamples-cartoonSampleDensity/2;
			Scalar widthFactor=betaStrandHeadWidth/Scalar(cartoonSampleDensity/2);
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=tailEnd;++i)
				{
				Vector x=sf.cartoonX[i]*betaStrandWidth;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				glNormal(sf.cartoonZ[i]);
				glVertex(sf.cartoonP[i]+z-x);
				glVertex(sf.cartoonP[i]+z+x);
				}
			for(int i=tailEnd;i<=sf.numCartoonSplineSamples;++i)
				{
				Scalar width=betaStrandWidth*Scalar(sf.numCartoonSplineSamples-i)*widthFactor;
				Vector x=sf.cartoonX[i]*width;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				glNormal(sf.cartoonZ[i]);
				glVertex(sf.cartoonP[i]+z-x);
				glVertex(sf.cartoonP[i]+z+x);
				}
			glEnd();
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=tailEnd;++i)
				{
				Vector x=sf.cartoonX[i]*betaStrandWidth;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				glNormal(-sf.cartoonZ[i]);
				glVertex(sf.cartoonP[i]-z+x);
				glVertex(sf.cartoonP[i]-z-x);
				}
			for(int i=tailEnd;i<=sf.numCartoonSplineSamples;++i)
				{
				Scalar width=betaStrandWidth*Scalar(sf.numCartoonSplineSamples-i)*widthFactor;
				Vector x=sf.cartoonX[i]*width;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				glNormal(-sf.cartoonZ[i]);
				glVertex(sf.cartoonP[i]-z+x);
				glVertex(sf.cartoonP[i]-z-x);
				}
			glEnd();
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=tailEnd;++i)
				{
				Vector x=sf.cartoonX[i]*betaStrandWidth;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				glNormal(sf.cartoonX[i]);
				glVertex(sf.cartoonP[i]+z+x);
				glVertex(sf.cartoonP[i]-z+x);
				if(i==tailEnd)
					{
					glNormal(-sf.cartoonY[i]);
					glVertex(sf.cartoonP[i]+z+x);
					glVertex(sf.cartoonP[i]-z+x);
					}
				}
			for(int i=tailEnd;i<=sf.numCartoonSplineSamples;++i)
				{
				Scalar width=betaStrandWidth*Scalar(sf.numCartoonSplineSamples-i)*widthFactor;
				Vector x=sf.cartoonX[i]*width;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				if(i==tailEnd)
					{
					glNormal(-sf.cartoonY[i]);
					glVertex(sf.cartoonP[i]+z+x);
					glVertex(sf.cartoonP[i]-z+x);
					}
				glNormal(sf.cartoonY[i]+sf.cartoonX[i]);
				glVertex(sf.cartoonP[i]+z+x);
				glVertex(sf.cartoonP[i]-z+x);
				}
			glEnd();
			glBegin(GL_QUAD_STRIP);
			for(int i=0;i<=tailEnd;++i)
				{
				Vector x=sf.cartoonX[i]*betaStrandWidth;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				glNormal(-sf.cartoonX[i]);
				glVertex(sf.cartoonP[i]-z-x);
				glVertex(sf.cartoonP[i]+z-x);
				if(i==tailEnd)
					{
					glNormal(-sf.cartoonY[i]);
					glVertex(sf.cartoonP[i]-z-x);
					glVertex(sf.cartoonP[i]+z-x);
					}
				}
			for(int i=tailEnd;i<=sf.numCartoonSplineSamples;++i)
				{
				Scalar width=betaStrandWidth*Scalar(sf.numCartoonSplineSamples-i)*widthFactor;
				Vector x=sf.cartoonX[i]*width;
				Vector z=sf.cartoonZ[i]*betaStrandThickness;
				if(i==tailEnd)
					{
					glNormal(-sf.cartoonY[i]);
					glVertex(sf.cartoonP[i]-z-x);
					glVertex(sf.cartoonP[i]+z-x);
					}
				glNormal(sf.cartoonY[i]-sf.cartoonX[i]);
				glVertex(sf.cartoonP[i]-z-x);
				glVertex(sf.cartoonP[i]+z-x);
				}
			glEnd();
			break;
			}
		
		default:
			; // Just to make compiler happy
		}
	
	/* Reset OpenGL state: */
	glDisable(GL_COLOR_MATERIAL);
	}

void ProteinRenderer::glDrawHydrogenBondSites(const ProteinRenderer::StructureFlags& sf,ProteinRenderer::DataItem& dataItem) const
	{
	/* Set rendering parameters for this structure: */
	glPointSize(sf.hydrogenBondSiteDiameter);
	glLineWidth(sf.hydrogenBondSiteWidth);
	
	/* Iterate through all amide groups in this secondary structure: */
	glColor(sf.amideColor);
	for(std::vector<Protein::Dipole>::const_iterator nh=sf.structure->amidesBegin;nh!=sf.structure->amidesEnd;++nh)
		{
		#if USE_GL_PICKING
		glPushName(nh->getMajorAtom()->residue->residueIndex);
		#endif
		
		Point bondSite=nh->getBondSite();
		glBegin(GL_POINTS);
		glVertex(bondSite);
		glEnd();
		glBegin(GL_LINES);
		glVertex(nh->getMajorAtom()->getPosition());
		glVertex(bondSite);
		glEnd();
		
		#if USE_GL_PICKING
		glPopName();
		#endif
		}
	
	/* Iterate through all carboxyl groups in this secondary structure: */
	glColor(sf.carboxylColor);
	for(std::vector<Protein::Dipole>::const_iterator co=sf.structure->carboxylsBegin;co!=sf.structure->carboxylsEnd;++co)
		{
		#if USE_GL_PICKING
		glPushName(co->getMajorAtom()->residue->residueIndex);
		#endif
		
		Point bondSite=co->getBondSite();
		glBegin(GL_POINTS);
		glVertex(bondSite);
		glEnd();
		glBegin(GL_LINES);
		glVertex(co->getMajorAtom()->getPosition());
		glVertex(bondSite);
		glEnd();
		
		#if USE_GL_PICKING
		glPopName();
		#endif
		}
	}

void ProteinRenderer::glDrawHydrogenCages(const ProteinRenderer::StructureFlags& sf,ProteinRenderer::DataItem& dataItem) const
	{
	/* Set point rendering parameters for this structure: */
	glPointSize(sf.hydrogenBondSiteDiameter);
	glBegin(GL_POINTS);
	
	/* Iterate through all carboxyl groups in this secondary structure: */
	glColor(sf.carboxylColor);
	for(std::vector<Protein::Dipole>::const_iterator co=sf.structure->carboxylsBegin;co!=sf.structure->carboxylsEnd;++co)
		glVertex(co->getMinorAtom()->getPosition());
	
	glEnd();
	
	/* Set line rendering parameters for this structure: */
	glLineWidth(sf.hydrogenBondSiteWidth);
	glBegin(GL_LINES);
	
	/* Iterate through all carboxyl groups in this secondary structure: */
	glColor(sf.carboxylColor);
	for(std::vector<Protein::Dipole>::const_iterator co=sf.structure->carboxylsBegin;co!=sf.structure->carboxylsEnd;++co)
		{
		glVertex(co->getMajorAtom()->getPosition());
		glVertex(co->getMinorAtom()->getPosition());
		}
	
	glEnd();
	
	/* Set hydrogen cage rendering parameters for this structure: */
	glLineWidth(sf.hydrogenCageWidth);
	glColor(sf.hydrogenCageColor);
	for(std::vector<Protein::Dipole>::const_iterator nh=sf.structure->amidesBegin;nh!=sf.structure->amidesEnd;++nh)
		{
		/* Update modelview matrix to position cage: */
		glPushMatrix();
		Scalar matrix[4][4]; // Column-major transformation matrix
		Vector y=nh->getMinorAtom()->getPosition()-nh->getMajorAtom()->getPosition();
		y.normalize();
		Vector x=nh->getAlphaCarbon()->getPosition()-nh->getMajorAtom()->getPosition();
		Vector z=Geometry::cross(x,y);
		z.normalize();
		x=Geometry::cross(y,z);
		for(int i=0;i<3;++i)
			{
			matrix[0][i]=x[i];
			matrix[1][i]=y[i];
			matrix[2][i]=z[i];
			matrix[3][i]=nh->getMinorAtom()->getPosition()[i];
			}
		matrix[0][3]=matrix[1][3]=matrix[2][3]=Scalar(0);
		matrix[3][3]=Scalar(1);
		glMultMatrix(&matrix[0][0]);
		
		/* Render the cage: */
		if(sf.hydrogenCageLarge)
			glCallList(dataItem.hydrogenCageLargeDisplayListId);
		else
			glCallList(dataItem.hydrogenCageSmallDisplayListId);
		
		glPopMatrix();
		}
	}

void ProteinRenderer::glDrawHydrogenBonds(GLContextData& contextData) const
	{
	/* Set hydrogen bond rendering parameters: */
	glLineWidth(hydrogenBondWidth);
	glColor(hydrogenBondColor);
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(4,0xAAAA);
	glBegin(GL_LINES);
	
	/* Iterate through all possible combinations of opposite dipoles: */
	for(std::vector<Protein::Dipole>::const_iterator nh=protein->amides.begin();nh!=protein->amides.end();++nh)
		for(std::vector<Protein::Dipole>::const_iterator co=protein->carboxyls.begin();co!=protein->carboxyls.end();++co)
			{
			if(formHydrogenBond(*nh,*co))
				{
				glVertex(nh->getMinorAtom()->getPosition());
				glVertex(co->getMinorAtom()->getPosition());
				}
			}
	
	glEnd();
	glDisable(GL_LINE_STIPPLE);
	}

void ProteinRenderer::glDrawCollisions(GLContextData& contextData) const
	{
	/* Set collision sphere rendering parameters: */
	glMaterial(GLMaterialEnums::FRONT,collisionSphereMaterial);
	
	/* Enable the collision sphere renderer: */
	collisionSphereRenderer.enable(modelScale,contextData);
	
	/* Calculate offsets for the 27 neighbours of a cell: */
	int cellOffsets[27];
	int* coPtr=cellOffsets;
	for(int z=-1;z<=1;++z)
		for(int y=-1;y<=1;++y)
			for(int x=-1;x<=1;++x,++coPtr)
				*coPtr=(z*numCollisionCells[1]+y)*numCollisionCells[0]+x;
	
	/* Iterate through all atoms: */
	glBegin(GL_POINTS);
	const AtomListItem* aliPtr=atomListItems;
	for(int i=0;i<protein->getNumAtoms();++i,++aliPtr)
		{
		const Protein::ChainAtom* aPtr1=aliPtr->atom;
		
		/* Check all 27 neighbouring cells: */
		for(int j=0;j<27;++j)
			{
			for(AtomListItem* cellAliPtr=aliPtr->cell[cellOffsets[j]].atoms;cellAliPtr!=0;cellAliPtr=cellAliPtr->cellSucc)
				{
				if(cellAliPtr>aliPtr)
					{
					const Protein::ChainAtom* aPtr2=cellAliPtr->atom;
					double radiusSum=(aPtr1->getCovalentRadius()+aPtr2->getCovalentRadius())*0.75;
					if(Geometry::sqrDist(aPtr1->getPosition(),aPtr2->getPosition())<Math::sqr(radiusSum))
						{
						/* Draw a big red sphere at the collision's center point: */
						double penetrationDepth=(radiusSum-Geometry::dist(aPtr1->getPosition(),aPtr2->getPosition()))*2.0;
						Point center=Geometry::mid(aPtr1->getPosition(),aPtr2->getPosition());
						glVertex4f(center[0],center[1],center[2],penetrationDepth);
						}
					}
				}
			}
		}
	glEnd();
	
	/* Disable the collision sphere renderer: */
	collisionSphereRenderer.disable(contextData);
	}

ProteinRenderer::ProteinRenderer(const Misc::ConfigurationFileSection& sConfigFileSection,const Protein* sProtein)
	:configFileSection(sConfigFileSection),protein(sProtein),
	 atomListItems(new AtomListItem[protein->getNumAtoms()]),collisionCells(0),
	 boundingBox(Geometry::Box<Scalar,3>::empty),
	 globalSettingsVersion(1),
	 drawAtoms(configFileSection.retrieveValue("./drawAtoms",false)),
	 atomMaterial(configFileSection.retrieveValue("./atomMaterial",GLMaterial(Color(0.1,0.1,0.1),Color(1.0,1.0,1.0),Color(1.0,1.0,1.0),25.0))),
	 atomTesselation(configFileSection.retrieveValue("./atomTesselation",4)),
	 mapAtomValues(false),
	 atomValues(new float[protein->getNumAtoms()]),
	 drawBonds(configFileSection.retrieveValue<bool>("./drawBonds",false)),
	 bondMaterial(configFileSection.retrieveValue<GLMaterial>("./bondMaterial",GLMaterial(Color(0.1,0.1,0.1),Color(1.0,0.0,0.0),Color(1.0,1.0,1.0),25.0))),
	 numBondVertices(configFileSection.retrieveValue<int>("./numBondVertices",24)),
	 bondRadius(configFileSection.retrieveValue<float>("./bondRadius",0.25f)),
	 alphaHelixColor(configFileSection.retrieveValue<Color>("./alphaHelixColor",Color(1.0,0.0,0.7))),
	 betaStrandColor(configFileSection.retrieveValue<Color>("./betaStrandColor",Color(0.7,0.0,1.0))),
	 coilColor(configFileSection.retrieveValue<Color>("./coilColor",Color(1.0,0.0,1.0))),
	 highlightColor(configFileSection.retrieveValue<Color>("./highlightColor",Color(1.0,1.0,0.0))),
	 drawBackbone(configFileSection.retrieveValue<bool>("./drawBackbone",true)),
	 drawBackboneRibbon(configFileSection.retrieveValue<bool>("./drawBackboneRibbon",false)),
	 backboneRibbonMaterial(configFileSection.retrieveValue<GLMaterial>("./backboneRibbonMaterial",GLMaterial(Color(0.05,0.05,0.05),Color(1.0,0.0,1.0),Color(1.0,1.0,1.0),50.0))),
	 backboneRibbonUseAllAtoms(configFileSection.retrieveValue<bool>("./backboneRibbonUseAllAtoms",false)),
	 backboneRibbonDegree(configFileSection.retrieveValue<int>("./backboneRibbonDegree",3)),
	 backboneRibbonSampleDensity(configFileSection.retrieveValue<int>("./backboneRibbonSampleDensity",16)),
	 drawCartoon(configFileSection.retrieveValue<bool>("./drawCartoon",true)),
	 cartoonMaterial(configFileSection.retrieveValue<GLMaterial>("./cartoonMaterial",GLMaterial(Color(0.05,0.05,0.05),Color(1.0,0.0,1.0),Color(1.0,1.0,1.0),50.0))),
	 cartoonSampleDensity(configFileSection.retrieveValue<int>("./cartoonSampleDensity",16)),
	 alphaHelixWidth(configFileSection.retrieveValue<Scalar>("./alphaHelixWidth",Scalar(1.5))),
	 alphaHelixThickness(configFileSection.retrieveValue<Scalar>("./alphaHelixThickness",Scalar(0.5))),
	 betaStrandWidth(configFileSection.retrieveValue<Scalar>("./betaStrandWidth",Scalar(1.5))),
	 betaStrandThickness(configFileSection.retrieveValue<Scalar>("./betaStrandThickness",Scalar(0.5))),
	 betaStrandHeadWidth(configFileSection.retrieveValue<Scalar>("./betaStrandHeadWidth",Scalar(1.333))),
	 numCoilVertices(configFileSection.retrieveValue<int>("./numCoilVertices",12)),
	 coilRadius(configFileSection.retrieveValue<Scalar>("./coilRadius",Scalar(0.333))),
	 drawHydrogenBonds(configFileSection.retrieveValue<bool>("./drawHydrogenBonds",true)),
	 hydrogenBondWidth(configFileSection.retrieveValue<float>("./hydrogenBondWidth",3.0f)),
	 hydrogenBondColor(configFileSection.retrieveValue<Color>("./hydrogenBondColor",Color(1.0,0.9,0.0))),
	 drawHydrogenBondSites(configFileSection.retrieveValue<bool>("./drawHydrogenBondSites",false)),
	 drawHydrogenCages(configFileSection.retrieveValue<bool>("./drawHydrogenCages",false)),
	 drawCollisions(configFileSection.retrieveValue<bool>("./drawCollisions",false)),
	 collisionSphereMaterial(configFileSection.retrieveValue<GLMaterial>("./collisionSphereMaterial",GLMaterial(Color(0.1,0.1,0.1),Color(1.0,0.0,0.0),Color(1.0,1.0,1.0),25.0))),
	 collisionSphereTesselation(configFileSection.retrieveValue<int>("./collisionSphereTesselation",4))
	{
	/* Initialize atom list: */
	AtomListItem* aliPtr=atomListItems;
	for(const Protein::ChainAtom* aPtr=protein->atoms;aPtr!=0;aPtr=aPtr->succ,++aliPtr)
		aliPtr->atom=aPtr;
	
	/* Initialize atom renderer: */
	atomRenderer.setVariableRadius();
	atomRenderer.setColorMaterial(true);
	
	/* Create atom color map: */
	static GLColorMap::Color mapColors[3]={GLColorMap::Color(0.0,1.0,0.0),GLColorMap::Color(1.0,1.0,0.0),GLColorMap::Color(1.0,0.0,0.0)};
	atomColorMap=GLColorMap(3,mapColors,0.0,1.0);
	
	/* Initialize atom mapping values: */
	for(int i=0;i<protein->getNumAtoms();++i)
		atomValues[i]=0.0f;
	
	/* Initialize bond renderer: */
	bondAtomRenderer.setFixedRadius(bondRadius);
	bondAtomRenderer.setColorMaterial(true);
	bondRenderer.setFixedRadius(bondRadius);
	bondRenderer.setCapped(false);
	bondRenderer.setColorMaterial(true);
	bondRenderer.setBicolor(true);
	
	/* Initialize collision sphere renderer: */
	collisionSphereRenderer.setVariableRadius();
	collisionSphereRenderer.setColorMaterial(false);
	
	/* Initialize per-structure rendering parameters: */
	updateStructureFlags();
	
	/* Initialize bond illuminator: */
	bondIlluminator.setMaterial(bondMaterial);
	bondIlluminator.setSceneCenter(GLVector<GLfloat,3>(protein->calcCentroid().getComponents()));
	bondIlluminator.enableAutoView();
	bondIlluminator.enableAutoLight(GL_LIGHT0);
	}

ProteinRenderer::~ProteinRenderer(void)
	{
	delete[] atomListItems;
	delete[] collisionCells;
	delete[] atomValues;
	}

void ProteinRenderer::initContext(GLContextData& contextData) const
	{
	/* Create a new context entry: */
	DataItem* dataItem=new DataItem();
	contextData.addDataItem(this,dataItem);
	
	/* Create hydrogen cage template: */
	Scalar a=Protein::Dipole::NHdist;
	Scalar b=Protein::Dipole::hydrogenBondRangeMinor.getMax();
	Scalar c=Protein::Dipole::hydrogenBondRangeMajor.getMin();
	const int nb=6;
	Scalar x[nb],y[nb];
	x[0]=Scalar(0);
	y[0]=Protein::Dipole::hydrogenBondRangeMinor.getMin();
	Scalar tan2AngleMax=(Scalar(1)-Math::sqr(Protein::Dipole::cosAngleMin))/Math::sqr(Protein::Dipole::cosAngleMin);
	Scalar denom=Scalar(1)+tan2AngleMax;
	y[1]=(-a+Math::sqrt(Math::sqr(c)*denom-Math::sqr(a)*tan2AngleMax))/denom;
	x[1]=y[1]*Math::sqrt(tan2AngleMax);
	x[nb-1]=Scalar(0);
	y[nb-1]=b;
	Scalar angleMax=Math::acos(Protein::Dipole::cosAngleMin);
	for(int i=2;i<nb-1;++i)
		{
		Scalar angle=Scalar(nb-i)*angleMax/Scalar(nb-2);
		x[i]=b*Math::sin(angle);
		y[i]=b*Math::cos(angle);
		}
	const int na=8;
	Scalar sj[na],cj[na];
	for(int j=0;j<na;++j)
		{
		Scalar angle=Scalar(j)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(na);
		sj[j]=Math::sin(angle);
		cj[j]=Math::cos(angle);
		}

	/* Create a display list for small hydrogen cages: */
	glNewList(dataItem->hydrogenCageSmallDisplayListId,GL_COMPILE);
	for(int i=2;i<nb-1;++i)
		{
		glBegin(GL_LINE_LOOP);
		for(int j=0;j<na;++j)
			glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
		glEnd();
		}
	for(int j=0;j<na/2;++j)
		{
		glBegin(GL_LINE_STRIP);
		for(int i=2;i<nb;++i)
			glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
		for(int i=nb-2;i>=2;--i)
			glVertex(x[i]*cj[j+na/2],y[i],x[i]*sj[j+na/2]);
		glEnd();
		}
	glEndList();
	
	/* Create a display list for large hydrogen cages: */
	glNewList(dataItem->hydrogenCageLargeDisplayListId,GL_COMPILE);
	for(int i=1;i<nb-1;++i)
		{
		glBegin(GL_LINE_LOOP);
		for(int j=0;j<na;++j)
			glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
		glEnd();
		}
	for(int j=0;j<na/2;++j)
		{
		glBegin(GL_LINE_LOOP);
		for(int i=1;i<nb;++i)
			glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
		for(int i=nb-2;i>=0;--i)
			glVertex(x[i]*cj[j+na/2],y[i],x[i]*sj[j+na/2]);
		glEnd();
		}
	glEndList();
	}

void ProteinRenderer::updateStructureFlags(void)
	{
	/* Synchronize structure parameter list with protein's secondary structure sequence: */
	structureFlags.clear();
	int startResidue=0;
	for(const Protein::SecondaryStructure* sPtr=protein->secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
		{
		Misc::ConfigurationFileSection structureSection=configFileSection;
		switch(sPtr->getStructureType())
			{
			case Protein::SecondaryStructure::COIL:
				structureSection.setSection("./Coil");
				break;
			
			case Protein::SecondaryStructure::ALPHA_HELIX:
				structureSection.setSection("./AlphaHelix");
				break;
			
			case Protein::SecondaryStructure::BETA_STRAND:
				structureSection.setSection("./BetaStrand");
				break;
			
			default:
				; // Just to make compiler happy
			}
		int endResidue=startResidue+sPtr->getNumResidues();
		structureFlags.push_back(StructureFlags(structureSection,sPtr,startResidue,endResidue));
		startResidue=endResidue;
		}
	resetAllBackboneColor();
	
	/* Create backbone ribbon spline: */
	createBackboneRibbonSpline();
	
	/* Create cartoon splines: */
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		createCartoonSpline(*sfIt);
	
	/* Fill in protein data structures: */
	updateProtein();
	}

void ProteinRenderer::updateProtein(void)
	{
	/* Update the protein's bounding box: */
	Point min,max;
	min=max=protein->atoms->getPosition();
	for(const Protein::ChainAtom* aPtr=protein->atoms->succ;aPtr!=0;aPtr=aPtr->succ)
		{
		for(int i=0;i<3;++i)
			{
			if(min[i]>aPtr->getPosition()[i])
				min[i]=aPtr->getPosition()[i];
			else if(max[i]<aPtr->getPosition()[i])
				max[i]=aPtr->getPosition()[i];
			}
		}
	for(int i=0;i<3;++i)
		{
		min[i]-=maxAtomRadius*2.0;
		max[i]+=maxAtomRadius*2.0;
		}
	boundingBox=Geometry::Box<Scalar,3>(min,max);
	
	if(drawCartoon)
		{
		/* Create and evaluate cartoon splines: */
		for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
			updateCartoonSpline(*sfIt);
		}
	
	if(drawCollisions)
		{
		/* Create the collision detection data structure: */
		double optCollisionCellSize=Math::pow(boundingBox.getSize(0)*boundingBox.getSize(1)*boundingBox.getSize(2)/double(protein->getNumAtoms()),1.0/3.0);
		double maxCollisionRadius=0.0;
		for(int i=0;i<118;++i)
			{
			double collisionRadius=Atom::getCovalentRadius((Atom::Element)i)*0.75;
			if(maxCollisionRadius<collisionRadius)
				maxCollisionRadius=collisionRadius;
			}
		if(optCollisionCellSize<2.0*maxCollisionRadius)
			optCollisionCellSize=2.0*maxCollisionRadius;
		for(int i=0;i<3;++i)
			{
			numCollisionCells[i]=int(Math::floor(boundingBox.getSize(i)/optCollisionCellSize));
			collisionCellSize[i]=boundingBox.getSize(i)/double(numCollisionCells[i]);
			numCollisionCells[i]+=2;
			}
		delete[] collisionCells;
		collisionCells=new CollisionCell[numCollisionCells[2]*numCollisionCells[1]*numCollisionCells[0]];
		collisionCellBase=&collisionCells[(1*numCollisionCells[1]+1)*numCollisionCells[0]+1];

		/* Insert all atoms into the collision detection structure: */
		AtomListItem* aliPtr=atomListItems;
		for(int i=0;i<protein->getNumAtoms();++i,++aliPtr)
			{
			int cellIndex[3];
			for(int i=0;i<3;++i)
				cellIndex[i]=int(Math::floor((aliPtr->atom->getPosition()[i]-boundingBox.getOrigin()[i])/collisionCellSize[i]));
			aliPtr->cell=&collisionCellBase[(cellIndex[2]*numCollisionCells[1]+cellIndex[1])*numCollisionCells[0]+cellIndex[0]];
			aliPtr->cellSucc=aliPtr->cell->atoms;
			aliPtr->cell->atoms=aliPtr;
			}
		}
	
	/* Create the backbone ribbon spline: */
	createBackboneRibbonSpline();
	}

void ProteinRenderer::setModelScale(GLfloat newModelScale)
	{
	modelScale=newModelScale;
	}

void ProteinRenderer::setViewDirection(const GLLineIlluminator::Vector& viewDirection)
	{
	bondIlluminator.disableAutoView();
	bondIlluminator.setViewDirection(viewDirection);
	}

void ProteinRenderer::setLightDirection(const GLLineIlluminator::Vector& lightDirection)
	{
	bondIlluminator.disableAutoLight();
	bondIlluminator.setLightDirection(lightDirection);
	}

void ProteinRenderer::glRenderAction(GLContextData& contextData) const
	{
	/* Get a pointer to the context entry: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	/* Check if the secondary structure list needs to be updated: */
	if(dataItem->secondaryStructureListVersion!=protein->getSecondaryStructureVersion())
		{
		/* Delete old secondary structure state: */
		if(dataItem->numSecondaryStructures>0)
			{
			glDeleteLists(dataItem->secondaryStructureListBaseId,dataItem->numSecondaryStructures);
			delete[] dataItem->secondaryStructureSettingsVersions;
			delete[] dataItem->secondaryStructureVersions;
			}
		
		/* Initialize new secondary structure state: */
		dataItem->numSecondaryStructures=protein->getNumStructures();
		dataItem->secondaryStructureListBaseId=glGenLists(dataItem->numSecondaryStructures);
		dataItem->secondaryStructureSettingsVersions=new unsigned int[dataItem->numSecondaryStructures];
		dataItem->secondaryStructureVersions=new unsigned int[dataItem->numSecondaryStructures];
		for(int i=0;i<dataItem->numSecondaryStructures;++i)
			{
			dataItem->secondaryStructureSettingsVersions[i]=0;
			dataItem->secondaryStructureVersions[i]=0;
			}
		
		/* Update the secondary structure list version number: */
		dataItem->secondaryStructureListVersion=protein->getSecondaryStructureVersion();
		}
	
	/* Invalidate all secondary structures if global rendering settings changed: */
	if(dataItem->globalSettingsVersion!=globalSettingsVersion)
		{
		for(int i=0;i<dataItem->numSecondaryStructures;++i)
			dataItem->secondaryStructureSettingsVersions[i]=0;
		dataItem->globalSettingsVersion=globalSettingsVersion;
		}
	
	/* Check if atom's van-der-Waals spheres need to be drawn: */
	if(drawAtoms)
		{
		/* Set atom rendering parameters: */
		glEnable(GL_LIGHTING);
		glEnable(GL_NORMALIZE);
		glMaterial(GLMaterialEnums::FRONT,atomMaterial);
		
		/* Enable the atom renderer: */
		atomRenderer.enable(modelScale,contextData);
		
		/* Render each secondary structure independently: */
		#if !USE_GL_PICKING
		glBegin(GL_POINTS);
		#endif
		GLuint structureIndex=0;
		for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt,++structureIndex)
			if(sfIt->drawAtoms)
				{
				#if USE_GL_PICKING
				/* Push the structure's index onto the name stack: */
				glPushName(structureIndex);
				#endif
				
				/* Iterate through all residues in this secondary structure: */
				for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
					{
					#if USE_GL_PICKING
					/* Push residue (PDB) index onto the selection name stack: */
					glPushName(rPtr->residueIndex);
					#endif
					
					/* Iterate through all atoms in this residue: */
					#if USE_GL_PICKING
					glBegin(GL_POINTS);
					#endif
					for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
						{
						if(mapAtomValues)
							glColor(atomColorMap(atomValues[aPtr->getAtomIndex()]));
						else
							glColor(elementColors[aPtr->getType()]);
						glVertex4f(aPtr->getPosition()[0],aPtr->getPosition()[1],aPtr->getPosition()[2],aPtr->getVanDerWaalsRadius());
						}
					#if USE_GL_PICKING
					glEnd();
					#endif
					
					#if USE_GL_PICKING
					glPopName();
					#endif
					}
				
				#if USE_GL_PICKING
				glPopName();
				#endif
				}
		#if !USE_GL_PICKING
		glEnd();
		#endif
		
		/* Disable the atom renderer: */
		atomRenderer.disable(contextData);
		}
	
	/* Check if side chain bond cylinders need to be drawn: */
	if(drawBonds)
		{
		/* Set bond cylinder rendering parameters: */
		glMaterial(GLMaterialEnums::FRONT,bondMaterial);
		
		/* Enable the bond atom renderer: */
		bondAtomRenderer.enable(modelScale,contextData);
		
		/* Render each secondary structure independently: */
		#if !USE_GL_PICKING
		glBegin(GL_POINTS);
		#endif
		GLuint structureIndex=0;
		for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt,++structureIndex)
			if(sfIt->drawBonds)
				{
				#if USE_GL_PICKING
				/* Push the structure's index onto the name stack: */
				glPushName(structureIndex);
				#endif
				
				/* Iterate through all residues in this secondary structure: */
				for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
					{
					#if USE_GL_PICKING
					/* Push residue (PDB) index onto the selection name stack: */
					glPushName(rPtr->residueIndex);
					#endif
					
					/* Iterate through all atoms in this residue: */
					#if USE_GL_PICKING
					glBegin(GL_POINTS);
					#endif
					for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
						{
						glColor(elementColors[aPtr->getType()]);
						glVertex(aPtr->getPosition());
						}
					#if USE_GL_PICKING
					glEnd();
					#endif
					
					#if USE_GL_PICKING
					glPopName();
					#endif
					}
				
				#if USE_GL_PICKING
				glPopName();
				#endif
				}
		#if !USE_GL_PICKING
		glEnd();
		#endif
		
		/* Disable the bond atom renderer: */
		bondAtomRenderer.disable(contextData);
		
		/* Enable the bond cylinder renderer: */
		bondRenderer.enable(modelScale,contextData);
		
		/* Render each secondary structure independently: */
		#if !USE_GL_PICKING
		glBegin(GL_LINES);
		#endif
		structureIndex=0;
		for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt,++structureIndex)
			if(sfIt->drawBonds)
				{
				#if USE_GL_PICKING
				/* Push the structure's index onto the name stack: */
				glPushName(structureIndex);
				#endif
				
				/* Iterate through all residues in this secondary structure: */
				for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
					{
					#if USE_GL_PICKING
					/* Push residue (PDB) index onto the selection name stack: */
					glPushName(rPtr->residueIndex);
					#endif
					
					/* Iterate through all atoms in this residue: */
					#if USE_GL_PICKING
					glBegin(GL_LINES);
					#endif
					for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
						{
						/* Iterate through all bonds for this atom: */
						for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
							if(static_cast<const Atom*>(aPtr)<*bIt)
								{
								glColor(elementColors[aPtr->getType()]);
								glVertex(aPtr->getPosition());
								glColor(elementColors[(*bIt)->getType()]);
								glVertex((*bIt)->getPosition());
								}
						}
					#if USE_GL_PICKING
					glEnd();
					#endif
					
					#if USE_GL_PICKING
					glPopName();
					#endif
					}
				
				#if USE_GL_PICKING
				glPopName();
				#endif
				}
		#if !USE_GL_PICKING
		glEnd();
		#endif
		
		/* Disable the bond cylinder renderer: */
		bondRenderer.disable(contextData);
		}
	
	/* Render each secondary structure independently: */
	GLuint structureIndex=0;
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt,++structureIndex)
		{
		#if USE_GL_PICKING
		/* Push the structure's index onto the name stack: */
		glPushName(structureIndex);
		#endif
		
		/* Check if the secondary structure needs to be re-rendered: */
		if(dataItem->secondaryStructureSettingsVersions[structureIndex]!=sfIt->settingsVersion||dataItem->secondaryStructureVersions[structureIndex]!=sfIt->structure->getVersion())
			{
			/* Render and cache the secondary structure: */
			glNewList(dataItem->secondaryStructureListBaseId+structureIndex,GL_COMPILE_AND_EXECUTE);
			
			/* Render line parts of secondary structure visualization: */
			glDisable(GL_LIGHTING);
			if(drawBackbone&&sfIt->drawBackbone)
				glDrawBackbone(*sfIt,*dataItem);
			#if 0
			if(drawBonds&&sfIt->drawBonds)
				{
				/* Enable illuminated lines: */
				bondIlluminator.enableLighting(contextData);
				
				glDrawBonds(*sfIt,*dataItem);
				
				/* Disable illuminated lines: */
				bondIlluminator.disableLighting(contextData);
				}
			#endif
			if(drawHydrogenBondSites&&sfIt->drawHydrogenBondSites)
				glDrawHydrogenBondSites(*sfIt,*dataItem);
			if(drawHydrogenCages&&sfIt->drawHydrogenCages)
				glDrawHydrogenCages(*sfIt,*dataItem);
			
			/* Render polygon parts of secondary structure visualization: */
			glEnable(GL_LIGHTING);
			glEnable(GL_NORMALIZE);
			if(drawBackboneRibbon&&sfIt->drawBackboneRibbon)
				glDrawBackboneRibbon(*sfIt,*dataItem);
			if(drawCartoon&&sfIt->drawCartoon)
				glDrawCartoon(*sfIt,*dataItem);
			
			glEndList();
			
			/* Update secondary structure's version number: */
			dataItem->secondaryStructureSettingsVersions[structureIndex]=sfIt->settingsVersion;
			dataItem->secondaryStructureVersions[structureIndex]=sfIt->structure->getVersion();
			}
		else
			{
			/* Call the secondary structure's cache display list: */
			glCallList(dataItem->secondaryStructureListBaseId+structureIndex);
			}
		
		#if USE_GL_PICKING
		glPopName();
		#endif
		}
	
	/* Render line parts of structure-independent parts of protein visualization: */
	glDisable(GL_LIGHTING);
	if(drawHydrogenBonds)
		glDrawHydrogenBonds(contextData);
	
	/* Render polygon parts of structure-independent parts of protein visualization: */
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);
	if(drawCollisions)
		glDrawCollisions(contextData);
	}

void ProteinRenderer::highlightResidue(GLContextData& contextData,const Protein::Residue* rPtr) const
	{
	/* Set bond rendering parameters: */
	glMaterial(GLMaterialEnums::FRONT,bondMaterial);
	
	/* Enable the bond renderer: */
	bondRenderer.enable(modelScale,contextData);
	
	/* Iterate through all atoms in this residue: */
	glBegin(GL_LINES);
	for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
		{
		/* Iterate through all bonds for this atom: */
		for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
			if(static_cast<const Atom*>(aPtr)<*bIt)
				{
				glColor(elementColors[aPtr->getType()]);
				glVertex(aPtr->getPosition());
				glColor(elementColors[(*bIt)->getType()]);
				glVertex((*bIt)->getPosition());
				}
		}
	glEnd();
	
	/* Disable the bond renderer: */
	bondRenderer.disable(contextData);
	}

int ProteinRenderer::glPick(GLContextData& contextData,int which) const
	{
	/* Start OpenGL selection mode: */
	GLsizei hitBufferSize=structureFlags.size()*7*5; // Enough entries for every structure hitting on every rendering pass with 2 names on the stack
	GLuint* hitBuffer=new GLuint[hitBufferSize];
	glSelectBuffer(hitBufferSize,hitBuffer);
	glRenderMode(GL_SELECT);
	glInitNames();
	
	/* Render protein as usual: */
	glRenderAction(contextData);
	
	/* Process hit records: */
	int numHits=glRenderMode(GL_RENDER);
	GLuint* hPtr=hitBuffer;
	GLuint minZ=0xffffffffU;
	int minStructure=-1;
	for(int i=0;i<numHits;++i)
		{
		GLuint numNames=hPtr[0];
		GLuint z1=hPtr[1];
		hPtr+=3;
		if(int(numNames)>which&&z1<minZ)
			{
			minStructure=int(hPtr[which]);
			minZ=z1;
			}
		hPtr+=numNames;
		}
	delete[] hitBuffer;
	
	return minStructure;
	}

ProteinRenderer::Color ProteinRenderer::getResidueBackboneColor(const Protein::Residue* rPtr) const
	{
	/* Find structure flags for the residue's secondary structure: */
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==rPtr->secondaryStructure)
			return sfIt->backboneColor;
		}
	
	return Color(0,0,0);
	}

void ProteinRenderer::lockResidueRange(int firstResidueIndex,int numResidues,bool locked)
	{
	/* Check each structure whether it is included in the residue range: */
	int structureFirstResidueIndex=0;
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		int structureNumResidues=sfIt->structure->getNumResidues();
		if(structureFirstResidueIndex>=firstResidueIndex&&structureFirstResidueIndex+structureNumResidues<=firstResidueIndex+numResidues)
			{
			if(locked)
				{
				/* Set the structure's backbone color to red: */
				sfIt->backboneColor=Color(1.0f,0.0f,0.0f);
				++sfIt->settingsVersion;
				}
			else
				{
				/* Reset the structure's backbone color: */
				Color newColor;
				switch(sfIt->structure->structureType)
					{
					case Protein::SecondaryStructure::COIL:
						newColor=coilColor;
						break;
					
					case Protein::SecondaryStructure::ALPHA_HELIX:
						newColor=alphaHelixColor;
						break;
					
					case Protein::SecondaryStructure::BETA_STRAND:
						newColor=betaStrandColor;
						break;
					
					default:
						; // Just to make compiler happy
					}
				if(sfIt->backboneColor!=newColor)
					{
					sfIt->backboneColor=newColor;
					++sfIt->settingsVersion;
					}
				}
			}
		structureFirstResidueIndex+=structureNumResidues;
		}
	}

void ProteinRenderer::setDrawAtoms(bool newDrawAtoms)
	{
	drawAtoms=newDrawAtoms;
	++globalSettingsVersion;
	}

void ProteinRenderer::setMapAtomValues(bool newMapAtomValues)
	{
	mapAtomValues=newMapAtomValues;
	++globalSettingsVersion;
	}

void ProteinRenderer::setDrawBonds(bool newDrawBonds)
	{
	drawBonds=newDrawBonds;
	++globalSettingsVersion;
	}

void ProteinRenderer::setDrawBackbone(bool newDrawBackbone)
	{
	drawBackbone=newDrawBackbone;
	++globalSettingsVersion;
	}

void ProteinRenderer::setDrawBackboneRibbon(bool newDrawBackboneRibbon)
	{
	drawBackboneRibbon=newDrawBackboneRibbon;
	++globalSettingsVersion;
	}

void ProteinRenderer::setDrawCartoon(bool newDrawCartoon)
	{
	drawCartoon=newDrawCartoon;
	++globalSettingsVersion;
	if(drawCartoon)
		{
		/* Create and evaluate cartoon splines: */
		for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
			updateCartoonSpline(*sfIt);
		}
	}

void ProteinRenderer::setDrawHydrogenBonds(bool newDrawHydrogenBonds)
	{
	drawHydrogenBonds=newDrawHydrogenBonds;
	}

void ProteinRenderer::setDrawHydrogenBondSites(bool newDrawHydrogenBondSites)
	{
	drawHydrogenBondSites=newDrawHydrogenBondSites;
	++globalSettingsVersion;
	}

void ProteinRenderer::setDrawHydrogenCages(bool newDrawHydrogenCages)
	{
	drawHydrogenCages=newDrawHydrogenCages;
	++globalSettingsVersion;
	}

void ProteinRenderer::setDrawCollisions(bool newDrawCollisions)
	{
	drawCollisions=newDrawCollisions;
	++globalSettingsVersion;
	if(drawCollisions)
		updateProtein(); // Might have to re-calculate collision data structure
	}

bool ProteinRenderer::getDrawAtoms(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==selector.structure)
			return sfIt->drawAtoms;
		}
	return false;
	}

void ProteinRenderer::setDrawAtoms(const Protein::StructureSelector& selector,bool newDrawAtoms)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->drawAtoms=newDrawAtoms;
			++sfIt->settingsVersion;
			}
		}
	}

bool ProteinRenderer::getDrawBonds(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==selector.structure)
			return sfIt->drawBonds;
		}
	return false;
	}

void ProteinRenderer::setDrawBonds(const Protein::StructureSelector& selector,bool newDrawBonds)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->drawBonds=newDrawBonds;
			++sfIt->settingsVersion;
			}
		}
	}

bool ProteinRenderer::getDrawBackbone(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==selector.structure)
			return sfIt->drawBackbone;
		}
	return false;
	}

void ProteinRenderer::setDrawBackbone(const Protein::StructureSelector& selector,bool newDrawBackbone)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->drawBackbone=newDrawBackbone;
			++sfIt->settingsVersion;
			}
		}
	}

bool ProteinRenderer::getDrawBackboneRibbon(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==selector.structure)
			return sfIt->drawBackboneRibbon;
		}
	return false;
	}

void ProteinRenderer::setDrawBackboneRibbon(const Protein::StructureSelector& selector,bool newDrawBackboneRibbon)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->drawBackboneRibbon=newDrawBackboneRibbon;
			++sfIt->settingsVersion;
			}
		}
	}

bool ProteinRenderer::getDrawCartoon(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==selector.structure)
			return sfIt->drawCartoon;
		}
	return false;
	}

void ProteinRenderer::setDrawCartoon(const Protein::StructureSelector& selector,bool newDrawCartoon)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->drawCartoon=newDrawCartoon;
			++sfIt->settingsVersion;
			if(drawCartoon&&sfIt->drawCartoon)
				updateCartoonSpline(*sfIt);
			}
		}
	}

void ProteinRenderer::resetAllBackboneColor(void)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		Color newColor;
		switch(sfIt->structure->structureType)
			{
			case Protein::SecondaryStructure::COIL:
				newColor=coilColor;
				break;
			
			case Protein::SecondaryStructure::ALPHA_HELIX:
				newColor=alphaHelixColor;
				break;
			
			case Protein::SecondaryStructure::BETA_STRAND:
				newColor=betaStrandColor;
				break;
			
			default:
				newColor=Color(0.5f,0.5f,0.5f);
			}
		if(sfIt->backboneColor!=newColor)
			{
			sfIt->backboneColor=newColor;
			++sfIt->settingsVersion;
			}
		}
	}

void ProteinRenderer::resetBackboneColor(const Protein::StructureSelector& selector)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			Color newColor;
			switch(sfIt->structure->structureType)
				{
				case Protein::SecondaryStructure::COIL:
					newColor=coilColor;
					break;
				
				case Protein::SecondaryStructure::ALPHA_HELIX:
					newColor=alphaHelixColor;
					break;
				
				case Protein::SecondaryStructure::BETA_STRAND:
					newColor=betaStrandColor;
					break;
				
				default:
					newColor=Color(0.5f,0.5f,0.5f);
				}
			if(sfIt->backboneColor!=newColor)
				{
				sfIt->backboneColor=newColor;
				++sfIt->settingsVersion;
				}
			}
		}
	}

void ProteinRenderer::setBackboneColor(const Protein::StructureSelector& selector,const Color& newBackboneColor)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->backboneColor=newBackboneColor;
			++sfIt->settingsVersion;
			}
		}
	}

bool ProteinRenderer::getDrawHydrogenBondSites(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			return sfIt->drawHydrogenBondSites;
		}
	return false;
	}

void ProteinRenderer::setAllDrawHydrogenBondSites(bool newDrawHydrogenBondSites)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->drawHydrogenBondSites!=newDrawHydrogenBondSites)
			{
			sfIt->drawHydrogenBondSites=newDrawHydrogenBondSites;
			++sfIt->settingsVersion;
			}
		}
	}

void ProteinRenderer::setDrawHydrogenBondSites(const Protein::StructureSelector& selector,bool newDrawHydrogenBondSites)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->drawHydrogenBondSites=newDrawHydrogenBondSites;
			++sfIt->settingsVersion;
			}
		}
	}

bool ProteinRenderer::getDrawHydrogenCages(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==selector.structure)
			return sfIt->drawHydrogenCages;
		}
	return false;
	}

void ProteinRenderer::setDrawHydrogenCages(const Protein::StructureSelector& selector,bool newDrawHydrogenCages)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->drawHydrogenCages=newDrawHydrogenCages;
			++sfIt->settingsVersion;
			}
		}
	}

bool ProteinRenderer::getDrawLargeHydrogenCages(const Protein::StructureSelector& selector) const
	{
	for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(sfIt->structure==selector.structure)
			return sfIt->hydrogenCageLarge;
		}
	return false;
	}

void ProteinRenderer::setDrawLargeHydrogenCages(const Protein::StructureSelector& selector,bool newDrawLargeHydrogenCages)
	{
	for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
		{
		if(selector.structure==sfIt->structure)
			{
			sfIt->hydrogenCageLarge=newDrawLargeHydrogenCages;
			++sfIt->settingsVersion;
			}
		}
	}

}
