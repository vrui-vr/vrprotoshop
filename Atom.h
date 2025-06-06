/***********************************************************************
Atom - Class to represent single atoms, not necessarily as parts of
proteins.
Copyright (c) 2001-2020 Oliver Kreylos

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

#ifndef ATOM_INCLUDED
#define ATOM_INCLUDED

#include <vector>
#include <Geometry/Sphere.h>

#include "MDGeometry.h"

namespace MD {

class Atom
	{
	/* Embedded classes: */
	public:
	enum Element // Enumerated type for elements (standard names)
		{
		H=0,He,
		Li,Be,B,C,N,O,F,Ne,
		Na,Mg,Al,Si,P,S,Cl,Ar,
		K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,
		Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,
		Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hq,Tl,Pb,Bi,Po,At,Rn,
		Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Uun,Uuu,Uub,UUt,UUq,Uup,Uuh,Uus,Uuo
		};
	
	/* Static elements: */
	private:
	static const char elementNames[118][4]; // Table of standard element names
	static const Scalar elementRadii[118]; // Table of atomic radii (in Angstrom)
	static const Scalar elementCovalentRadii[118]; // Table of covalent atomic radii (in Angstrom)
	static const Scalar elementVanDerWaalsRadii[118]; // Table of van-der-Waals atomic radii (in Angstrom)
	
	/* Elements: */
	protected:
	Element type; // The atom's type
	std::vector<Atom*> bonds; // List of pointers to other atoms representing covalent bonds
	Position position; // The atom's position (in Angstrom)
	
	/* Constructors and destructors: */
	public:
	Atom(Element sType,const Position& sPosition)
		:type(sType),position(sPosition)
		{
		};
	
	/* Methods: */
	static const char* getElementName(Element element) // Returns the name of an element
		{
		return elementNames[element];
		};
	static Scalar getRadius(Element element) // Returns radius of an element
		{
		return elementRadii[element];
		};
	static Scalar getCovalentRadius(Element element) // Returns radius of an element
		{
		return elementCovalentRadii[element];
		};
	static Scalar getVanDerWaalsRadius(Element element) // Returns radius of an element
		{
		return elementVanDerWaalsRadii[element];
		};
	static Element parseType(const char* elementName); // Converts an element name into an element type
	Element getType(void) const // Returns an atom's element type
		{
		return type;
		};
	const char* getElementName(void) const // Returns the element type of an atom as a string
		{
		return elementNames[type];
		};
	Scalar getRadius(void) const // Returns the radius of an atom
		{
		return elementRadii[type];
		};
	Scalar getCovalentRadius(void) const // Returns the covalent radius of an atom
		{
		return elementCovalentRadii[type];
		};
	Scalar getVanDerWaalsRadius(void) const // Returns the van-der-Waals radius of an atom
		{
		return elementVanDerWaalsRadii[type];
		};
	Position getPosition(void) const // Returns the position of an atom
		{
		return position;
		};
	void setPosition(const Position& newPosition) // Sets the position of an atom
		{
		position=newPosition;
		};
	const std::vector<Atom*>& getBonds(void) const // Returns the list of covalent bonds
		{
		return bonds;
		};
	std::vector<Atom*>& getBonds(void) // Ditto
		{
		return bonds;
		};
	bool isBonded(void) const // Returns true if the atom has at least one covalent bond
		{
		return !bonds.empty();
		};
	friend double dist(const Atom& atom1,const Atom& atom2) // Returns the Euklidean distance between two atoms
		{
		return Geometry::dist(atom1.position,atom2.position);
		};
	friend void bond(Atom& atom1,Atom& atom2); // Creates a covalent bond between two atoms
	Scalar intersectRay(const Ray& ray) const // Intersects atom with ray
		{
		return Geometry::Sphere<Scalar,3>(position,elementVanDerWaalsRadii[type]).intersectRay(ray).getParameter();
		};
	};

}

#endif
