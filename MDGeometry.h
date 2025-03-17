/***********************************************************************
MDGeometry - Basic definitions for geometry types in molecular dynamics
and kinematics.
Copyright (c) 2002-2020 Oliver Kreylos

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

#ifndef MDGEOMETRY_INCLUDED
#define MDGEOMETRY_INCLUDED

#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/Ray.h>

namespace MD {

typedef double Scalar; // Basic scalar type for geometry representations
typedef Geometry::Point<Scalar,3> Position; // Type for atom positions
typedef Geometry::Point<Scalar,3> Point; // Type for generic points
typedef Geometry::Vector<Scalar,3> Vector; // Type for generic vectors
typedef Geometry::Ray<Scalar,3> Ray; // Type for rays (used for picking and such)
}

#endif
