/***********************************************************************
CatmullRomSpline - Class to represent and evaluate Catmull-Rom splines.
Copyright (c) 2020 Oliver Kreylos
***********************************************************************/

#ifndef CATMULLROMSPLINE_INCLUDED
#define CATMULLROMSPLINE_INCLUDED

#include <stddef.h>
#include <Math/Math.h>
#include <Geometry/Point.h>
#include <Geometry/Vector.h>

template <class ScalarParam,int dimensionParam>
class CatmullRomSpline
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // Scalar type
	static const int dimension=dimensionParam; // Affine space dimension
	typedef Geometry::Point<ScalarParam,dimensionParam> Point; // Type for affine points
	typedef Geometry::Vector<ScalarParam,dimensionParam> Vector; // Type for vectors
	
	struct ControlPoint // Structure representing a Catmull-Rom control point
		{
		/* Elements: */
		public:
		Point position; // Control point position
		Vector tangent; // Control point tangent vector
		
		/* Constructors and destructors: */
		ControlPoint(void) // Dummy constructor
			{
			}
		ControlPoint(const Point& sPosition,const Vector& sTangent) // Elementwise constructor
			:position(sPosition),tangent(sTangent)
			{
			}
		};
	
	/* Elements: */
	private:
	size_t numControlPoints; // Number of control points in the spline
	ControlPoint* controlPoints; // Array of control points
	
	/* Constructors and destructors: */
	public:
	CatmullRomSpline(void) // Creates an invalid Catmull-Rom spline
		:numControlPoints(0),controlPoints(0)
		{
		}
	CatmullRomSpline(size_t sNumControlPoints) // Creates a Catmull-Rom spline with the given number of control points
		:numControlPoints(sNumControlPoints),
		 controlPoints(numControlPoints!=0?new ControlPoint[numControlPoints]:0)
		{
		}
	CatmullRomSpline(size_t sNumControlPoints,ControlPoint sControlPoints[]) // Creates a Catmull-Rom spline with the given control point array
		:numControlPoints(sNumControlPoints),
		 controlPoints(numControlPoints!=0?new ControlPoint[numControlPoints]:0)
		{
		/* Copy the given control points: */
		for(size_t i=0;i<numControlPoints;++i)
			controlPoints[i]=sControlPoints[i];
		}
	CatmullRomSpline(const CatmullRomSpline& source) // Copy constructor
		:numControlPoints(source.numControlPoints),
		 controlPoints(numControlPoints!=0?new ControlPoint[numControlPoints]:0)
		{
		/* Copy the source's control points: */
		for(size_t i=0;i<numControlPoints;++i)
			controlPoints[i]=source.controlPoints[i];
		}
	~CatmullRomSpline(void) // Destroys the Catmull-Rom spline
		{
		delete[] controlPoints;
		}
	
	/* Methods: */
	CatmullRomSpline& operator=(const CatmullRomSpline& source) // Assignment operator
		{
		/* Check if the source has a different number of control points: */
		if(numControlPoints!=source.numControlPoints)
			{
			/* Re-allocate the control point array: */
			delete[] controlPoints;
			numControlPoints=source.numControlPoints;
			controlPoints=numControlPoints!=0?new ControlPoint[numControlPoints]:0;
			}
		
		/* Copy the source's control points: */
		for(size_t i=0;i<numControlPoints;++i)
			controlPoints[i]=source.controlPoints[i];
		
		return *this;
		}
	size_t getNumControlPoints(void) const // Returns the number of control points
		{
		return numControlPoints;
		}
	void setNumControlPoints(size_t newNumControlPoints) // Sets the spline's number of control points
		{
		if(numControlPoints!=newNumControlPoints)
			{
			/* Re-allocate the control point array: */
			delete[] controlPoints;
			numControlPoints=newNumControlPoints;
			controlPoints=numControlPoints!=0?new ControlPoint[numControlPoints]:0;
			}
		}
	const ControlPoint& operator[](size_t index) const // Returns a control point
		{
		return controlPoints[index];
		}
	ControlPoint& operator[](size_t index) // Ditto
		{
		return controlPoints[index];
		}
	size_t getNumSegments(void) const // Returns the number of spline segments
		{
		return numControlPoints-1;
		}
	Point operator()(Scalar parameter) const // Evaluates the Catmull-Rom spline for the given spline parameter in [0, numControlPoints-1]
		{
		/* Find the spline segment containing the given parameter: */
		size_t segment;
		if(parameter<=Scalar(0))
			segment=0;
		else if(parameter>=Scalar(numControlPoints-1))
			segment=numControlPoints-2;
		else
			segment=size_t(Math::floor(parameter));
		parameter-=Scalar(segment);
		
		/* Evaluate the spline segment for the fractional parameter: */
		const ControlPoint& cp0=controlPoints[segment];
		Point c1=cp0.position+cp0.tangent;
		const ControlPoint& cp1=controlPoints[segment+1];
		Point c2=cp1.position-cp1.tangent;
		Point m0=Geometry::affineCombination(cp0.position,c1,parameter);
		Point m1=Geometry::affineCombination(c1,c2,parameter);
		Point m2=Geometry::affineCombination(c2,cp1.position,parameter);
		Point m3=Geometry::affineCombination(m0,m1,parameter);
		Point m4=Geometry::affineCombination(m1,m2,parameter);
		return Geometry::affineCombination(m3,m4,parameter);
		}
	Vector diff(Scalar parameter) const // Evaluates the Catmull-Rom spline's derivative for the given spline parameter in [0, numControlPoints-1]
		{
		/* Find the spline segment containing the given parameter: */
		size_t segment;
		if(parameter<=Scalar(0))
			segment=0;
		else if(parameter>=Scalar(numControlPoints-1))
			segment=numControlPoints-2;
		else
			segment=size_t(Math::floor(parameter));
		parameter-=Scalar(segment);
		
		/* Evaluate the spline segment for the fractional parameter: */
		const ControlPoint& cp0=controlPoints[segment];
		Point c1=cp0.position+cp0.tangent;
		const ControlPoint& cp1=controlPoints[segment+1];
		Point c2=cp1.position-cp1.tangent;
		Point m0=Geometry::affineCombination(cp0.position,c1,parameter);
		Point m1=Geometry::affineCombination(c1,c2,parameter);
		Point m2=Geometry::affineCombination(c2,cp1.position,parameter);
		Point m3=Geometry::affineCombination(m0,m1,parameter);
		Point m4=Geometry::affineCombination(m1,m2,parameter);
		return (m4-m3)*Scalar(3);
		}
	};

#endif
