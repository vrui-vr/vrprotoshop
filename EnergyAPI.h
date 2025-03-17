/***********************************************************************
EnergyAPI - Classes to dynamically load energy computation libraries at
run-time.
Copyright (2) 2003 Oliver Kreylos
***********************************************************************/

#ifndef ENERGYAPI_INCLUDED
#define ENERGYAPI_INCLUDED

/* Forward declarations: */
namespace MD {
class Protein;
}
class EnergyCalculator;

class EnergyLibrary // Class representing dynamically loaded energy computation libraries
	{
	/* Embedded classes: */
	public:
	typedef EnergyLibrary* (*CreateEnergyLibraryFunction)(void*); // Function signature for energy library creation function inside DSO
	
	/* Elements: */
	private:
	void* dsoHandle; // Handle of dynamic shared object containing energy computation library
	
	/* Constructors and destructors: */
	protected:
	EnergyLibrary(void* sDsoHandle); // Constructs an energy library object; called from static DSO loader function
	public:
	static EnergyLibrary* load(const char* dsoName); // Opens a DSO containing an energy library and returns library object
	virtual ~EnergyLibrary(void); // Destroys library object and closes DSO
	
	/* Methods: */
	virtual int getNumEnergyComponents(void) const =0; // Returns number of energy components reported by energy library
	virtual const char* getEnergyComponentName(int energyComponentIndex) const =0; // Returns user-readable name for each energy component reported by energy library
	virtual EnergyCalculator* createEnergyCalculator(const MD::Protein* protein) =0; // Creates energy calculator for given protein
	};

class EnergyCalculator // Class representing specific energy calculator for a single protein
	{
	/* Elements: */
	private:
	const MD::Protein* protein; // Pointer to protein whose energy is to be calculated
	char* pdbFileName; // Name of temporary PDB file used to initialize energy calculation code
	
	/* Protected methods: */
	protected:
	int getNumAtoms(void) const; // Returns number of atoms in protein
	int getNumResidues(void) const; // Returns number of amino acid residues in protein
	const char* getPdbFileName(void) const // Returns name of temporary PDB file
		{
		return pdbFileName;
		};
	void getAtomCoordinates(float x[],float y[],float z[]) const; // Copies cartesian coordinates of protein's atoms into three float arrays
	void getAtomCoordinates(double x[],double y[],double z[]) const; // Same thing for double arrays
	void getAtomCoordinates(float xyz[][3]) const; // Copies cartesian coordinates of protein's atoms into 2D float array
	void getAtomCoordinates(double xyz[][3]) const; // Same thing for double array
	
	/* Constructors and destructors: */
	public:
	EnergyCalculator(const MD::Protein* sProtein); // Creates energy calculator for given protein
	virtual ~EnergyCalculator(void); // Destroys energy calculator and all resources it allocated
	
	/* Methods: */
	virtual void setEnergyComponentState(int energyComponentIndex,bool enable) =0; // Enables or disables an energy component for subsequent per-atom queries
	virtual bool getEnergyComponentState(int energyComponentIndex) const =0; // Returns current state of energy component
	virtual void updateProtein(void) =0; // Tells the energy calculator that atom coordinates have changed
	virtual double calcEnergy(void) const =0; // Calculates and returns total energy value for current atom coordinates
	virtual double getAtomEnergy(int atomIndex) const =0; // Returns per-atom energy of given atom, using current energy component states
	};

#endif
