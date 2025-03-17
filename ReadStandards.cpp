/***************************
* ReadStandards  function  *
* read residues from files *
*  of standard amino acids *
* (c)2001 Oliver Kreylos   *
*  and Nelson Max          *
***************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "Protein.h"
#include "CreateProtein.h"

#define PI 3.14159265358979323

// #include <Geometry/MatrixRotations.h>

// #include <Geometry/AffineTransformation.h>

double residueAtomPos[25][25][3];
char residueAtomName[25][25][5];
double residued2[25];
double residued3[25];
double StandardPhi[25];
double StandardPsi[25];
double StandardAlpha[25];
double StandardBeta[25];
double StandardGamma[25];
int numbAtoms[25];
double v1[3], v2[3], v3[3], v4[3], v5[3], n1[3], n2[3], n3[3], a1[3], a2[3];
double a[3][3], b[3][3], c[3][3], d[3][3], va[3];
double temp[25][3];
double ps, ph, d2, d3, t1, t2, t3, t4, t5;

double dot(double a[3], double b[3]) {
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
  }

void normalize( double a[3]) {
  double d;
  d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  a[0] /= d;
  a[1] /= d;
  a[2] /= d;
  }
  
void cross( double a[3], double b[3], double c[3]) {
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
   }

void matmult (double a[3][3], double b[3][3], double c[3][3]) {

   int i, j, k;

   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 3; ++k)
            c[i][j] += a[i][k] * b[k][j];
         }
   }

void matmult_transp (double a[3][3], double b[3][3], double c[3][3]) {

   int i, j, k;

   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 3; ++k)
            c[i][j] += a[i][k] * b[j][k];
         }
   }

void matrix_vector (double a[3][3], double b[3], double c[3]) {

   int i, k;

   for (i = 0; i < 3; ++i) {
      c[i] = 0.;
      for (k = 0; k < 3; ++k)
         c[i] += a[i][k] *b[k];
      }
   }

void rotX(double angle, double a[3][3]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         a[i][j] = 0;
   a[1][1] = a[2][2] = cosa;
   a[1][2] = -sina;
   a[2][1] = sina;
   a[0][0] = 1.;
   }

void rotY(double angle, double a[3][3]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         a[i][j] = 0;
   a[0][0] = a[2][2] = cosa;
   a[2][0] = -sina;
   a[0][2] = sina;
   a[1][1] = 1.;
   }

void rotZ(double angle, double a[3][3]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         a[i][j] = 0;
   a[1][1] = a[0][0] = cosa;
   a[0][1] = -sina;
   a[1][0] = sina;
   a[2][2] = 1.;
   }


namespace MD {

void ReadStandards(const char* standardsDirectory) 
    {
    int i, j, k, l;
    char filename[256];
    char residueName[4], head[12];
    
    /* Build filename */
    for (i = 0; i < 25; ++i) {
	if( i == 4 || i == 8 || i == 21) continue;
	strcpy(filename, standardsDirectory);
	char buildResidueName[4];
	if(i==23)
		strcpy(buildResidueName,"Ace");
	else if(i==24)
		strcpy(buildResidueName,"Nme");
	else
		strcpy(buildResidueName,Protein::Residue::abbreviatedLcNames[i]);
	strcat(filename, buildResidueName);
	if(i == 10) strcat(filename, ".NDProtonated");
	if(i ==  5) strcat(filename, ".Thiol");
	strcat(filename, ".new.Standard");
	   
	/* Open the input file and check for validity: */
	FILE* file=fopen(filename,"rt");
	if(file==0) {
		printf("unable to open file %s\n", filename);
		continue;
		}
	
	char line[1024];
        j = -1;
	while(fgets(line,sizeof(line),file)!=0)
		{
		/* Parse the line just read: */
		sscanf(&line[0],"%11s",head);
		if(strcmp("ATOM", head)) continue;
		++j;
		sscanf(&line[12],"%4s",residueAtomName[i][j]);
		sscanf(&line[17],"%3s",residueName);
		sscanf(&line[30],"%lg",&residueAtomPos[i][j][0]);
		sscanf(&line[38],"%lg",&residueAtomPos[i][j][1]);
		sscanf(&line[46],"%lg",&residueAtomPos[i][j][2]);
		if (strcasecmp(residueName, buildResidueName)) {
			printf("Residue names do not agree: %s %s\n%d %s %s\n",
				residueName, buildResidueName,
				j, residueAtomName[i][j], filename);
			}
		}
	numbAtoms[i] = j+1;
	double d2 = 0., d3 = 0.;
	StandardGamma[i] = 59.*PI/180.;
	l = 2;
	if(i == 16) l = 1;
	
	for (k = 0; k < 3; ++k) {

// v1 is H-N  vector; trouble for proline! 
// v2 is N-CA vector; 
// v3 is CA-C vector; 
// v4 is C-O  vector; 

		v1[k] = residueAtomPos[i][1][k] - residueAtomPos[i][0][k];
		v2[k] = residueAtomPos[i][l][k] - residueAtomPos[i][0][k];
		v3[k] = residueAtomPos[i][l+2][k] - residueAtomPos[i][l][k];
		v4[k] = residueAtomPos[i][l+3][k] - residueAtomPos[i][l+2][k];
		d2 += v2[k]*v2[k];
		d3 += v3[k]*v3[k];
		v5[k] = residueAtomPos[i][0][k];
		va[k] = residueAtomPos[i][l][k];
		}
	d2 = sqrt(d2);
	d3 = sqrt(d3);
	residued2[i] = d2;
	residued3[i] = d3;

// Compute dihedral angles

		cross(v1, v2, n1);
		cross(v2, v3, n2);
		cross(v3, v4, n3);
		normalize(n1);
		normalize(n2);
		normalize(n3);
		cross(n1, n2, a1);
		cross(n2, n3, a2);
		t1 = dot(n1, n2);
		t2 = dot(n2, n3);
		t3 = dot(v2, v3);
		t4 = dot(a1, v2);
		t5 = dot(a2, v3);
		ph = acos(t1);
		ps = acos(t2);
		StandardAlpha[i] = acos(t3/(d2*d3));
		
	StandardBeta[i] = 65.*PI/180.;
	
// For proteins, dihedral angles are clockwise when viewed along bond.

	if(t4 > 0.) ph = -ph;
	if(t5 > 0.) ps = -ps;

// ph and ps differ by 180 degrees from the backbone angles.
	
	StandardPhi[i] = ph + PI;
	StandardPsi[i] = ps + PI;

#if 0
	if(i == -1)
         cout << Protein::Residue::abbreviatedNames[i]  << " " << StandardPhi[i]
           *180/PI << " " << StandardPsi[i]*180/PI  << " t2  " << t2 <<
           "  t5  " << t5 << endl;
#endif

	double angx, angy, angz;
/*
	angy = atan2(-v2[2], -v2[0]);
	angz = asin(-v2[1]/d2);
*/
	angy = atan2(v2[2], v2[0]);
	angz = asin(v2[1]/d2);

	rotY(angy, a);
	rotZ(-angz, b);
	matmult(b, a, c);
//	printf("%f  %f  v2 %f  %f  %f\n", angy, angz, v2[0], v2[1], v2[2]);

	for (j = 0; j < numbAtoms[i]; ++j) {
		for (k = 0; k < 3; ++k)
			residueAtomPos[i][j][k] -= va[k];
		matrix_vector(c, residueAtomPos[i][j], temp[j]);
		if (j < 0) printf("%d %4s %3s %f %f %f\n",
			j, residueAtomName[i][j], residueName,
			temp[j][0], temp[j][1], temp[j][2]);
	}
	for (k = 0; k < 3; ++k) 
		v4[k] = temp[l+2][k] - temp[l][k];
	angx = atan2(v4[2], v4[1]);
	rotX(-angx, d);
	for (j = 0; j < numbAtoms[i]; ++j) {
		matrix_vector(d, temp[j], v3);
		for (k = 0; k < 3; ++k)
		   residueAtomPos[i][j][k] = v3[k] - va[k];
		matrix_vector(d, temp[j], residueAtomPos[i][j]);
		if (i == -1) 
                  printf("ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
			j, residueAtomName[i][j], residueName,
			1+1, residueAtomPos[i][j][0], residueAtomPos[i][j][1],
			residueAtomPos[i][j][2], 0., 0.);
	}

	/* Clean up and return the constructed protein: */
	fclose(file);
	fflush(stdout);
	}
    }
}
