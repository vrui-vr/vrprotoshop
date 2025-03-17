/*********************************
* ReadPredictionFile - Functions *
* to read prediction files and   *
* create their protein           *
* structures                     *
* (c)2002 Nelson Max and         *
*         Oliver Kreylos         *
*********************************/

#include <stdio.h>

#include "Protein.h"
#include "CreateProtein.h"

static int strand = 0, nstrand = 0;
static int coil = 0, ncoil = 0;
static int strandstart[100], strandend[100];
static int coilstart[100], coilend[100];

#define MAXRES 1400
static int type[MAXRES] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static char pred[MAXRES] = {'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'} ;
static int numResidues = 10;
static double phi[MAXRES];
static double psi[MAXRES];
static int conf[MAXRES];

namespace MD {

extern double proteinCreatorHphi;
extern double proteinCreatorHpsi;
extern double proteinCreatorEphi;
extern double proteinCreatorEpsi;
extern double proteinCreatorCphi;
extern double proteinCreatorCpsi;
extern double proteinCreatorProlinePhi;

Protein* ReadPredictionFile(const char* predictionFilename)
	{
	int c,d,i,j;
	
	FILE* fp = fopen(predictionFilename, "r");
	if (!fp)
		{
		printf("unable to open prediction input file %s\n", predictionFilename);
		return 0;
		}
	strand = 0;
	nstrand = 0;
	coil = 0;
	ncoil = 0;
	numResidues = 0;
	while((c = fgetc(fp)) != ':')
		;
	while((c = fgetc(fp)) == ' ')
		;
	while(c >= '0' && c <= '9' && numResidues < MAXRES)
		{
		conf[numResidues] = c - '0';
		c = fgetc(fp);
		++numResidues;
		}
	i = 0;
	while((c = fgetc(fp)) != ':')
		;
	while((c = fgetc(fp)) == ' ')
		;
	while((c == 'C' || c == 'E' || c == 'H') && i < MAXRES)
		{
		pred[i] = c;
		if(c == 'C')
			{
			phi[i] = proteinCreatorCphi;
			psi[i] = proteinCreatorCpsi;
			strand = 0;
			if(!coil)
				{
				coil = 1;
				++ncoil;
				coilstart[ncoil] = i;
				}
			else
				coilend[ncoil] = i;
			}
		else if(c == 'E')
			{
			phi[i] = proteinCreatorEphi;
			psi[i] = proteinCreatorEpsi;
			coil = 0;
			if(!strand)
				{
				strand = 1;
				++nstrand;
				strandstart[nstrand] = i;
				}
			else
				strandend[nstrand] = i;
			}
		else if(c == 'H')
			{
			phi[i] = proteinCreatorHphi;
			psi[i] = proteinCreatorHpsi;
			coil = 0;
			strand = 0;
			}
		c = fgetc(fp);
		++i;
		}
	if (i != numResidues)
		{
		printf("numbers of confidences %d and predictions %d do not agree\n", numResidues, i);
		return 0;
		}
	i = 0;
	while((c = fgetc(fp)) != ':')
		;
	while((c = fgetc(fp)) == ' ')
		;
	while(c > 64 && i < MAXRES)
		{
		for( j = 0; j < 23; ++j)
			{
			d = MD::Protein::Residue::singleLetterNames[j];
			if (c == d || c + ('a' - 'A') == d)
				{
				type[i] = j;
				if(j == 16) phi[i] = proteinCreatorProlinePhi;
				break;
				}
			}
		if(j == 23)
			{
			printf("unknown single letter residue name %c\n", c);
			return 0;
			}
		c = fgetc(fp);
		++i;
		}
	if (i != numResidues)
		{
		printf("numbers of confidences %d and single letter names %d do not agree\n", numResidues, i);
		return 0;
		}
	return SetDihedrals( numResidues, type, pred, phi, psi);
	}

}
