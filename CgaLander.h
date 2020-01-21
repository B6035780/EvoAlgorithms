#ifndef CGALANDER_H
#define CGALANDER_H

//-----------------------------------------------------------------------
//  
//  Name: CgaLander.h
//  
//  Author: Mat Buckland 2002
//
//  Desc: This is the genetic algorithm class to evolve the genomes
//        used in the Lunar Lander project.
//
//------------------------------------------------------------------------


#include <vector>
#include <sstream>
#include <algorithm>

#include "defines.h"
#include "CLander.h"
#include "utils.h"

using namespace std;



//--------------------------------------------------------------
//	define our genome structure
//--------------------------------------------------------------

//--------------------------------------------------------------
//	define our genetic algorithm class
//---------------------------------------------------------------
class CgaLander
{
private:

	//our population of genomes
	
	//size of population
	
	//how many actions per chromosome

	
	void			      Mutate(vector<SGene> &vecBits);

  void            CrossoverMulti( const vector<SGene> &mum,
						                      const vector<SGene> &dad,
						                      vector<SGene>       &baby1,
						                      vector<SGene>	      &baby2);
	 
	SGenome&	      RouletteWheelSelection();

	void			      GrabNBest(int				      NBest,
							              const int	      NumCopies,
							              vector<SGenome>	&vecPop);


	void			      CreateStartPopulation();

  vector<SGenome>	Epoch();

  void            Reset();

  void            CalculateBestWorstAvTot();



public:
	
	CgaLander(	double cross_rat,
			        double mut_rat,
			        int	   pop_size,
			        int	   num_actions,
			        int	   cxClient,
			        int	   cyClient):m_dCrossoverRate(cross_rat),
							                 m_dMutationRate(mut_rat),
							                 m_iPopSize(pop_size),
							                 m_iChromoLength(num_actions),
							                 m_dTotalFitness(0),
							                 m_iGeneration(0),
                               m_dAverageFitness(0),
                               m_dWorstFitness(0),
                               m_dBestFitness(0)
					
	{
		//create the population of genomes
		CreateStartPopulation();   
	}

  //takes a population of genomes and performs one epoch.
  //and creates a new population
  void  UpdatePop(vector<SGenome> &vOldPop);

  
	
	//accessor methods
	int				       Generation()const{return m_iGeneration;}
	int				       FittestGenome()const{return m_iFittestGenome;}
  double           Best()const{return m_dBestFitness;}
  double           Average()const{return m_dAverageFitness;}
  vector<SGenome>  GrabPop()const{return m_vecPop;}
  
};



#endif