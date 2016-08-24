#include <iostream>
#include <random>
#include <math.h>
#include <vector>
#include <algorithm>
#include <map>
#include "Particle.hpp"
#include "common.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::map;

#define Pi 3.1415926535897L

template <typename ParticleType>
class PSO
{
public:
	//Default Constructor
	PSO<ParticleType>() {}

	/*Set up the PSO, problem independent, parameters and objects using the values passed in to the derived class*/
	void constructPSO(const int _swarmSize, const int _particleSize, const int _numIterations, const bool _min,
		const float _c1, const float _c2, const char _nhtype)
	{
		swarmSize = _swarmSize;
		numIterations = _numIterations;
		min = _min;
		particleSize = _particleSize;
		c1 = _c1;
		c2 = _c2;
		neighborhoodType = _nhtype;
		if (neighborhoodType == 'g' || neighborhoodType == 'G')
		{
			neighborhoodSize = 1; //No need to have a neighborhood. They all can see each other
		}
		else if (neighborhoodType == 'l' || neighborhoodType == 'L')
		{
			neighborhoodSize = 2; //For all instances of a local neighborhood
		}
		else if (neighborhoodType == 'v' || neighborhoodType == 'V')
		{
			//Von Neuman topologies have neighbors to the NSEW
			neighborhoodSize = 4;
			//Now I need to ensure that the swarm size is a multiple of 4.
			if (swarmSize % 4 != 0)
			{
				swarmSize = (neighborhoodSize - (swarmSize%neighborhoodSize)) + swarmSize;
			}
		}
		else
		{
			cerr << "Invalid neighborhood type.  Using global best neighborhood" << endl;
			neighborhoodSize = 1;
			neighborhoodType = 'g';
		}
		swarm.resize(swarmSize);
		neighborhoodBests.resize(swarmSize);
		constructNeighborhood();
		initialization();//Called from derived class
		for (int i = 0; i < swarmSize; i++)
		{
			computeFitness(swarm[i]); //Called from derived class
		}
		updateNeighborhoodBests<ParticleType>();
		//Set the initial mean of the population to 0
		meanPop = (float)0.0;
	}

	//Virtual function to be implemented by derived classes that initialize the swarm population.
	virtual void initialization() = 0;  //Inherited classes must implment the initialization function

	//Virtual function to be implemented by derived classes that computes the fitness of each particle in the swarm.
	virtual void computeFitness(Particle<ParticleType> &part) = 0; //Inherited classes must implement this fitness function

	/*Update the neighborhood best for each particle in the swarm*/
	template<typename PType>
	void updateNeighborhoodBests()
	{
		if (neighborhoodType == 'g' || neighborhoodType == 'G')
		{
			PType nhoodbest = swarm[0].fitness; //Because gbest particles share all information
			int index = 0;
			//Can't just set the initial index to 0.  I have to keep track of whether or not it's been changed
			bool foundBetter = false;
			for (int i = 0; i < swarmSize; i++)
			{
				PType fitval = swarm[i].fitness;
				if (min) //If this is a minimization problem
				{
					if (fitval < nhoodbest)
					{
						nhoodbest = fitval;
						index = i;
						foundBetter = true;
					}
				}
				else //If this is a maximization problem
				{
					if (fitval > nhoodbest)
					{
						nhoodbest = fitval;
						index = i;
						foundBetter = true;
					}
				}
			}
			if (foundBetter)
			{
				Particle<PType> gbest = swarm[index];
				neighborhoodBests[0] = gbest; //Because there's only one neighborhood best in gbest
			}
			else
			{
				//No better particle was found, set to the initial particle being checked.
				neighborhoodBests[0] = swarm[0];
			}
		}
		else if (neighborhoodType == 'l' || neighborhoodType == 'L')
		{
			//Neighborhood size is always 2 for lbest.
			//Loop through the swarm, compute the neighborhood best for each particle using the particle neighborhood map
			for (int i = 0; i < swarmSize; i++)
			{
				int bestIdx = i; //Initially, the particle being checked is assumed best
				PType bestFitVal = swarm[i].fitness; //Initially, the particle being checked has best fitness
				int *nhood = particleNeighborhood[i]; //Always 2 neighbors

				//Get the fitness values of the neighborhood and compare them to the idx fitness value
				for (int j = 0; j < neighborhoodSize; j++)
				{
					int idx = nhood[j];
					PType idxFitVal = swarm[idx].fitness;
					if (min) //If this is a minimization problem
					{
						if (idxFitVal < bestFitVal)
						{
							bestIdx = idx;
							bestFitVal = idxFitVal;
						}
					}
					else //If this is a maximization problem
					{
						if (idxFitVal > bestFitVal)
						{
							bestIdx = idx;
							bestFitVal = idxFitVal;
						}
					}
				}
				//Set the particle at bestIdx to be the neighborhood best for the particle at index i
				Particle<PType> nbest = swarm[bestIdx];
				neighborhoodBests[i] = nbest;
			}
		}
		else if (neighborhoodType == 'v' || neighborhoodType == 'V')
		{
			//Neighborhood size is always 4 for Von Neuman.
			//Skip any values of -1 in the neighborhood, as they denote no neighborhood connections in that direction
			//4 Neighbors in nhood: [N][S][E][W]
			//Loop through the swarm, compute the neighborhood best for each particle
			for (int i = 0; i < swarmSize; i++)
			{
				//Initially, the particle being checked is assumed best
				int bestIdx = i;
				PType bestFitVal = swarm[i].fitness;
				int *nhood = particleNeighborhood[i];
				//Get the fitness values of the neighboring particles and compare to the particle in question
				for (int j = 0; j < neighborhoodSize; j++)
				{
					int idx = nhood[j];
					if (idx == -1)//No neighbor in that direction, skip
					{
						continue;
					}
					PType idxFitVal = swarm[idx].fitness;
					if (min) //If this is a minimization problem
					{
						if (idxFitVal < bestFitVal)
						{
							bestIdx = idx;
							bestFitVal = idxFitVal;
						}
					}
					else //If this is a maximization problem
					{
						if (idxFitVal > bestFitVal)
						{
							bestIdx = idx;
							bestFitVal = idxFitVal;
						}
					}
				}
				//Set the particle at bestIdx to be the neighborhood best for the particle at index i
				Particle<PType> nbest = swarm[bestIdx];
				neighborhoodBests[i] = nbest;
			}
		}
	}
	/*Update particle velocity using the following equation:
	V[i][j](t+1) = V[i][j](t) + c1U1[j](t)(y[i][j](t) - x[i][j](t)) + c2U2[j](t)(y'[j](t) - x[i][j](t))
	This function can be overridden in derived classes.*/
	virtual void updateVelocity()
	{
		tao1 = computeRand<ParticleType>(0, 1);
		tao2 = computeRand<ParticleType>(0, 1);
		for (int i = 0; i < swarmSize; i++)
		{
			//Update velocity at element ij, i.e. velocity[i][j], using the equation above
			for (int j = 0; j < particleSize; j++)
			{
				tao1 = computeRand<ParticleType>(0, 1);
				tao2 = computeRand<ParticleType>(0, 1);
				ParticleType initVelocityVal = swarm[i].velocity[j];
				ParticleType localBestVal = c1*tao1*(swarm[i].personalBest[j] - swarm[i].position[j]);
				ParticleType neighborhoodBestVal;

				if (neighborhoodType == 'g' || neighborhoodType == 'G')
				{
					ParticleType nhbestVecVal = neighborhoodBests[0].personalBest[j];
					neighborhoodBestVal = c2*tao2*(nhbestVecVal - swarm[i].position[j]);
					int m = 0;
				}
				else
				{
					ParticleType nhbestVecVal = neighborhoodBests[i].personalBest[j];
					neighborhoodBestVal = c2*tao2*(nhbestVecVal - swarm[i].position[j]);
				}
				//Now aggregate the values...
				swarm[i].velocity[j] = initVelocityVal + localBestVal + neighborhoodBestVal;
			}
		}
	}

	/*Update particle position using the following equation:
	x[i][j](t+1) = x[i][j](t) + v[i][j](t+1)
	This function can be overridden in derived classes.*/
	virtual void updatePosition()
	{
		for (int i = 0; i < swarmSize; i++)
		{
			//Update each particle position, i.e., position[i][j]
			for (int j = 0; j < particleSize; j++)
			{
				swarm[i].position[j] = swarm[i].position[j] + swarm[i].velocity[j];
			}
		}
	}

	/*Test for convergence using the mean of the swarm.  If there is no change between iterations, the swarm is converged.
	If no isConverged() is defined in the derived class, this one will be called.*/
	virtual bool isConverged()
	{
		//Check the mean of the fitness.  If it is < x return true
		bool converged = false;
		float mean = 0;
		for (int i = 0; i < swarmSize; i++)
		{
			mean += (float)swarm[i].fitness;
		}
		mean = mean / swarmSize;
		float diff = std::abs(mean - meanPop);
		std::cout << "Population Mean: " << meanPop << " Iteration Mean: " << mean << " Mean Difference: " << diff << endl;
		if (diff < 0.0001)
		{
			converged = true;
		}
		else
		{
			meanPop = mean;
		}
		return converged;
	}

	/*Construct the neighborhood structure based on neighborhood type passed in to the derived class*/
	void constructNeighborhood()
	{
		//Set up the neighborhood based on neighborhood size and type.
		//If using global best, there is no need to construct a neighborhood since every particle can see every other particle
		if (neighborhoodType == 'g' || neighborhoodType == 'G')
			return;
		for (int i = 0; i < swarmSize; i++)
		{
			//Each particle will have an int[] of neighborhood size.
			int* nhood = new int[neighborhoodSize];
			for (int j = 0; j < particleSize; j++)
			{
				nhood[j] = -1; //Set all values to -1 upon creation.
			}
			particleNeighborhood.insert(std::pair<int, int*>(i, nhood));
			//Access each element with: particleNeighborhood[idx], will return the int*
		}
		//Now I need to set up the structure based on the neighborhood type
		if (neighborhoodType == 'l' || neighborhoodType == 'L')
		{
			//There are specific cases based on the index of the swarm.  Handle accordingly.
			//LBest Topology: (no multiples restrictions). Always neighborhood size of 2.
			/*
			*---*
			|   |
			|   |
			*---*
			*/
			for (int i = 0; i < swarmSize; i++)
			{
				//Case 1: index 0 - {i+1, n-1}, where n is swarm size
				if (i == 0)
				{
					particleNeighborhood[i][0] = (i + 1);
					particleNeighborhood[i][1] = (swarmSize - 1);
				}
				//Case n-1: index n-1 - {0, i-1}, where n is swarm size
				else if (i == swarmSize - 1)
				{
					particleNeighborhood[i][0] = 0;
					particleNeighborhood[i][1] = (i - 1);
				}
				//All vertices in between {i-1, i+1}
				else
				{
					particleNeighborhood[i][0] = (i - 1);
					particleNeighborhood[i][1] = (i + 1);
				}
			}
		}
		else if (neighborhoodType == 'v' || neighborhoodType == 'V')
		{
			//There are specific cases based on the index of the swarm.  Handle accordingly.
			//Von Neuman Topology: (Must be multiples of 4).  Always neighborhood size of 4.
			/*
			*---*---*---*
			|   |   |   |
			*---*---*---*
			|   |   |   |
			*---*---*---*
			*/
			int noEastInc = 1;
			for (int i = 0; i < swarmSize; i++)
			{
				//There are 4 cases: N S E W, nhood = [N][S][E][W]
				//Compute their values for each particle.
				//Any index of -1 indicates no neighbor in that position

				//Compute North
				if (i < neighborhoodSize - 1) //Indicates first row, no north neighbor
				{
					particleNeighborhood[i][0] = -1;
				}
				else //If the above condition is false, the particle has a north neighbor
				{
					particleNeighborhood[i][0] = i - 4;
				}
				//Compute South
				if (i >(swarmSize - neighborhoodSize - 1)) //Indicates bottom row, no south neighbor
				{
					particleNeighborhood[i][1] = -1;
				}
				else //If the above condition is false, the particle has a south neighbor
				{
					particleNeighborhood[i][1] = i + 4;
				}
				//Compute East
				if (i == (noEastInc*neighborhoodSize - 1)) //Indicates far right column, no east neighbor
				{
					particleNeighborhood[i][2] = -1;
					noEastInc++;
				}
				else //If the above condition is false, the particle has an east neighbor
				{
					particleNeighborhood[i][2] = i + 1;
				}
				//Compute West
				if ((i%neighborhoodSize) == 0) //Indicates far left column, no west neighbor
				{
					particleNeighborhood[i][3] = -1;
				}
				else //If the above condition is false, the particle has a west neighbor
				{
					particleNeighborhood[i][3] = i - 1;
				}
			}
		}
	}

	/*Run the PSO algorithm until numIterations reached, or to convergence*/
	void execute()
	{
		for (int i = 0; i < numIterations; i++)
		{
			if (i > 0)
			{
				//Compute the fitness for each particle, and update the personal best if necessary.
				for (int j = 0; j < swarmSize; j++)
				{
					computeFitness(swarm[j]);
				}
				if (isConverged())
				{
					//print the converged swarm and exit
					writeSwarm(i, swarm, particleSize, swarmSize);
					cout << "Converged at iteration " << i << " with mean " << meanPop << endl;
					cleanMemory();
					return;
				}
				updateNeighborhoodBests<ParticleType>();
			}
			updateVelocity();
			updatePosition();
			//Write the updated swarm to file
			writeSwarm(i, swarm, particleSize, swarmSize);
			writeVelocities(i, swarm, particleSize, swarmSize);
		}
		cout << "Exceeded iteration " << numIterations << ".  Final mean: " << meanPop << endl;
		cleanMemory();
	}

	/*Clean up the dynamic memory allocated*/
	void cleanMemory()
	{
		if (neighborhoodType == 'g' || neighborhoodType == 'G')
			return; //No memory was allocated
		for (int i = 0; i < swarmSize; i++)
		{
			delete[] particleNeighborhood[i];
		}
	}
	
protected:
	int swarmSize;
	int numIterations;
	int particleSize;
	bool min;
	float c1, c2, tao1, tao2;
	char neighborhoodType;
	int neighborhoodSize;
	float meanPop; //the mean will always be a floating point value
	vector<Particle<ParticleType> > swarm;
	vector<Particle<ParticleType> > neighborhoodBests;

	//The neighborhood structure is stored as a map<int,int[]> where the first int is the index
	//of the particle in question in the swarm, and the int[] are the connections to its neighbors.
	map<int, int*> particleNeighborhood;
};
