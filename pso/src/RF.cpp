#include <iostream>
#include <iomanip>
#include <vector>
#include "PSO.hpp"
#include "time.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

class RF : public PSO<float>
{
public:
	/*Construct an RF PSO instance
	_swarmSize - number of particles in the swarm
	_particleSize - number of dimensions of each particle
	_numIterations - number of times to execute the update operation
	_min - true signifies a minimization problem, false signifies a maximization problem
	_c1 - acceleration constant for the cognitive vector
	_c2 - acceleration constant for the social vector
	_nhtype - neighborhood type - valid values 'g', 'l', 'v'*/
	RF(int _swarmSize, int _particleSize, int _numIterations, bool _min, float _c1, float _c2, char _nhtype)
	{
		constructPSO(_swarmSize, _particleSize, _numIterations, _min, _c1, _c2, _nhtype);
	}


	/*Initialize the swarm of particles to random numbers in the range -5.0f to 5.0f*/
	void initialization()
	{
		//Create swarmsize particles to put in swarm
		//values should range between -5.0f and 5.0f
		float rangeLo = -5.0f, rangeHi = 5.0f;
		for (int i = 0; i < swarmSize; i++)
		{
			Particle<float> part(particleSize);
			part.personalBest.resize(particleSize);
			for (int j = 0; j < particleSize; j++)
			{
				float randPos = computeRand<float>(rangeLo, rangeHi);
				part.position[j] = randPos;
				part.personalBest[j] = part.position[j];
				part.velocity[j] = 0.0;
			}
			swarm[i] = part;
		}
	}

	/*Compute the fitness of a single particle and update its personal best if necessary.*/
	void computeFitness(Particle<float> &part)
	{
		float prevFit = part.fitness; //Will be used in the check to update the personal best
		//Compute the fitness of the particle
		float rval = 0.0;
		for (int i = 0; i < particleSize; i++)
		{
			float xval = part.position[i];
			rval += (xval * xval) - 10 * cos(2 * Pi * xval);
		}
		rval += abs(10 * particleSize);
		part.fitness = rval;

		//Update the personal best of the particle if necessary
		if (min) //If this is a minimization problem
		{
			if (rval < prevFit)
			{
				//Set the new positions of the personal best
				for (int i = 0; i < particleSize; i++)
				{
					part.personalBest[i] = part.position[i];
				}
			}
		}
	}

	//Sorting functions for the fitness values. Must be static in order to use with std::sort
	/*Sort the fitness values in increasing order*/
	template<typename PType>
	static bool sortMinFitnessValues(Particle<PType> x, Particle<PType> y)
	{
		return x.fitness > y.fitness;
	}
	/*Sort the fitness values in descreasing order*/
	template<typename PType>
	static bool sortMaxFitnessValues(Particle<PType> x, Particle<PType> y)
	{
		return x.fitness < y.fitness;
	}
};

int main(int argc, char **argv)
{
	createOutputFolders();

	const int swarmSize = 30;
	const int particleSize = 2;
	const int numIterations = 1000;
	const bool min = true; //RF is a minimization problem
	const float c1 = 0.5f;
	const float c2 = 0.5f;
	const char nhType = 'v';

	RF rf(swarmSize, particleSize, numIterations, min, c1, c2, nhType);
	clock_t start, end; //Time the code.  This will not work in Windows.
	start = clock();
	rf.execute();
	end = clock();
	double elapsedSeconds = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Finished in " << elapsedSeconds << " seconds." << endl;
}
