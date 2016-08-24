#include <vector> 
using std::vector;

template<typename PType> 
struct Particle
{
	Particle<PType>() {} //default constructor

	Particle<PType>(int _particleSize)
	{
		particleSize = _particleSize;
		velocity.resize(particleSize);
		position.resize(particleSize);
		fitness = 0;
	}
	vector<PType> velocity;
	vector<PType> position;
	vector<PType> personalBest;
	int particleSize;
	PType fitness;
};