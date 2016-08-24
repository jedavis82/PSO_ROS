#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctime>

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <random>
#include <typeinfo>

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::getline;
using std::istringstream;
using std::ofstream;
using std::ios;

//If I make these variables static and functions inline, I can use common.hpp in multiple locations.
static string swarmOutputDir = "output/swarms/";
static string velocityOutputDir = "output/velocities/";

/**Get the present working directory **/
inline string getPWD()
{
	char *cwd;
	char buff[PATH_MAX + 1];
	cwd = getcwd(buff, PATH_MAX + 1);
	string pwd(cwd);
	return pwd;
}

/*Create the output folders*/
inline void createOutputFolders()
{
	string pwd = getPWD();
	string outputDir = pwd + "/output/";
	int dir = mkdir(outputDir.c_str(), 0777);
	swarmOutputDir = outputDir + "swarms/";
	velocityOutputDir = outputDir + "velocities/";
	dir = mkdir(swarmOutputDir.c_str(), 0777);
	dir = mkdir(velocityOutputDir.c_str(), 0777);
}

//Template functions must be defined inline when using a headless class apparently
//Convert a generic number to a string
template <typename type>
string toString(type num)
{
	std::stringstream ret;
	ret << num;
	return ret.str();
}

//Compute a random float, int, or double
template <typename type>
type computeRand(type rangeLo, type rangeHi)
{
	std::random_device rd;
	std::mt19937 generator(rd());

	if (typeid(rangeLo) == typeid(float) && typeid(rangeHi )== typeid(float))
	{
		std::uniform_real_distribution<> dis(rangeLo, rangeHi);
		return dis(generator);
	}
	else if (typeid(rangeLo) == typeid(double) && typeid(rangeHi) == typeid(double))
	{
		std::uniform_real_distribution<> dis(rangeLo, rangeHi);
		return dis(generator);
	}
	else if (typeid(rangeLo) == typeid(int) && typeid(rangeHi) == typeid(int))
	{
		std::uniform_int_distribution<> dis(rangeLo, rangeHi);
		return dis(generator);
	}
	else //Should never happen, but just in case
	{
		std::cerr << "Invalid parameters to RNG. Returning -1";
		return -1;
	}
}

/*Write the swarm vector to output file*/
template <typename PType>
void writeSwarm(int fno, std::vector<Particle<PType> > swarm, int particleSize, int swarmSize)
{
	/**Swarm File Syntax**/
	/**particle0[0], ..., particle0[n], particle0.fitness**/
	/**...**/
	/**particlen[0], ..., particlen[n], particlen.fitness**/
	ofstream outfile;
	string fn = swarmOutputDir + toString(fno) + ".txt";
	outfile.open(fn, ios::trunc | ios::out);

	for (int i = 0; i < swarmSize; i++)
	{
		for (int j = 0; j < particleSize; j++)
		{
			outfile << swarm[i].position[j] << ", ";
		}
		outfile << swarm[i].fitness;
		outfile << endl;
	}
	outfile.close();
}

/*Write the velocities out for testing*/
template <typename PType>
void writeVelocities(int fno, vector<Particle<PType> > swarm, int particleSize, int swarmSize)
{
	ofstream outfile;
	string fn = velocityOutputDir + toString(fno) + ".txt";
	outfile.open(fn, ios::trunc | ios::out);

	for (int i = 0; i < swarmSize; i++)
	{
		for (int j = 0; j < particleSize; j++)
		{
			outfile << swarm[i].velocity[j];
			if (j != particleSize - 1)
			{
				outfile << ", ";
			}
		}
		outfile << endl;
	}
}
