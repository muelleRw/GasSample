// cppSample.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Sample.h"
#include "Units.h"

int main()
{
	MolsData mol = {
		83.020,
		0.320,
		2.020,
		7.450,
		4.390,
		0.830,
		1.080, 
		0.310, 
		0.250
	};
	mol.c6p(0.300);
	mol.he = 0.030;
	Sample x(mol);
	x.SetBase(14.696, 60);
	//std::cout << "Mols: " << *(x.unnorm_mols) << std::endl;
	//std::cout << "Norms: " << x.unnorm_mols->NormalizedData() << std::endl;
	std::cout << "BTU:" << x.GetBTUData() << std::endl;
	std::cout << "Gr: " << x.GetGravityData() << std::endl;
	std::cout << "GPM: " << x.GetGPMData() << std::endl;

	mol.ch4 = 0.8;
	std::cout << "\n\n";
	std::cout << "BTU:" << x.GetBTUData() << std::endl;

    std::cout << "Hello World!\n"; 
}
