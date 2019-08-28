#ifndef SAMPLE_H
#define SAMPLE_H
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include "MolsData.h"

class Sample
{
    struct Gpa{
            std::string name;
            std::string symbol;
            double btu;
            double gr;
            double summation;
            double cubeft;
       };

    public:
        Sample(MolsData &mols);
        virtual ~Sample();

        MolsData GetBTUData();
        MolsData GetGravityData();
        MolsData GetGPMData();
        
        double GetZAir();
        double GetGr();
        double GetBTU();
        double GetPB();
        double GetTB();

        //void AdjustWater(const double& wv);
        //void AdjustWater(const double& press, const double& temp);
		//double GetWaterVapor(const double& press, const double& temp);

        void SetBase(const double &press, const double &temp);
        //void SetGr(const double &gravity);
        //void SetBtu(const double &hv);
        //void SetZ(const double &compress);
        void SetGpa(const std::string &gpa_standard);
        std::vector<std::string> GetVersions();

    protected:

    private:
		//const MolsData &unnorm_mols;
		//std::unique_ptr<MolsData> unnorm_mols;
		MolsData *unnorm_mols;

		double ZSimple();
		MolsData GetSummationData();

		double zair, pb, tb;
		std::string standard;

        //std::map<std::string, std::vector<struct Gpa>> gpa;
		std::unique_ptr<std::vector<Gpa>> gpa;
        std::vector<Gpa> gpa_2145_09 = {
            { "Methane",            "CH4",      1010.0, 0.5539,     0.0116,     59.138},
            { "Nitrogen",           "N2",       0.0,    0.9672,     0.00442,    91.128},
            { "Carbon Dioxide",     "CO2",      0.0,    1.5195,     0.0195,     58.746},
            { "Ethane",             "C2H6",     1769.7, 1.0382,     0.0238,     37.488},
            { "Propane",            "C3H8",     2516.1, 1.5225,     0.0347,     36.391},
            { "i-Butane",           "iC4H10",   3251.9, 2.0068,     0.0441,     30.637},
            { "n-Butane",           "C4H10",    3262.3, 2.0068,     0.0470,     31.801},
            { "i-Pentane",          "iC5H12",   4000.9, 2.4911,     0.0576,     27.414},
            { "n-Pentane",          "C5H12",    4008.7, 2.4911,     0.0606,     27.658},
            { "n-Hexane",           "C6H14",    4755.9, 2.9754,     0.0776,     24.380},
            { "n-Heptane",          "C7H16",    5502.6, 3.4597,     0.0951,     21.730},
            { "n-Octane",           "C8H18",    6249.0, 3.9440,     0.1128,     19.570},
            { "n-Nonane",           "C9H20",    6996.3, 4.4283,     0.1307,     17.816},
            { "n-Decane",           "C10H22",   7742.9, 4.9126,     0.1556,     16.334},
            { "Hydrogen",           "H2",       0.0,    0.0,        0.0,        0.0   },
            { "Oxygen",             "O2",       0.0,    1.1048,     0.0072,     112.95},
            { "Carbon Monoxide",    "CO",       0.0,    0.0,        0.0,        0.0   },
            { "Water",              "H2O",      50.310, 0.62202,    0.06510,    175.62},
            { "Hydrogen Sulfide",   "H2S",      637.10, 1.1767,     0.0239,     74.16 },
            { "Helium",             "He",       0.0,    0.1382,     0.0,        98.693},
            { "Argon",              "Ar",       0.0,    0.0,        0.0,        0.0   }
        };
        std::vector<Gpa> gpa_2145_03 = {//No GPM for Helium, using 2145_09 value
            { "Methane",            "CH4",      1010.0, 0.55397,    0.0116,     59.137},
            { "Nitrogen",           "N2",       0.0,    0.9673,     0.00442,    91.129},
            { "Carbon Dioxide",     "CO2",      0.0,    1.5197,     0.0195,     59.095},
            { "Ethane",             "C2H6",     1769.7, 1.0383,     0.0238,     37.503},
            { "Propane",            "C3H8",     2516.2, 1.5227,     0.0349,     36.404},
            { "i-Butane",           "iC4H10",   3252.0, 2.0071,     0.0444,     30.644},
            { "n-Butane",           "C4H10",    3262.4, 2.0071,     0.0471,     31.794},
            { "i-Pentane",          "iC5H12",   4000.9, 2.4914,     0.0572,     27.390},
            { "n-Pentane",          "C5H12",    4008.7, 2.4914,     0.0603,     27.676},
            { "n-Hexane",           "C6H14",    4756.0, 2.9758,     0.0792,     24.380},
            { "n-Heptane",          "C7H16",    5502.5, 3.4601,     0.0953,     21.729},
            { "n-Octane",           "C8H18",    6248.9, 3.9445,     0.1214,     19.582},
            { "n-Nonane",           "C9H20",    6996.4, 4.4289,     0.1350,     17.807},
            { "n-Decane",           "C10H22",   7743.0, 4.9132,     0.1516,     16.323},
            { "Hydrogen",           "H2",       0.0,    0.0,        0.0,        0.0   },
            { "Oxygen",             "O2",       0.0,    1.1050,     0.0072,     112.94},
            { "Carbon Monoxide",    "CO",       0.0,    0.0,        0.0,        0.0   },
            { "Water",              "H2O",      50.312, 0.62210,    0.05557,    175.62},
            { "Hydrogen Sulfide",   "H2S",      637.11, 1.1769,     0.0242,     74.514},
            { "Helium",             "He",       0.0,    0.1382,     0.0,        98.693},
            { "Argon",              "Ar",       0.0,    0.0,        0.0,        0.0   }
        };

};


#endif // SAMPLE_H
