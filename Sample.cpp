#include "Sample.h"
#include <numeric>
#include <cmath>
#include <cstring>
#include <string>
#ifdef __GNUC__
#include <boost/make_unique.hpp>
#elif __MINGW32__
#include <boost/make_unique.hpp>
#endif

constexpr double GPAPRESSUREBASE = 14.696;
constexpr int WATERPOSITION = 17;
constexpr double k1 = 25.36794227;
constexpr double k2 = 7170.42747964;
constexpr double k3 = 389.5293906;
constexpr double k4 = 15.97666211;
constexpr double k5 = 7737.37631961;
constexpr double k6 = 483.28778105;
constexpr double RMULT = 10.7316;
constexpr double RANKINE = 459.67;
constexpr double MWATER = 18.0153;
constexpr double MILLION = 1000000.0;
constexpr double WVPRESSUREBASE = 14.73;
constexpr double BIAIR = 0.00537;
//constexpr double BIAIR = 0.005;

constexpr double PRESSUREBASE = 14.73;
constexpr double TEMPERATUREBASE = 60.0;
constexpr auto GPADEFAULT = "2145-09";

Sample::~Sample(){
//	delete unnorm_mols;
}

Sample::Sample(MolsData &mols){
#ifdef __GNUC__
    gpa = boost::make_unique<std::vector<Gpa>>(gpa_2145_09);
#elif _MSC_VER
    gpa = std::make_unique<std::vector<Gpa>>(gpa_2145_09);
#elif __MINGW32__
    gpa = boost::make_unique<std::vector<Gpa>>(gpa_2145_09);
#endif

    unnorm_mols = &mols;

	SetBase(PRESSUREBASE, TEMPERATUREBASE);
}

void Sample::SetBase(const double &press, const double &temp){
    pb = press;
    tb = temp;
    zair = 1 - (pb * (BIAIR * BIAIR));
}

//Create vector of individual mol BTU values and adjust for pressure base.
MolsData Sample::GetBTUData(){
	MolsData btu = unnorm_mols->NormalizedData();
	btu.ch4 = btu.ch4 * gpa->at(0).btu * (pb / GPAPRESSUREBASE);
	btu.n2 = btu.n2 * gpa->at(1).btu * (pb / GPAPRESSUREBASE);
	btu.co2 = btu.co2 * gpa->at(2).btu * (pb / GPAPRESSUREBASE);
	btu.c2h6 = btu.c2h6 * gpa->at(3).btu * (pb / GPAPRESSUREBASE);
	btu.c3h8 = btu.c3h8 * gpa->at(4).btu * (pb / GPAPRESSUREBASE);
	btu.ic4h10 = btu.ic4h10 * gpa->at(5).btu * (pb / GPAPRESSUREBASE);
	btu.c4h10 = btu.c4h10 * gpa->at(6).btu * (pb / GPAPRESSUREBASE);
	btu.ic5h12 = btu.ic5h12 * gpa->at(7).btu * (pb / GPAPRESSUREBASE);
	btu.c5h12 = btu.c5h12 * gpa->at(8).btu * (pb / GPAPRESSUREBASE);
	btu.c6h14 = btu.c6h14 * gpa->at(9).btu * (pb / GPAPRESSUREBASE);
	btu.c7h16 = btu.c7h16 * gpa->at(10).btu * (pb / GPAPRESSUREBASE);
	btu.c8h18 = btu.c8h18 * gpa->at(11).btu * (pb / GPAPRESSUREBASE);
	btu.c9h20 = btu.c9h20 * gpa->at(12).btu * (pb / GPAPRESSUREBASE);
	btu.c10h22 = btu.c10h22 * gpa->at(13).btu * (pb / GPAPRESSUREBASE);
	btu.h2 = btu.h2 * gpa->at(14).btu * (pb / GPAPRESSUREBASE);
	btu.o2 = btu.o2 * gpa->at(15).btu * (pb / GPAPRESSUREBASE);
	btu.co = btu.co * gpa->at(16).btu * (pb / GPAPRESSUREBASE);
	btu.h2o = btu.h2o * gpa->at(17).btu * (pb / GPAPRESSUREBASE);
	btu.h2s = btu.h2s * gpa->at(18).btu * (pb / GPAPRESSUREBASE);
	btu.he = btu.he * gpa->at(19).btu * (pb / GPAPRESSUREBASE);
	btu.ar = btu.ar * gpa->at(20).btu * (pb / GPAPRESSUREBASE);
    
    return btu;
}

MolsData Sample::GetGravityData(){
	MolsData gr = unnorm_mols->NormalizedData();
	gr.ch4 = gr.ch4 * gpa->at(0).gr;
	gr.n2 = gr.n2 * gpa->at(1).gr;
	gr.co2 = gr.co2 * gpa->at(2).gr;
	gr.c2h6 = gr.c2h6 * gpa->at(3).gr;
	gr.c3h8 = gr.c3h8 * gpa->at(4).gr;
	gr.ic4h10 = gr.ic4h10 * gpa->at(5).gr;
	gr.c4h10 = gr.c4h10 * gpa->at(6).gr;
	gr.ic5h12 = gr.ic5h12 * gpa->at(7).gr;
	gr.c5h12 = gr.c5h12 * gpa->at(8).gr;
	gr.c6h14 = gr.c6h14 * gpa->at(9).gr;
	gr.c7h16 = gr.c7h16 * gpa->at(10).gr;
	gr.c8h18 = gr.c8h18 * gpa->at(11).gr;
	gr.c9h20 = gr.c9h20 * gpa->at(12).gr;
	gr.c10h22 = gr.c10h22 * gpa->at(13).gr;
	gr.h2 = gr.h2 * gpa->at(14).gr;
	gr.o2 = gr.o2 * gpa->at(15).gr;
	gr.co = gr.co * gpa->at(16).gr;
	gr.h2o = gr.h2o * gpa->at(17).gr;
	gr.h2s = gr.h2s * gpa->at(18).gr;
	gr.he = gr.h2 * gpa->at(19).gr;
	gr.ar = gr.ar * gpa->at(20).gr;

	return gr;
}

MolsData Sample::GetSummationData(){
	MolsData summation = unnorm_mols->NormalizedData();
	summation.ch4 = summation.ch4 * gpa->at(0).summation;
	summation.n2 = summation.n2 * gpa->at(1).summation;
	summation.co2 = summation.co2 * gpa->at(2).summation;
	summation.c2h6 = summation.c2h6 * gpa->at(3).summation;
	summation.c3h8 = summation.c3h8 * gpa->at(4).summation;
	summation.ic4h10 = summation.ic4h10 * gpa->at(5).summation;
	summation.c4h10 = summation.c4h10 * gpa->at(6).summation;
	summation.ic5h12 = summation.ic5h12 * gpa->at(7).summation;
	summation.c5h12 = summation.c5h12 * gpa->at(8).summation;
	summation.c6h14 = summation.c6h14 * gpa->at(9).summation;
	summation.c7h16 = summation.c7h16 * gpa->at(10).summation;
	summation.c8h18 = summation.c8h18 * gpa->at(11).summation;
	summation.c9h20 = summation.c9h20 * gpa->at(12).summation;
	summation.c10h22 = summation.c10h22 * gpa->at(13).summation;
	summation.h2 = summation.h2 * gpa->at(14).summation;
	summation.o2 = summation.o2 * gpa->at(15).summation;
	summation.co = summation.co * gpa->at(16).summation;
	summation.h2o = summation.h2o * gpa->at(17).summation;
	summation.h2s = summation.h2s * gpa->at(18).summation;
	summation.he = summation.he * gpa->at(19).summation;
	summation.ar = summation.ar * gpa->at(20).summation;

	return summation;
}

MolsData Sample::GetGPMData(){
	MolsData cubeft = unnorm_mols->NormalizedData();
	double z = ZSimple();

	cubeft.ch4 = std::isnan(((cubeft.ch4 * 1000.0 / gpa->at(0).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.ch4 * 1000.0 / gpa->at(0).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.n2 = std::isnan(((cubeft.n2 * 1000.0 / gpa->at(1).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.n2 * 1000.0 / gpa->at(1).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.co2 = std::isnan(((cubeft.co2 * 1000.0 / gpa->at(2).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.co2 * 1000.0 / gpa->at(2).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c2h6 = std::isnan(((cubeft.c2h6 * 1000.0 / gpa->at(3).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c2h6 * 1000.0 / gpa->at(3).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c3h8 = std::isnan(((cubeft.c3h8 * 1000.0 / gpa->at(4).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c3h8 * 1000.0 / gpa->at(4).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.ic4h10 = std::isnan(((cubeft.ic4h10 * 1000.0 / gpa->at(5).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.ic4h10 * 1000.0 / gpa->at(5).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c4h10 = std::isnan(((cubeft.c4h10 * 1000.0 / gpa->at(6).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c4h10 * 1000.0 / gpa->at(6).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.ic5h12 = std::isnan(((cubeft.ic5h12 * 1000.0 / gpa->at(7).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.ic5h12 * 1000.0 / gpa->at(7).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c5h12 = std::isnan(((cubeft.c5h12 * 1000.0 / gpa->at(8).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c5h12 * 1000.0 / gpa->at(8).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c6h14 = std::isnan(((cubeft.c6h14 * 1000.0 / gpa->at(9).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c6h14 * 1000.0 / gpa->at(9).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c7h16 = std::isnan(((cubeft.c7h16 * 1000.0 / gpa->at(10).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c7h16 * 1000.0 / gpa->at(10).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c8h18 = std::isnan(((cubeft.c8h18 * 1000.0 / gpa->at(11).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c8h18 * 1000.0 / gpa->at(11).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c9h20 = std::isnan(((cubeft.c9h20 * 1000.0 / gpa->at(12).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c9h20 * 1000.0 / gpa->at(12).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.c10h22 = std::isnan(((cubeft.c10h22 * 1000.0 / gpa->at(13).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.c10h22 * 1000.0 / gpa->at(13).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.h2 = std::isnan(((cubeft.h2 * 1000.0 / gpa->at(14).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.h2 * 1000.0 / gpa->at(14).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.o2 = std::isnan(((cubeft.o2 * 1000.0 / gpa->at(15).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.o2 * 1000.0 / gpa->at(15).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.co = std::isnan(((cubeft.co * 1000.0 / gpa->at(16).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.co * 1000.0 / gpa->at(16).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.h2o = std::isnan(((cubeft.h2o * 1000.0 / gpa->at(17).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.h2o * 1000.0 / gpa->at(17).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.h2s = std::isnan(((cubeft.h2s * 1000.0 / gpa->at(18).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.h2s * 1000.0 / gpa->at(18).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.he = std::isnan(((cubeft.he * 1000.0 / gpa->at(19).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.he * 1000.0 / gpa->at(19).cubeft) / z) * (pb / GPAPRESSUREBASE);
	cubeft.ar = std::isnan(((cubeft.ar * 1000.0 / gpa->at(20).cubeft) / z) * (pb / GPAPRESSUREBASE)) ? 0 : ((cubeft.ar * 1000.0 / gpa->at(20).cubeft) / z) * (pb / GPAPRESSUREBASE);

	return cubeft;
}

double Sample::ZSimple() {
	MolsData summation = GetSummationData();
	return 1.0 - (pb * pow(summation.TotalData(), 2));
}
/*
void Sample::AdjustWater(const double& wv){
    if(accum > 0.0){
        if(mols[WATERPOSITION] == 0.0 && wv > 0.0){
            accum += wv;
            for(unsigned int i = 0; i < mols.size(); ++i){
                mols[i] *= (1.0 - wv / 100.0);
            }
            mols[WATERPOSITION] = wv / 100.0;
        } else if(mols[WATERPOSITION] > 0.0 && wv == 0.0){
            accum -= mols[WATERPOSITION];
            for(unsigned int i = 0; i < mols.size(); ++i){
                mols[i] /= (1.0 - mols[WATERPOSITION]);
            }
            mols[WATERPOSITION] = 0.0;

        } else if(mols[WATERPOSITION] != 0.0 && wv > 0.0){
            //Zero out Water Mol
            accum -= mols[WATERPOSITION];
            for(unsigned int i = 0; i < mols.size(); ++i){
                mols[i] /= (1.0 - mols[WATERPOSITION]);
            }
            mols[WATERPOSITION] = 0.0;
            //Add in wv
            accum += wv;
            for(unsigned int i = 0; i < mols.size(); ++i){
                mols[i] *= (1.0 - wv / 100.0);
            }
            mols[WATERPOSITION] = wv / 100.0;
        }

        ZSimple();

    } else{
        std::cout << "Mols not Set -- Unable to add Water Vapor\n";
    }
}

void Sample::AdjustWater(const double& press, const double& temp){
    AdjustWater(GetWaterVapor(press, temp) * 100.0);
}


double Sample::GetWaterVapor(const double& press, const double& temp){
    double a = exp(k1 - (k2 / (k3 + temp)));
    double b = exp(k4 - (k5 / (k6 + temp)));
    //press must be in psia -- with true atmospheric
    double w = (a / press + b) * pb / WVPRESSUREBASE;
    return w * RMULT / pb * (RANKINE + tb) / MILLION / MWATER;
}
*/

double Sample::GetZAir(){
    return zair;
}
double Sample::GetGr(){
	MolsData gr = GetGravityData();
	double z = ZSimple();
    return gr.TotalData() * zair / z;
}

double Sample::GetBTU(){
	MolsData btu = GetBTUData();
	double z = ZSimple();
	return btu.TotalData() * zair / z;
}

void Sample::SetGpa(const std::string &gpa_standard){
	if ((gpa_standard != "2145-03") && (gpa_standard != "2145-09")) {
		throw 10;
	}

	if (gpa_standard == "2145-03") {
#ifdef __GNUC__
    gpa = boost::make_unique<std::vector<Gpa>>(gpa_2145_03);
#elif _MSC_VER
    gpa = std::make_unique<std::vector<Gpa>>(gpa_2145_03);
#elif __MINGW32__
    gpa = boost::make_unique<std::vector<Gpa>>(gpa_2145_03);
#endif
	}
	else {
#ifdef __GNUC__
    gpa = boost::make_unique<std::vector<Gpa>>(gpa_2145_09);
#elif _MSC_VER
    gpa = std::make_unique<std::vector<Gpa>>(gpa_2145_09);
#elif __MINGW32__
    gpa = boost::make_unique<std::vector<Gpa>>(gpa_2145_09);
#endif
	}
}

std::vector<std::string> Sample::GetVersions(){
    std::vector<std::string> versions;
	versions.push_back("2145-09");
	versions.push_back("2145-03");
    return versions;
}

double Sample::GetPB(){
    return pb;
}
double Sample::GetTB(){
    return tb;
}


