#pragma once
struct MolsData {
	double ch4 = 1.0;
	double n2 = 0.0;
	double co2 = 0.0;
	double c2h6 = 0.0;
	double c3h8 = 0.0;
	double ic4h10 = 0.0;
	double c4h10 = 0.0;
	double ic5h12 = 0.0;
	double c5h12 = 0.0;
	double c6h14 = 0.0;
	double c7h16 = 0.0;
	double c8h18 = 0.0;
	double c9h20 = 0.0;
	double c10h22 = 0.0;
	double h2 = 0.0;
	double o2 = 0.0;
	double co = 0.0;
	double h2o = 0.0;
	double h2s = 0.0;
	double he = 0.0;
	double ar = 0.0;

	void c6p(double c6, double p_c6h14 = 0.60, double p_c7h16 = 0.30, double p_c8h18 = 0.10) {
		c6h14 = p_c6h14 * c6;
		c7h16 = p_c7h16 * c6;
		c8h18 = p_c8h18 * c6;
	}
	double TotalData() const {
		return ch4 + n2 + co2 + c2h6 + c3h8 +
			ic4h10 + c4h10 + ic5h12 + c5h12 +
			c6h14 + c7h16 + c8h18 +
			c9h20 + c10h22 +
			h2 + o2 + co +
			h2o + h2s + he + ar;
	}

	MolsData NormalizedData() {
		double total_data = TotalData();
		MolsData normalized_data = (*this);

		normalized_data.ch4 /= total_data;
		normalized_data.n2 /= total_data;
		normalized_data.co2 /= total_data;
		normalized_data.c2h6 /= total_data;
		normalized_data.c3h8 /= total_data;
		normalized_data.ic4h10 /= total_data;
		normalized_data.c4h10 /= total_data;
		normalized_data.ic5h12 /= total_data;
		normalized_data.c5h12 /= total_data;
		normalized_data.c6h14 /= total_data;
		normalized_data.c7h16 /= total_data;
		normalized_data.c8h18 /= total_data;
		normalized_data.c9h20 /= total_data;
		normalized_data.c10h22 /= total_data;
		normalized_data.h2 /= total_data;
		normalized_data.o2 /= total_data;
		normalized_data.co /= total_data;
		normalized_data.h2o /= total_data;
		normalized_data.h2s /= total_data;
		normalized_data.he /= total_data;
		normalized_data.ar /= total_data;

		return normalized_data;
	}

	friend std::ostream& operator<<(std::ostream &out, const MolsData& mols) {
		out << "Total: " << mols.TotalData() << " (";
		out << " CH4=" << mols.ch4 << " N2=" << mols.n2 << " CO2=" << mols.co2 << " C2H6=" << mols.c2h6 << " C3H8=" << mols.c3h8;
		out << " )";
		return out;
	}
};