#include "functions.h"
#include "stdin.h"

	//Chapter 1
	string funcs::ITB(int num, int numBits) {
		bool done = false;
		string ret = "";

		for (int i = numBits - 1; i >= 0; i--) {
			int binNum = (int)(pow(2.0, i));
			if (num - binNum >= 0) {
				ret += "1";
				num = num - binNum;
			}
			else {
				ret += "0";
			}
		}
		return ret;
	}


	void funcs::DAC(double VFS, string inputCode)
	{
		cout << "Digital to Analog Converter (DAC)" << "\n";
		string binary = inputCode;


		//bitset<10> bs;
		std::list<int> bits;




		double total = 0;
		int count = -1;

		for (char c : binary) {
			//bits.push_back((int)c - '0');
			total += (int)(c - '0') * pow(2.0, count--);
		}


		cout << "Vo: " << total * VFS << "V\n"; //Voltage OUT (V)
		cout << "LSB: " << VFS * pow(2.0, count + 1) << "V\n"; //Least significant bit (V)
		cout << "MSB: " << VFS * pow(2.0, -1) << "V\n"; //Most significant bit (V)
	}



	void funcs::ADC(int numBits, double VFS, double inputVoltage) {
		cout << "Analog to Digital Converter (ADC)" << "\n";
		string output = "";

		int count = pow(2, numBits) - 1;
		double minError = 100;
		for (int i = count; i >= 0; i--) {
			string bits = ITB(i, numBits);
			double error;
			double total = inputVoltage;
			double total2 = 0;
			for (int j = 1; j <= numBits; j++) {
				total2 += (int)(bits[j - 1] - '0') * pow(2.0, j * -1);
			}
			total2 *= VFS;
			total -= total2;
			error = abs(total);
			if (error < minError)
			{
				minError = error;
				output = bits;

			}
		}
		cout << "Minimum error: " << minError << "\n";
		cout << "Output: " << output << "\n";
		cout << "LSB: " << VFS * pow(2.0, numBits * -1) << "\n";

	}


	//Chapter 2

	//Drift hole elocity
	void funcs::electronHoleVelocity(double electronField, double holeMobility) {
		// electric field is in V/cm
		cout << "The velocity of the hole is : ";
		cout << scientific << holeMobility * electronField << " cm/s" << endl;
	}
	//drift electron velocity
	void funcs::electronVelocity(double electronField, double electronMobility) {
		// electric field is in V/cm
		cout << "The velocity of the electron is : ";
		cout << scientific << electronMobility * electronField << " cm/s" << endl;
	}

	void funcs::electricField(double voltage, double length) {
		//voltage in V, length in cm
		cout << "E = " << scientific << voltage / length << " V/cm" << endl;
	}

	float funcs::intrinsicCarrierDensity(float T, float B, float EG) {
		//T is absolute temperature, K
		//B is a material dependent parameter, 1.08(10^31 K^-3*cm^-6 for Si
		//EG is semiconductor bandgap energy measured in eV (electronVolts)
		long float total2 = (-EG / (k * T));
		//cout << total2 << endl;
		long float total1 = pow(T, 3.0);
		//cout << total1 << endl;
		long float n = B * total1 * exp(total2);
		cout << fixed << "Intrinsic carrier density for T= " << T << scientific << " B = " << B << fixed << " and EG = " << EG << " = " << scientific;
		cout << sqrt(n) << " /cm^3" << endl;
		return sqrt(n);
	}

	float funcs::conductivity(float carrierDensity, float electronVelocity) {
		float conductivity = funcs::electronCharge * carrierDensity * electronVelocity;
		cout << scientific << conductivity << " (ohm * cm)^-1" << endl;
		return conductivity;
	}

	float funcs::electricalConductivity(float carrierDensity, float electronMobility, float holeMobility) {
		float EC = funcs::electronCharge * (carrierDensity * (electronMobility + holeMobility));
		cout << "Electrical conductivity: " << EC << " (ohm*cm)^-1" << endl;
		return EC;
	}

	//There's two different scenarios for each of these equations-
	//hole drift density, and electron drift density.
	//The only major difference (since p = n) is the electron / hole mobility which is different for each semiconductor.
	//Plugging in 500 is standard for silicon hole mobility, for example, while plugging in 1350 would be standard for silicon electron mobility.
	float funcs::driftCurrentDensity(float chargeDensity, float chargeVelocity) {
		float j = chargeDensity * chargeVelocity;
		cout << scientific << "Drift current density: " << j << " A/cm^2" << endl;
		return j;
	}

	float funcs::driftCurrentDensity(float carrierDensity, float electronMobility, float electricField) {
		float j = carrierDensity * funcs::electronCharge * electronMobility * electricField;
		cout << scientific << "Drift current density: " << j << " A/cm^2" << endl;
		return j;
	}

	void funcs::driftTotal(float driftDensityHole, float driftDensityElectron) {
		float TD = driftDensityHole + driftDensityElectron;
		cout << "Total drift density: " << TD << " A/cm^2" << endl;
	}


	void funcs::resistance(float length, float conductivity, float area) {
		float R = length / (conductivity * area);
		cout << scientific << "Resistance: " << R << " ohm" << endl;
	}

	float funcs::resistivity(float conductivity) {
		return 1.0 / conductivity;
	}

	//Electron doping

	//ND>NA where ni is intrinsic carrier density
	void funcs::electronConcentration(float ND, float NA, float ni) {
		float n = ND - NA;
		cout << "Electron concentration: " << n << " electrons / cm^3" << endl;
		cout << "Hole concentration: " << (pow(ni, 2) / n) << " holes / cm^3" << endl;
	}

	//NA>ND where ni is intrinsic carrier density
	void funcs::holeConcentration(float NA, float ND, float ni) {
		float p = NA - ND;
		cout << "Hole concentration: " << p << " holes / cm^3" << endl;
		cout << "Electron concentration: " << (pow(ni, 2) / p) << " electrons / cm^3" << endl;
	}

	//Chapter 3
	//Diodes and diode circuits
	void funcs::differentiate(string sequence) {

	}

	//Junction potential: potential energy which exists acros the pn junction space charge
	//thermal voltage VT, acceptor electrons NA, donor electrons ND, intrinsic carrier density ni.
	float funcs::junctionPotential(float thermalVoltage, float NA, float ND, float ni) {
		float out = thermalVoltage * log((NA * ND) / pow(ni, 2));
		cout << "Junction potential: " << out << "V" << endl;
		return out;
	}

	float funcs::thermalVoltage(int T) {
		float out = (kJK * T) / electronCharge;
		cout << "Thermal Voltage: " << out << "V" << endl;
		return out;
	}

	//this is w_do
	void funcs::depletionRegionWidth(float NA, float ND, float junctionPotential) {
		float out = sqrt(((2 * relativePermittivity) / electronCharge) * (1 / NA + 1 / ND) * junctionPotential);
		cout << "Depletion region width: " << out << "m" << endl;
	}

	//input in amperes and volts
	//Variables: IS, ID, VT
	void funcs::diodeVoltage(float reverseSaturationCurrent, float diodeCurrent, float thermalVoltage) {
		float out = nonIdealityFactor * thermalVoltage * log(1 + (diodeCurrent) / (reverseSaturationCurrent));
		cout << "Diode voltage: " << out << "V" << endl;
	}

	float funcs::diodeCurrent(float reverseSaturationCurrent, float diodeVoltage, int T) {
		float out = reverseSaturationCurrent * (exp((electronCharge * diodeVoltage) / (nonIdealityFactor * kJK * T)) - 1);
		cout << "Diode current: " << out << "A" << endl;
		return out;
	}

	float funcs::diodeCurrentVT(float reverseSaturationCurrent, float diodeVoltage, float thermalVoltage) {
		float out = reverseSaturationCurrent * (exp(diodeVoltage / (nonIdealityFactor * thermalVoltage)) - 1);
		cout << "Diode current: " << out << "A" << endl;
		return out;
	}

	void funcs::reverseSaturationCurrent(float diodeCurrent, float diodeVoltage, float thermalVoltage) {
		float out = diodeCurrent / (exp(diodeVoltage / (nonIdealityFactor * thermalVoltage)) - 1);
		cout << "Reverse saturation current: " << out << "A" << endl;
	}

	void funcs::junctionVoltage(float phi_j, float reverseVoltage) {
		if (reverseVoltage <= 0) {
			cout << "Reverse voltage must be greater than zero." << endl;
			return;
		}
		float out = phi_j + reverseVoltage;
		cout << "Junction voltage: " << out << "V" << endl;
	}

	//Ok: xn is x when y = q*N_D; -xp is x when y =-q*N_A on a space charge density graph
	void funcs::depletionLayerWidth(float xn, float xp) {
		float out = xn + xp;
		cout << "Depletion-layer width: " << out << "m" << endl;
	}
	//NA->Number of acceptor electrons - ND -> Num of donor electrons - Phi_j -> electrostatic potential - v_r->reverse voltage
	void funcs::depletionLayerWidth(float NA, float ND, float junctionPotential, float reverseVoltage) {
		float out = sqrt((2*relativePermittivity/electronCharge)*(1/NA+1/ND)*(junctionPotential+reverseVoltage));
		cout << "Depletion-layer width: " << out << "m" << endl;
	}
	
	void funcs::depletionLayerWidth(float depletionRegionWidth, float reverseVoltage, float junctionPotential) {
		float out = depletionRegionWidth * sqrt(1+reverseVoltage/junctionPotential);
		cout << "Depletion-layer width: " << out << "m" << endl;
	}

	
	float funcs::NMOSCurrentTriode(double knp, double vtn, double vgs, double vds, double WL) {
		if ((vgs - vtn) > vds && (vgs - vtn) > 0)
		{
			double out = knp * WL * (vgs - vtn - (vds / 2.0)) * vds;
			cout << "Drain current (triode): " << out << "A" << endl;
			return out;
		}
		else cout << "Condition not met" << endl;
		return 0;
	}

	float funcs::NMOSCurrentSaturation(double knp, double vtn, double vgs, double vds, double WL) {
		if (vds >= vgs - vtn)
		{
			double out = (1.0 / 2.0) * knp * WL * pow(vgs - vtn, 2);
			cout << "Drain current (Saturation): " << out << "A" << endl;
			return out;
		}
		else cout << "Condition not met" << endl;
		return 0;
	}

	float funcs::NMOSCurrentTriodekn(double kn, double vtn, double vgs, double vds) {
		if ((vgs - vtn) > vds && (vgs - vtn) > 0)
		{
			double out = kn * (vgs - vtn - (vds / 2.0)) * vds;
			cout << "Drain current (triode): " << out << "A" << endl;
			return out;
		}
		else cout << "Condition not met" << endl;
		return 0;
	}

	float funcs::NMOSCurrentSaturationkn(double kn, double vtn, double vgs, double vds) {
		if ( vds >= vgs-vtn)
		{
			double out = (1.0 / 2.0) * kn * pow(vgs - vtn, 2);
			cout << "Drain current (Saturation): " << out << "A" << endl;
			return out;
		}
		else cout << "Condition not met" << endl;
		return 0;
	}

	void funcs::FETDrainCurrent(float kn, float vgs, float vtn) {
					double out = (1.0 / 2.0) * kn * pow(vgs - vtn, 2);
	}


	float funcs::NMOSCurrentSmallVds(float VGS, float VTN, float VDS) {
		if ((VGS - VTN) > VDS) {
			float out = (VGS - VTN) * VDS;
			cout << "Drain current: " << out << "A" << endl;
			return out;
		}
		else cout << "Condition VGS - VTN > VDS not met for VGS - VTN = " << VGS - VTN << "V and VDS = " << VDS << "V" << endl;
		return 0;
	}

	float funcs::NMOSCurrentIncreasedVDS(float VGS, float VTN, float VDS) {
		if ((VGS - VTN) > VDS) {
			float out = ((VGS - VTN) * VDS) - pow(VDS, 2) / 2.0;
			cout << "Drain current: " << out << "A" << endl;
			return out;
		}
		else cout << "Condition VGS - VTN > VDS not met for VGS - VTN = " << VGS - VTN << "V and VDS = " << VDS << "V" << endl;
		return 0;
	}

	float funcs::NMOSCurrentLargerVDS(float VGS, float VTN, float VDS) {
		if (VDS >= VGS-VTN) {
			float out = pow(VGS - VTN, 2);
			cout << "Drain current: " << out << "A" << endl;
			return out;
		}
		else cout << "Condition VDS >= VGS-VTN not met for VGS - VTN = " << VGS - VTN << "V and VDS = " << VDS << "V" << endl;
		return 0;
	}

	void funcs::NMOSCurrentCLM(float knp, float WL, float VGS, float VTN, float lambda, float VDS) {
		float out = 1 / 2.0 * knp * WL * pow(VGS - VTN, 2) * (1 + lambda * VDS);
		cout << "Drain current (Channel-Length Modulation): " << out << "A" << endl;
	}

	void funcs::NMOSr0lm(float lambda, float iD) {
		float out = 1 / (lambda * iD);
		cout << "R0 = " << out << "V/A(I think) and slope = 1/r0 = " << (1 / out) << "V/A" << endl;
	}
	void funcs::NMOSr0VA(float VA, float iD) {
		float out = (VA / iD);
		cout << "R0 = " << out << "V/A(I think) and slope = 1/r0 = " << (1 / out) << "V/A" << endl;
	}

	void funcs::NMOS_BESS_VT(float VTO, float y, float VSB, float phi_f) {
		if (VSB == 0) {
			cout << "VSB cannot equal zero." << endl;
			return;
		}
		float out = VTO + y*(sqrt(VSB+2*phi_f)-sqrt(2*phi_f));
		cout << "VTN equals " << VTO << "V" << endl;
	}

	void funcs::NMOSTriodeRON(float KNP, float WL, float VGS, float VTN) {
		float out = 1 / (KNP * WL * (VGS - VTN));
		cout << "Triode region - voltage-controlled resistance (NMOS): " << out << "Ohms" << endl;
	}

	float funcs::PMOSCurrentTriode(float KP, float VGS,float VTP, float VDS) {
		if (abs(VGS - VTP) >= abs(VDS) && abs(VGS) >= 0) {
			float out = KP * (VGS - VTP - (VDS / 2.0)) * VDS;
			cout << "Drain current (PMOS): " << out << "A" << endl;
			return out;
		}
		else cout << "Condition not met" << endl;
		return 0;
	}

	float funcs::PMOSCurrentSaturation(float KP, float VGS, float VTP, float lambda, float VDS) {
		if (abs(VGS) >= abs(VGS - VTP) && abs(VGS) >= 0) {
			float out = 1/2.0*KP*pow(VGS-VTP, 2)*(1+lambda*abs(VDS));
			cout << "Drain current (PMOS): " << out << "A" << endl;
			return out;
		}
		else cout << "Condition not met" << endl;
		return 0;
	}

	void funcs::PMOS_BESS_VT(float VTO, float y, float VBS, float phi_f) {
		float out = VTO - y * (sqrt(VBS + 2 * (phi_f)) - sqrt(2 * phi_f));
		cout << "Triode region - voltage-controlled resistance (PMOS): " << out << "Ohms" << endl;
	}

	void funcs::NMOSLoadLine(float VDD, float ID, float RD) {
		float out = VDD - ID * RD;
		cout << "VDS = " << out << "Therefore, the q-point (I, V) = (" << ID << "A, " << out << "V)" << endl;
	}

	void funcs::npnForwardBaseCurrent(float BF, float IF) {
		if (BF >= 20 && BF <= 500) {
			float out = IF / BF;
			cout << "Base Current = " << out << "A" << endl;
		}
		else cout << "BF out of range (20 <= BF <= 500)" << endl;
	}

	float funcs::npnIF(float IS, float VBE, float VT) {
		float out = IS * (exp(VBE / VT) - 1);
		return out;
	}

	void funcs::npnEmitterCurrent(float IC, float IB) {
		float out = IC + IB;
		cout << "Emitter current = " << out << "A" << endl;
	}

	void funcs::npnEmitterCurrent(float IS, float alpha_F, float VBE, float VT) {
		if (alpha_F >= .95) {
			float out = (IS / alpha_F) * (exp(VBE / VT) - 1);
			cout << "Emitter current = " << out << "A" << endl;
		}
		else cout << "alpha_F out of range (must be greater than .95)" << endl;
	}

	void funcs::npnForwardTransportCurrent(float IF) {
		cout << "Forward transport current = " << IF << "A" << endl;
	}

	void funcs::npnForwardTransportCurrent1(float BF, float IB) {
		float out = BF * IB;
		cout << "Forward transport current = " << out << "A" << endl;
	}
	void funcs::npnForwardTransportCurrent2(float alpha_F, float IE) {
		float out = alpha_F * IE;
		cout << "Forward transport current = " << out << "A" << endl;
	}

	float funcs::alpha(float beta) {
		return beta / (beta + 1);
	}

	void funcs::npnReverseBaseCurrent(float BR, float IR) {
		if (BR >= 20 && BR <= 500) {
			float out = IR / BR;
			cout << "Base Current = " << out << "A" << endl;
		}
		else cout << "BR out of range (20 <= BR <= 500)" << endl;
	}

	float funcs::npnIR(float IS, float VBC, float VT) {
		float out = IS * (exp(VBC / VT) - 1);
		return out;
	}

	void funcs::npnReverseTransportCurrent(float iE) {
		cout << "Reverse transport current = " << -iE << "A" << endl;
	}
	void funcs::npnReverseTransportCurrent(float IS, float VBC, float VT) {
		float out = npnIR(IS, VBC, VT);
		cout << "Reverse transport current = " << out << "A" << endl;
	}

	void funcs::npnReverseCorrectorCurrent(float alpha_R, float IR) {
		
		if (alpha_R >= 0) {
			float out = -IR / alpha_R;
			cout << "Corrector current = " << out << "A" << endl;
		}
		else cout << "Alpha_R must be greater than or equal to zero." << endl;		
	}

	void funcs::npnTransportModelIC(float IS, float VBE, float VBC, float VT, float BR) {
		float out = IS * (exp(VBE/VT)-exp(VBC/VT)) - (IS / BR) * (exp(VBC / VT) - 1);
		cout << "Collector current = " << out << "A" << endl;
	}

	void funcs::npnTransportModelIE(float IS, float VBE, float VBC, float VT, float BF) {
		float out = IS* (exp(VBE / VT) - exp(VBC / VT)) + (IS / BF) * (exp(VBE/VT)-1);
		cout << "Emitter current = " << out << "A" << endl;
	}

	void funcs::npnTransportModelIB(float IS, float VBE, float VBC, float VT, float BF, float BR) {
		float out = (IS / BF) * (exp(VBE / VT) - 1) + (IS / BR) * (exp(VBC / VT) - 1);
		cout << "Base current = " << out << "A" << endl;
	}

	void funcs::npnTransportModelIB(float IF, float BF, float IR, float BR) {
		float out = (IF / BF) + (IR / BR);
		cout << "Base current = " << out << "A" << endl;
	}

	void funcs::pnpForwardBaseCurrent(float BF, float IF) {
		if (BF >= 20 && BF <= 500) {
			float out = IF / BF;
			cout << "Base Current = " << out << "A" << endl;
		}
		else cout << "BF out of range (20 <= BF <= 500)" << endl;
	}

	float funcs::pnpIF(float IS, float VEB, float VT) {
		float out = IS * (exp(VEB / VT) - 1);
		return out;
	}

	void funcs::pnpEmitterCurrent(float IS, float BF, float VEB, float VT) {
		float out = IS * (1 + 1 / BF) * (exp(VEB / VT) - 1);
		cout << "Emitter current = " << out << "A" << endl;
		
	}

	void funcs::pnpForwardTransportCurrent(float IF) {
		cout << "Forward transport current = " << IF << "A" << endl;
	}

	void funcs::pnpReverseTransportCurrent(float iE) {
		cout << "Reverse transport current = " << -iE << "A" << endl;
	}
	float funcs::pnpIE(float IS, float VCB, float VT) {
		float out = pnpIR(IS, VCB, VT);
		return out * -1;
	}
	void funcs::pnpReverseBaseCurrent(float IR, float BR) {
		if (BR >= 20 && BR <= 500) {
			float out = IR / BR;
			cout << "Base Current = " << out << "A" << endl;
		}
		else cout << "BR out of range (20 <= BR <= 500)" << endl;
	}

	float funcs::pnpIR(float IS, float VCB, float VT) {
		float out = IS * (exp(VCB / VT) - 1);
		return out;
	}
	void funcs::pnpReverseCollectorCurrent(float BR, float IR) {
		float out = -1 * (1 + (1 / BR)) * IR;
		cout << "Collector current = " << out << "A" << endl;
	}


	void funcs::pnpTransportModelIC(float IS, float VEB, float VCB, float VT, float BR) {
		float out = IS * (exp(VEB / VT) - exp(VCB / VT)) - (IS / BR) * (exp(VCB / VT) - 1);
		cout << "Collector current = " << out << "A" << endl;
	}

	void funcs::pnpTransportModelIE(float IS, float VEB, float VCB, float VT, float BF) {
		float out = IS * (exp(VEB / VT) - exp(VCB / VT)) + (IS / BF) * (exp(VEB / VT) - 1);
		cout << "Emitter current = " << out << "A" << endl;
	}
	
	void funcs::pnpTransportModelIB(float IF, float BF, float IR, float BR) {
		float out = (IF / BF) + (IR / BR);
		cout << "Base current = " << out << "A" << endl;
	}

	void funcs::npnBaseTraversalCurrent(float IF, float IR) {
		float out = IF - IR;
		cout << "Total current traversing base: " << out << "A" << endl;
	}
	
	void funcs::npnBaseCurrent(float IF, float BF, float IR, float BR) {
		float out = IF / BF + IR / BR;
		cout << "Base current (by adding up the total current leaving the base junction) = " << out << "A" << endl;
	}