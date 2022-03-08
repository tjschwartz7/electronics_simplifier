#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "stdin.h"


class funcs {
public:
	template < class T, class Alloc = allocator<T> > class list;
	template <size_t N> class bitset;

	const long float k = 8.62 * pow(10, -5);

	//boltzmann const k in J / K
	const long float kJK = 1.38e-23;
	const float EGsilicon = 1.12;
	const long float Bsilicon = 1.08 * pow(10, 31);
	const long float Bge = 2.31 * pow(10, 30);
	const float EGge = .66;
	const float electronMobilitySi = 1350; // cm^2/V
	const float holeMobilitySi = 500; // cm^2/V
	const float electronCharge = 1.6 * pow(10, -19);


	const float relativePermittivity = 11.7 * 8.85e-14;
	//Nonideality factor n. Assumed to be one in our textbook, but approaches 2 on devices with high current densities.
	const int nonIdealityFactor = 1;
	const float thermalVoltageRoomTemp = .025;

	//This is only at room temp
	const float defaultSiliconNI = 1e10;

	//Chapter 1
	string ITB(int num, int numBits);
	void DAC(double VFS, string inputCode);
	void ADC(int numBits, double VFS, double inputVoltage);


	//Chapter 2

	//Drift hole elocity
	void electronHoleVelocity(double electronField, double holeMobility);

	//drift electron velocity
	void electronVelocity(double electronField, double electronMobility);
	void electricField(double voltage, double length);
	float intrinsicCarrierDensity(float T, float B, float EG);
	float conductivity(float carrierDensity, float electronVelocity);
	float electricalConductivity(float carrierDensity, float electronMobility, float holeMobility);

	//There's two different scenarios for each of these equations-
	//hole drift density, and electron drift density.
	//The only major difference (since p = n) is the electron / hole mobility which is different for each semiconductor.
	//Plugging in 500 is standard for silicon hole mobility, for example, while plugging in 1350 would be standard for silicon electron mobility.
	float driftCurrentDensity(float chargeDensity, float chargeVelocity);
	float driftCurrentDensity(float carrierDensity, float electronMobility, float electricField);
	void driftTotal(float driftDensityHole, float driftDensityElectron);


	void resistance(float length, float conductivity, float area);

	float resistivity(float conductivity);

	//Electron doping

	//ND>NA where ni is intrinsic carrier density
	void electronConcentration(float ND, float NA, float ni);

	//NA>ND where ni is intrinsic carrier density
	void holeConcentration(float NA, float ND, float ni);

	//Chapter 3
	//Diodes and diode circuits
	void differentiate(string sequence);

	//Junction potential: potential energy which exists acros the pn junction space charge
	//thermal voltage VT, acceptor electrons NA, donor electrons ND, intrinsic carrier density ni.
	float junctionPotential(float thermalVoltage, float NA, float ND, float ni);
	float thermalVoltage(int T);
	void depletionRegionWidth(float NA, float ND, float junctionPotential);

	//input in amperes and volts
	//Variables: IS, ID, VT
	void diodeVoltage(float reverseSaturationCurrent, float diodeCurrent, float thermalVoltage);
	float diodeCurrent(float reverseSaturationCurrent, float diodeVoltage, int T);
	float diodeCurrentVT(float reverseSaturationCurrent, float diodeVoltage, float thermalVoltage);
	void reverseSaturationCurrent(float diodeCurrent, float diodeVoltage, float thermalVoltage);

	//3.6- Diodes under reverse bias

	//Reverse voltage vr applied across the diode terminals is dropped across the space charge region and adds directly to the built-in potential of the junction:
	//v_j=phi_j+v_r - phi_j is electrostaticPotential. - v_j is junction voltage, v_r is reverse voltage. -> works for v_r > 0
	void junctionVoltage(float electrostaticPotential, float reverseVoltage);

	//depletion-layer width is w_d
	//Ok: xn is x when y = q*N_D; -xp is x when y =-q*N_A on a space charge density graph
	void depletionLayerWidth(float xn, float xp);

	//NA->Number of acceptor electrons - ND -> Num of donor electrons - Phi_j -> electrostatic(junction) potential - v_r->reverse voltage
	void depletionLayerWidth(float NA, float ND, float junctionPotential, float reverseVoltage);

	//depletionRegionWidth is w_do, v_r is reversevoltage, phi_j is junction potential
	void depletionLayerWidth(float depletionRegionWidth, float reverseVoltage, float junctionPotential);

	//Chapter 4 FETS
	// 
	// NMOS-
	// n-type MOS - drain current flows from northern drain
	// positive VDS faces northern drain
	// PMOS-
	// p-type MOS- drain current ID flows towards southern Drain
	// positive VDS faces southern drain
	// Universal MOS-
	// Source is grounded
	// VGS is the voltage difference from the Source(-) to the Gate(+)
	// VDS is the voltage across the Drain(+) and Source(-) terminals 
	// The voltage VGS between G and S of the FET controls the current flow ID in the terminal D.

	// NMOS 4.2--
	// 
	// -Normal-Operation-
	// VGS>VTN>0 induces an N-type channel, which will be controlled by VDS.
	// iG=0,iB=0 and ID =iS
	// When VGS << VTN, only small leakage current flows - accumulation
	// When VGS < VTN, depletion region formed under gate merges with source and drain depletion regions. No current flows between Source and Drain - Depletion
	// When VGS>VTN, channel formed between source and drain. A voltage difference has to be set between D and S for charge transportation. - Inversion
	// 
	// iS=iD ; iG = 0 ; iD->from Drain
	// 
	// 
	// -Small-VDS-
	// VGS-VTN > VDS - The device acts as a resistance and iD is proportional to (VGS-VTN)* VDS
	// Electrons move from S to D and current from D to S
	// 
	// iS=iD ; iG = 0 ; iD->from Drain
	// 
	// All inputs in V output in A
	float NMOSCurrentSmallVds(float VGS, float VTN, float VDS);
	// 
	// -VDS-Slightly-Increased
	// As VDS is increased but VGS-VTN>VDS still holds.
	// The induced channel acquires a tapered shape, and its resistance increases. 
	// The current iD is now proportional to [(VGS-VTN)VDS-VDS^2/2].
	// All inputs in V output in A
	float NMOSCurrentIncreasedVDS(float VGS, float VTN, float VDS); 
	//
	// -|VDS>=VGS-VTN|-
	// As VDS increases until VDS >= VGS-VTN 
	// the channel is pinched off. Current iD is nearly a constant and proportional to [(VGS-VTN)^2]
	//Inputs in V and output in A
	//VDS is not used for any calculation except for the initial condition. Can plug in any value VDS if you know this is the correct state ahead of time.
	float NMOSCurrentLargerVDS(float VGS, float VTN, float VDS); 

	// -|iD-VDS Characteristics|-
	//knp : transconductance parameter prime; kn = knp * WL
	//vtn : threshold voltage nmos (V)
	//vgs: voltage gate-source (V)
	// vds : voltage drain-source (V)
	//WL = width / length
	//	
	// -Triode Region-
	// VGS -VTN > VDS > 0
	// iD=Kn'W/L(VGS-VTN-VDS/2)VDS
	float NMOSCurrentTriode(double knp, double vtn, double vgs, double vds, double WL);
	float NMOSCurrentTriodekn(double kn, double vtn, double vgs, double vds);
	// 
	// -Saturation Region-
	// For VDS >= VGS-VTN
	// iD=1/2(K'n)W/L(VGS-VTN)^2
	float NMOSCurrentSaturation(double knp, double vtn, double vgs, double vds, double WL);
	float NMOSCurrentSaturationkn(double kn, double vtn, double vgs, double vds);
	//
	// -|Channel-Length modulation|-
	// As VDS increases above VDsat, length of depleted channel beyond pinch-off point, delta L, increases and actual L decreases.
	// iD increases slightly with VDS instead of being constant.
	// iD = 1/2(K'n)W/L(VGS-VTN)^2(1+lambda*VDS)
	// where lambda = channel length modulation parameter
	void NMOSCurrentCLM(float knp, float WL, float VGS, float VTN, float lambda, float VDS);
	//
	// Slope of VDS vs. iD graph: 1/r0
	// lambda = channel length modulation parameter - iD is drain current
	void NMOSr0lm(float lambda, float iD);
	//Parameter VA depends on the process technology and, for a given process, is proportional to the channel length L.
	//
	void NMOSr0VA(float VA, float iD);
	//
	// -|Body effect or substrate sensitivity|-
	// VSB!=0
	// Non-zero VSB (Source-bulk voltage(V)) changes threshold voltage, causing substrate sensitivity modeled by:
	// VTN = VTO + y((VSB+2phi_f)^1/2-(2phi_f)^1/2)
	//Where VTO = threshold voltage with zero substrate bias, y = body-effect parameter (sqrtV), and 2phi_f=surface potential parameter(V). 
	//BESS = Body Effect or Substrate Sensitivity
	void NMOS_BESS_VT(float VTO , float y, float VSB, float phi_f);
	//
	//  -|NMOS SUMMARY|-
	// All Regions: Kn = Knp*W/L = (micro_n)*(C''_ox)*W/L; iG=0, iB=0, iD=iS;
	// Cutoff Region: iD = 0 for VGS <= VTN
	// Triode Region: iD=Kn'W/L(VGS-VTN-VDS/2)VDS
	// Saturation Region: iD=1/2(K'n)W/L(VGS-VTN)^2
	// Threshold Voltage: VTN = VTO + y((VSB+2phi_f)^1/2-(2phi_f)^1/2)
	// 
	//  -|Triode Region as voltage-controlled resistor|-
	//  Output characteristics appear to be linear
	//	FET behaves like a gate-source voltage-controlled resistor between source and drain with:
	//	R_ON=[dId/dVDS]^-1=1/(Knp*W/L*(VGS-VTN))
	void NMOSTriodeRON(float KNP, float WL, float VGS, float VTN);
	// -------------------- END OF 4.2

	//4.3 PMOS Transistors
	//
	//  -|Enhancement-Mode PMOS Transistors|-
	//	P-Type source and drain regions in n-type substrate.
	//	VGS < 0 required to create p-type inversion layer in channel region.
	//	For current flow, VGS < VTP
	//	To maintain reverse bias on source-substrate and drain-substrate junctions, VSB<0 and VDB > 0.
	//	Positive bulk-source potential causes VTP to become more negative.
	//	-Output Characteristics-
	//	For VGS >= VTP, transistor is off.
	//	For more negative VGS, drain current increases in magnitude.
	//	PMOS is in triode region for small values of VDS and in saturation for larger values.
	//
	//	-|PMOS Mathematical Model Summary|-
	//	All Regions: KP = (K'P)W/L=(micro_p)(C''_ox)W/L; iG=0, iB=0, iD=iS
	//	Cutoff Region: iD = 0 for VGS >= VTP
	//	Triode Region: iD=KP(VGS-VTP-VDS/2)VDS for |VGS-VTP|>= |VDS| >= 0 
	float PMOSCurrentTriode(float KP, float VGS, float VTP, float VDS);
	//	Saturation Region: iD=1/2KP(VGS-VTP)^2(1+lambda|VDS|) for |VDS| >= |VGS-VTP| >= 0
	float PMOSCurrentSaturation(float KP, float VGS, float VTP, float lambda, float VDS);
	//  Threshold Voltage: VTP = VTO-y(sqrt(VBS+2(phi_f))-sqrt(2(phi_f)))
	void PMOS_BESS_VT(float VTO, float y, float VBS, float phi_f);
	//
	//
	void FETDrainCurrent(float kn, float vgs, float vtn);
	//------------------------------ END OF 4.3

	//4.9 Biasing the NMOS Field Effect Transistor
	//Biasing will determine the quiescent operating point (Q-POINT) which falls on a particular region of operation.
	// 
	// At triode operation: VI=VGS=VH, VO=VDS=VL
	// At cutoff: VI=VGS=Vl, VO=VDD=VH
	// 
	// Load Line: VDS = VDD - ID*RD
	void NMOSLoadLine(float VDD, float ID, float RD);
	//
	//
	//



	//Chapter 5
	// The BJT (bipolar junction transistor) is a three-terminal electronic device-
	// the output (collector) terminal i-v characteristics is controlled by the current injected into the input port (Base).
	// Both electrons and holes are the charge carriers.
	// 
	// -|npn transistor|-
	// Three terminal: North (Collector), West(Base), and South(Emitter)
	// iE=iC+iB
	// vCE=vCB+vBE
	// -|pnp transistor|-
	// Three terminal: North (Emitter), West(Base), and South(Collector)
	// iE=iC+iB
	// vEC=vEB+vBC
	// 
	// 
	 
	//-------------------------------------------------
	// -|npn Transistor: Forward Characteristics|-
	// VBE >= 0, VBC = 0
	// 
	// Base current
	// where 20<=BF<=500
	// is forward common-emitter current gain:
	void npnForwardBaseCurrent(float IF, float BF);
	float npnIF(float IS, float VBE, float VT);
	//
	//Emitter current
	// where .95 <= alpha_F = BF/(BF+1) <= 1.0
	// is forward common-base current gain
	void npnEmitterCurrent(float IC, float IB);
	void npnEmitterCurrent(float IS, float alpha_F, float VBE, float VT);
	//
	// Forward transport current
	// iC=iF=iS[exp(VBE/VT)-1]
	// where 1e-18A<=IS<=1e-9A
	//Returns iC in a forward active operation region
	void npnForwardTransportCurrent(float IF); 
	//
	// In this forward active operation region,
	// iC=BF*iB
	// iC=alpha_F*iE
	//Returns iC in a forward active operation region
	void npnForwardTransportCurrent1(float BF, float IB);
	void npnForwardTransportCurrent2(float alpha_F, float IE);
	//------------------------------------------------
	
	//------------------------------------------------
	// -|npn Transistor: Reverse Characteristics|-
	// VBE = 0, VBC>= 0
	// 
	// 0 <= BR <= 20 is reverse common-emitter current gain
	// Base currents in forward and reverse modes are different due to asymmetric doping levels in emitter and collector regions.
	// 
	// Reverse Transport current
	// iR = -iE = IS[exp(VBC/VT)-1]
	void npnReverseTransportCurrent(float iE);
	void npnReverseTransportCurrent(float IS, float VBC, float VT);
	//
	// Base current
	// iB=iR/BR=IS/BR[exp(VBC/VT)-1)]
	void npnReverseBaseCurrent(float IR, float BR);
	float npnIR(float IS, float VBC, float VT);
	//  
	//  Corrector current
	// iC = -IS/alpha_R[exp(VBC/VT)-1]
	//  where 0 <= alpha_R = BR/(BR+1) <= .95
	// is reverse common-base current gain
	void npnReverseCorrectorCurrent(float alpha_R, float IR);
	//
	//--------------------------------- 

	//---------------------------------
	//-|npn Transistors: Complete transport model equations for any bias|-
	// iC for any bias
	void npnTransportModelIC(float IS, float VBE, float VBC, float VT, float BR);
	//
	// iE for any bias
	void npnTransportModelIE(float IS, float VBE, float VBC, float VT, float BF);
	//
	// iB for any bias
	void npnTransportModelIB(float IS, float VBE, float VBC, float VT, float BF, float BR);
	void npnTransportModelIB(float IF, float BF, float IR, float BR);
	// 
	// 
	//---------------------------------

	//---------------------------------
	// -|The pnp Transistor|-
	// The three currents are determined by the two junction voltages VEB and VCB.
	// Current iE flows towards base (South)
	// Current Ic flows away from base (South)
	// Current iB flows away from base (West)
	// Voltages VEB and VCB are positive when they forward bias their respective pn junctions.
	//  
	//---------------------------------

	//---------------------------------
	// -|pnp Transistor: forward characteristics|-
	// VEB >= 0, VCB = 0
	//
	//These equations are equivalent to the npn versions just with some slight variable name mixups and one change in the emitter equation (alpha_f swapped for BF)
	void pnpForwardBaseCurrent(float IF, float BF);
	float pnpIF(float IS, float VEB, float VT);
	void pnpEmitterCurrent(float IS, float BF, float VEB, float VT);
	void pnpForwardTransportCurrent(float IF);
	// 
	//---------------------------------

	//---------------------------------
	// -|pnp Transistor: reverse characteristics|-
	// VEB = 0, VCB >= 0
	//
	//These equations are nearly identical to their npn counterparts. The collector current alpha_R was again swapped with its BR equivalency
	void pnpReverseTransportCurrent(float iE);
	float pnpIE(float IS, float VCB, float VT);
	void pnpReverseBaseCurrent(float IR, float BR);
	float pnpIR(float IS, float VCB, float VT);
	void pnpReverseCollectorCurrent(float BR, float IR);
	// 
	//---------------------------------
	
	//---------------------------------
	//-|pnp Transistor: Complete transport model equations for any bias|-
	// 
	// iC for any bias
	void pnpTransportModelIC(float IS, float VEB, float VCB, float VT, float BR);
	//
	// iE for any bias
	void pnpTransportModelIE(float IS, float VEB, float VCB, float VT, float BF);
	//
	// iB for any bias
	void pnpTransportModelIB(float IF, float BF, float IR, float BR);
	// 
	// 
	//---------------------------------

	//---------------------------------
	// -|npn Transistor: Circuit Representation for Transport Models|-
	// In npn transistor, total current traversing base is modeled by current source given by:
	// iT = iF - iR = IS * (exp(VBE/VBT)-exp(VBC/VT))
	//This equation also works with pnp, but rather VBE switches with VEB and the like.
	void npnBaseTraversalCurrent(float IF, float IR);
	//
	//This is just the calculation of the total current leaving the junction in which:
	//iB enters junction and IR/BR and IF/BF leave
	//This equation also works with pnp, but rather VBE switches with VEB and the like.
	void npnBaseCurrent(float IF, float BF, float IR, float BR);
	//---------------------------------

	//---------------------------------
	// -|pnp Transistor: Circuit Representation for Transport Models|-
	//	These are the exact same equations.
	// Just use the npn equations here.
	// For practical use, here are some tips:
	// For npn VBE is equivalent to pnp VEB
	// For npn VBC is equivalent to pnp VCB
	// These are interchangeable
	//---------------------------------

	//Leaves us off on slide 28 (Section 5.5)






	// -------------------------------
	// Reverse active: VBE < 0 : VBC > 0
	// iE=-BR*IB : iE=aR*iC : iC=(BR+1)IB
	// ------------------------------- C <--> B
	// Forward Active: VBE>0: VBC < b
	// IC=BF*IB : IE=(BF+1)IB : IC=aF*Ie
	// --------------------------------
	//Where IS is saturation current, aF is ?, aR is ?,VBE is ?, VBC is ?.
	void terminalCurrent(float IS, float aF, float aR, float VBE, float VBC);
	//iC and iE
	//BF/aF= BF/BF+1 - know one -> know other
	//VBE is usually assumed to be .7V if it is not given
	const float VBE = .7;

	//VCE=VCB+VBE
	//VCEsat=VCBsat+VBEsat

	//alpha = beta/beta+1:
	float alpha(float beta);

	//This is only at room temperature and in mV.
	const int VT = 25;

};


#endif