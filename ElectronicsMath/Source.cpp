using namespace std;
#include <stdio.h>
#include<iostream>
#include <list>
#include <bitset>
#include <math.h>
#include <string>
#include <sstream>

template < class T, class Alloc = allocator<T> > class list;
template <size_t N> class bitset;



string ITB(int num, int numBits) {
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


void DAC(double VFS, string inputCode)
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



void ADC(int numBits, double VFS, double inputVoltage) {
	cout << "Analog to Digital Converter (ADC)" << "\n";
	string output = "";

	int count = 0b11111111;
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
	cout << "Minumum error: " << minError << "\n";
	cout << "Output: " << output << "\n";
	cout << "LSB: " << VFS * pow(2.0, numBits * -1) << "\n";

}


int main()
{
	DAC(5.12, "1100010001");
	ADC(8, 5.0, 1.2);
}