#include <iostream>
#include <complex>
#include <iomanip>
#include "matrix.h"
#include "examples.h"

using namespace std;

int main()
{	
	int k;
	cout << "//////////" << endl;
	cout << "Quantum Schrodinger Dynamics with FDTD method." << endl;
	cout << "//////////" << endl;
	
	examples::run_cn4();
	
	cout << "Simulation Complete..." << endl;
	
	cin >> k;
	return k;

}
