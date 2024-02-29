#include <iostream>
#include<sstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include<fstream>
using namespace std;

// MeV
const double Mn = 939.57; 
const double Mp = 938.28;
const double me = 0.511;

const double GeVtoMeV = 1000;
const double radian = M_PI / 180;
const double eps1 = 1e-3;
const double eps2 = 1e-4;

void CalcQuotients(double pnu, double theta, double* pa, double* pb, double* pc);
void CalcQuotientsMasslessCase(double pnu, double theta, double* pa0, double* pb0, double* pc0);

double Discriminant(double a, double b, double c);
void PrintQuotients(double a, double b, double c);
void Impulse(double pnu, double a, double b, double c, double* pe1, double* pe2);
void PrintImpulse(double pnu1, double pnu2);
double TotalPulse(double pp, double pe);

double ScatteringCrossSection(double pnu, double theta, double pe);
double ScatteringCrossSection(double pnu, double theta, double pe, double z);
double ScatteringCrossSectionMasslessCase(double pnu, double theta, double pe0);

double ProtonImpulse(double pnu, double pe, double theta_e);
double ProtonAngle(double pnu, double pe, double theta_e);

double func1(double pnu, double theta, double pe, double z);
double func2(double pnu, double theta, double pe, double z);
double Simpson(double pnu, double theta, double pe, double (*str)(double pnu, double theta, double pe),
	double end, int splits, double (*func)(double pnu, double theta, double pe, double z));
double SigmaLLL(double pnu, double theta, double pe);
double LowIntLim(double pnu, double theta, double pe);
double ZeroIntLim(double pnu, double theta, double pe);
double MaxImpulse(double pnu, double a, double b, double c);



//for testing
double TestStart(double a);
double TestEnd(double b);
double TestFunc(double x);
double TestSimpson(double (*str)(double a), double (*end)(double b), double n, double (*func)(double), double a, double b);




int main()
{
	// for the usual case
	double A, B, C = 0;
	double pe1, pe2 = 0;

	// for the massless case
	double A0, B0, C0 = 0;
	double pe10, pe20 = 0;
	double pe = 0; // this variable takes a positive root value in the for loop below (either pe1 or pe2)

	double pnu;

	double theta_deg;
	double theta;

	cout << "Enter the neutrino pulse (GeV): " << endl;
	cin >> pnu;

	pnu = pnu * GeVtoMeV;


	cout << "Enter the angle (degrees): " << endl;
	cin >> theta_deg;
	theta = theta_deg * radian;

	CalcQuotients(pnu, theta, &A, &B, &C);
	CalcQuotientsMasslessCase(pnu, theta, &A0, &B0, &C0);
	cout << "Discriminant is equal to (for the usual case) " << Discriminant(A, B, C) << endl << endl;
	cout << "Discriminant is equal to (for the massless case) " << Discriminant(A0, B0, C0) << endl << endl;


	//PrintQuotients(A, B, C);
	//PrintQuotients(A0, B0, C0);

	Impulse(pnu, A, B, C, &pe1, &pe2);
	Impulse(pnu, A0, B0, C0, &pe10, &pe20);

	PrintImpulse(pe1, pe2);
	//PrintImpulse(pe10, pe20);


	cout << ScatteringCrossSection(pnu, theta, pe1) << endl;

	double delta_theta = 0.25 * radian;

	int SIZE = 730;
	double* Angle = new double[730];

	double d_sigma, diff; // the difference between the normal and massless case & the fractional difference respectively
	double sigma1, sigma2, sigmalll;
	double totalEnergy;

	double protonImpulse;
	double protonAngle;

	string path1 = "DATA.txt";
	ofstream data;
	data.open(path1);

	string path2 = "PULSE.txt";
	ofstream pulse;
	pulse.open(path2);

	string path3 = "PROTON.txt";
	ofstream proton;
	proton.open(path3);

	string path4 = "ENERGY.txt";
	ofstream energy;
	energy.open(path4);

	string path5 = "DIFFERENCE.txt";
	ofstream difference;
	difference.open(path5);

	string path6 = "RAD_CORR.txt";
	ofstream rad_corr;
	rad_corr.open(path6);

	double theta_min = 0;

	double delta_theta1 = 0.25 * radian;
	double delta_theta2 = 3 * radian;

	/*for (int i = 0; i < 100; i++)
	{
		Angle[i] = 0;
	}

	for (int i = 0; Angle[i] < M_PI; i++)
	{
		if (Angle[i] <= 15 * radian)
			Angle[i] = i * delta_theta1;
		else
			Angle[i] = i * delta_theta2;
	}*/

	/*cout << i << endl;*/

	for (int i = 0; ; i++)
	{
		Angle[i] = i * delta_theta;

		CalcQuotients(pnu, Angle[i], &A, &B, &C);
		CalcQuotientsMasslessCase(pnu, Angle[i], &A0, &B0, &C0);

		Impulse(pnu, A, B, C, &pe1, &pe2);
		Impulse(pnu, A0, B0, C0, &pe10, &pe20);

		pulse << Angle[i] << "    " << pe1 << endl;
		//pulse << Angle[i] << "    " << pe1 << " | " << pe2 << endl;

		if (pe1 > 0 && pe2 < 0)
			pe = pe1;
		else
		{
			if (pe1 < 0 && pe2>0)
				pe = pe2;
			else
				cout << "Error! Both roots are either positive or negative!" << endl;
		}
		
		sigma2 = ScatteringCrossSection(pnu, Angle[i], pe);
		sigma1 = ScatteringCrossSectionMasslessCase(pnu, Angle[i], pe10);
		sigmalll = SigmaLLL(pnu, Angle[i], pe);

		d_sigma = (sigma2 - sigma1) / sigma1;
		//diff = (sigmalll - sigma2);
		diff = sigmalll;

		data << Angle[i] << " " << sigma2 << endl;
		difference << Angle[i] << " " << diff << endl;
		rad_corr << Angle[i] << " " << sigmalll << endl;

		protonImpulse = ProtonImpulse(pnu, pe, Angle[i]);
		protonAngle = ProtonAngle(pnu, pe, Angle[i]);

		totalEnergy = TotalPulse(protonImpulse, pe);
		energy << Angle[i] << "   " << totalEnergy << endl;
		// Angle[i] not protonAngle in order to draw both graphs (for an electron and a proton)
		//of the same x (in this case an angle of the eletron)
		proton << protonAngle << "    " << protonImpulse << endl;


		/*if (pnu <= pe1)
			theta_min = Angle[i];*/
		if (M_PI <= Angle[i])
			break;

		//data << Angle[i] << " " << sigma2 << endl;
	}

	pulse << endl << endl;
	/*pulse << "Mininal angle is equal to: " << theta_min << endl;*/

	data.close();
	pulse.close();
	proton.close();
	energy.close();
	difference.close();

	delete[]Angle;
	//cout << TestSimpson(TestStart, TestEnd, 100, TestFunc, -0.2, 1.82) << endl;

	return 0;
}

void PrintQuotients(double A, double B, double C)
{
	cout << "There's quotient A: " << A << endl;
	cout << "There's quotient B: " << B << endl;
	cout << "There's quotient C: " << C << endl;
}

void CalcQuotients(double pnu, double theta, double* pa, double* pb, double* pc)
{
	(*pa) = 4 * pnu * pnu - 4 * pnu * pnu * cos(theta) * cos(theta) + 8 * Mn * pnu + 4 * Mn * Mn;
	(*pb) = -4 * me * me * pnu * cos(theta) + 4 * Mp * Mp * pnu * cos(theta) - 8 * Mn * pnu * pnu * cos(theta) - 4 * Mn * Mn * pnu * cos(theta);
	(*pc) = 4 * me * me * pnu * pnu - pow(me, 4) + 2 * Mp * Mp * me * me - pow(Mp, 4) + 4 * Mn * me * me * pnu + 4 * Mn * Mp * Mp * pnu
		- 4 * Mn * Mn * pnu * pnu + 2 * Mn * Mn * me * me + 2 * Mn * Mn * Mp * Mp - 4 * pow(Mn, 3) * pnu - pow(Mn, 4);
}

void CalcQuotientsMasslessCase(double pnu, double theta, double* pa0, double* pb0, double* pc0)
{
	(*pa0) = 4 * pnu * pnu - 4 * pnu * pnu * cos(theta) * cos(theta) + 8 * Mn * pnu + 4 * Mn * Mn;
	(*pb0) = -8 * Mn * pnu * pnu * cos(theta);
	(*pc0) = -4 * Mn * Mn * pnu * pnu;
}

void Impulse(double pnu, double a, double b, double c, double* pnu1, double* pnu2)
{
	(*pnu1) = (-b + sqrt(Discriminant(a, b, c))) / (2 * a);
	(*pnu2) = (-b - sqrt(Discriminant(a, b, c))) / (2 * a);
}

void PrintImpulse(double pnu1, double pnu2)
{
	cout << "p1 = " << pnu1 << endl;
	cout << "p2 = " << pnu2 << endl;
}

double TotalPulse(double pp, double pe)
{
	double pulseTotal = pp + pe;
	return pulseTotal;
}

double Discriminant(double a, double b, double c)
{
	double discriminant = b * b - 4 * a * c;
	return discriminant;
}


double ScatteringCrossSection(double pnu, double theta, double pe) //  ds/d(cos[theta])
{
	const double Fa = 1.25, Fv = 1, Fm = 3.71, hc2 = 0.389*10E3;
	const double Gf = 1.166 * 1E-11;

	double theta_C = 13 * radian; // the Cabbibo angle

	double Ee = sqrt(pe * pe + me * me);
	double q2 = -me * me + 2 * pnu * Ee - 2 * pnu * pe * cos(theta);
	double y = q2 / (2 * Mn * pnu);

	double jacobian = ((2 * pnu * ((-1 + 2 * cos(theta)) * pow(Mn, 3) * pnu * pnu + pow((-1 +
		cos(theta)), 2) * pnu * pnu *
		Mn * pnu * (Mn + pnu) +
		Mn * Mn * ((-2 + 4 * cos(theta)) * pow(pnu, 3) + Mn * pnu * (Mn + pnu)) +
		Mn * pnu * (-pow((-1 + cos(theta)), 2) * pow(pnu, 3) +
			2 * Mn * pnu * (Mn + pnu)))) / (pow(Mn * Mn + 2 * Mn * pnu - (-1 + cos(theta) * cos(theta)) * pnu * pnu, 2)));

	double sigma = jacobian * (Gf * Gf / M_PI * cos(theta_C) * cos(theta_C) *
		(pow((Fv + Fa) / 2, 2) + pow((1 - y), 2) * pow((Fv - Fa) / 2, 2)
			+ Mn * y / (4 * pnu) * (Fa * Fa - Fv * Fv)
			+ 0.5 * y * Fm * ((1 - y) * pnu / (2 * Mn) * Fm + y * (Fv + 0.25 * Fm - Fa) + 2 * Fa)));

	return sigma * hc2;
}

double ScatteringCrossSection(double pnu, double theta, double pe, double z)
{
	const double Fa = 1.25, Fv = 1, Fm = 3.71, hc2 = 0.389 * 10E3;
	const double Gf = 1.166 * 1E-11;

	double theta_C = 13 * radian; // the Cabbibo angle

	double pe_cap = pe / z;
	double Ee_cap = sqrt(pe_cap * pe_cap + me * me); //that's capped energy
	double q2_cap = (-me * me + 2 * pnu * Ee_cap - 2 * pnu * pe_cap * cos(theta));
	double y = q2_cap / (2 * Mn * pnu);

	double jacobian = -((2 * pnu * ((-1 + 2 * cos(theta)) * pow(Mn, 3) * pnu * pnu + pow((-1 +
		cos(theta)), 2) * pnu * pnu *
		Mn * pnu * (Mn + pnu) +
		Mn * Mn * ((-2 + 4 * cos(theta)) * pow(pnu, 3) + Mn * pnu * (Mn + pnu)) +
		Mn * pnu * (-pow((-1 + cos(theta)), 2) * pow(pnu, 3) +
			2 * Mn * pnu * (Mn + pnu)))) / (pow(Mn * Mn + 2 * Mn * pnu - (-1 + cos(theta) * cos(theta)) * pnu * pnu, 2)));

	double sigma = jacobian * (Gf * Gf / M_PI * cos(theta_C) * cos(theta_C) *
		(pow((Fv + Fa) / 2, 2) + pow((1 - y), 2) * pow((Fv - Fa) / 2, 2)
			+ Mn * y / (4 * pnu) * (Fa * Fa - Fv * Fv)
			+ 0.5 * y * Fm * ((1 - y) * pnu / (2 * Mn) * Fm + y * (Fv + 0.25 * Fm - Fa) + 2 * Fa)));

	return sigma * hc2;
}

double ScatteringCrossSectionMasslessCase(double pnu, double theta, double pe0)
{
	const double Fa = 1.25, Fv = 1, Fm = 3.71, hc2 = 0.389 * 10E3; //MeV2*barn not 5
	const double Gf = 1.166 * 1E-11; //MeV-2

	double q2 = 2 * pnu * pe0 - 2 * pnu * pe0 * cos(theta);
	double y = q2 / (2 * Mn * pnu);

	double theta_C = M_PI / 180 * 13; // the Cabbibo angle

	double jacobian = -((2 * pnu * ((-1 + 2 * cos(theta)) * pow(Mn, 3) * pnu * pnu + pow((-1 +
		cos(theta)), 2) * pnu * pnu *
		Mn * pnu * (Mn + pnu) +
		Mn * Mn * ((-2 + 4 * cos(theta)) * pow(pnu, 3) + Mn * pnu * (Mn + pnu)) +
		Mn * pnu * (-pow((-1 + cos(theta)), 2) * pow(pnu, 3) +
			2 * Mn * pnu * (Mn + pnu)))) / (pow(Mn * Mn + 2 * Mn * pnu - (-1 + cos(theta) * cos(theta)) * pnu * pnu, 2)));

	double sigma = jacobian * (Gf * Gf / M_PI * cos(theta_C) * cos(theta_C) *
		(pow((Fv + Fa) / 2, 2) + pow((1 - y), 2) * pow((Fv - Fa) / 2, 2)
			+ Mn * y / (4 * pnu) * (Fa * Fa - Fv * Fv)
			+ 0.5 * y * Fm * ((1 - y) * pnu / (2 * Mn) * Fm + y * (Fv + 0.25 * Fm - Fa) + 2 * Fa)));

	return sigma * hc2;
}

double ProtonImpulse(double pnu, double pe, double theta_e)
{
	double pp = sqrt(pnu * pnu + pe * pe - 2 * pnu * pe * cos(theta_e));
	return pp;
}

double ProtonAngle(double pnu, double pe, double theta_e)
{
	double pp = ProtonImpulse(pnu, pe, theta_e);
	double theta_p = acos((pnu - pe * cos(theta_e)) / pp);
	return theta_p;
}

double func1(double pnu, double theta, double pe, double z)
{
	double jacobian = -((2 * pnu * ((-1 + 2 * cos(theta)) * pow(Mn, 3) * pnu * pnu + pow((-1 +
		cos(theta)), 2) * pnu * pnu *
		Mn * pnu * (Mn + pnu) +
		Mn * Mn * ((-2 + 4 * cos(theta)) * pow(pnu, 3) + Mn * pnu * (Mn + pnu)) +
		Mn * pnu * (-pow((-1 + cos(theta)), 2) * pow(pnu, 3) +
			2 * Mn * pnu * (Mn + pnu)))) / (pow(Mn * Mn + 2 * Mn * pnu - (-1 + cos(theta) * cos(theta)) * pnu * pnu, 2)));
	double f1 =  (ScatteringCrossSection(pnu, theta, pe, z) / z) * (1 + z * z) / (1 - z);
	return f1;
}

double func2(double pnu, double theta, double pe, double z)
{
	//jacobian sign +/-??
	double jacobian = -((2 * pnu * ((-1 + 2 * cos(theta)) * pow(Mn, 3) * pnu * pnu + pow((-1 +
		cos(theta)), 2) * pnu * pnu *
		Mn * pnu * (Mn + pnu) +
		Mn * Mn * ((-2 + 4 * cos(theta)) * pow(pnu, 3) + Mn * pnu * (Mn + pnu)) +
		Mn * pnu * (-pow((-1 + cos(theta)), 2) * pow(pnu, 3) +
			2 * Mn * pnu * (Mn + pnu)))) / (pow(Mn * Mn + 2 * Mn * pnu - (-1 + cos(theta) * cos(theta)) * pnu * pnu, 2)));
	double f2 = ScatteringCrossSection(pnu, theta, pe) * (1 + z * z) / (1 - z);
	return f2;
}


double Simpson(double pnu, double theta, double pe, double (*str)(double pnu, double theta, double pe), double end, int splits,
	double (*func)(double pnu, double theta, double pe, double z))
{
	double sum = 0;
	if (str(pnu, theta, pe) < end)
	{
		double h = (end - str(pnu, theta, pe)) / splits;

		sum = func(pnu, theta, pe, str(pnu, theta, pe)) + func(pnu, theta, pe, end);
		int k;

		for (int i = 1; i <= splits - 1; i++)
		{
			k = 2 + 2 * (i % 2);
			sum += k * func(pnu, theta, pe, str(pnu, theta, pe) + i * h);
		}
		sum *= h / 3;
	}
	return sum;
}




double SigmaLLL(double pnu, double theta, double pe)
{
	const double alpha = 0.0073; //fine structure constant
	double s = me * me + 2 * pnu * Mn + Mn * Mn; //total energy of the system in the SCM

	double E_astr = (s + me * me - Mn * Mn) / (2 * sqrt(s)); //electron energy in the SCM


	double sigmalll = ScatteringCrossSection(pnu, theta, pe) + alpha / (2 * M_PI) * log((4 * E_astr) / (me * me)) *
		(Simpson(pnu, theta, pe, LowIntLim, 0.99/*1-eps2*/, 100,
			func1) - Simpson(pnu, theta, pe, ZeroIntLim, 0.99/*1-eps2*/, 100,
				func2));
	return sigmalll;
}

double LowIntLim(double pnu, double theta, double pe) // pe1 not pe
{
	double a, b, c;
	CalcQuotients(pnu, M_PI, &a, &b, &c); //max
	double pe_max = MaxImpulse(pnu, a, b, c);

	double low = (2 * pnu * pe * (1 - cos(theta)) - me * me) / (4 * pnu * pe_max - me * me);
	//return eps2;
	return low;
}

double ZeroIntLim(double pnu, double theta, double pe)
{
	//return eps2;
	return 0;
}

double MaxImpulse(double pnu, double a, double b, double c)
{
	double pe_max = (-b + sqrt(Discriminant(a, b, c))) / (2 * a);
	return pe_max;
}

double TestStart(double a)
{
	return a;
}

double TestEnd(double b)
{
	return b;
}

double TestFunc(double x)
{
	//return cos(x);
	return sin(x);
}

double TestSimpson(double(*str)(double a), double(*end)(double b), double n, double(*func)(double), double a, double b)
{
	double h = ((*end)(b) - (*str)(a)) / n;

	double sum = func(a) + func(b);
	int k;

	for (int i = 1; i <= n - 1; i++)
	{
		k = 2 + 2 * (i % 2);
		sum += k * func(a + i * h);
	}
	sum *= h / 3;

	return sum;
}

