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
const double mmu = 105;

const double GeVtoMeV = 1000;
const double radian = M_PI / 180;
const double eps1 = 1e-4;
const double eps2 = 1e-2;
const string sort = "norm"; //neutrino or antineutrino


void Pulse(double pnu, double q2, double* pe);
//double ScatteringCrossSection(double pnu, double q2);
double ScatteringCrossSection(double pnu, double q2, string sort); // for testing

double ScatteringCrossSectionMasslessCase(double pnu, double q2, string sort);


double func1(double pnu, double q2, double z);
double func2(double pnu, double q2, double z);
//double Romberg(double pnu, double q2, double (*str)(double pnu, double q2),
//	double end, double (*func)(double pnu, double q2, double z));

//double Romberg_Q2(double pnu, double q2_max, string sort, string gen, double str, double end, double (*func)(double pnu, double q2, string sort, string gen));
double Simpson(double pnu, double q2, double (*str)(double pnu, double q2),
	double end, double splits, double (*func)(double pnu, double q2, double z));

//double Trapezium(double pnu, double q2[], double q2_max_curr, string sort, string gen, int corr_ord, int n);
double Trapezium(double pnu, double q2_max, string sort, string gen, int corr_ord, int n);
double Delta_diff(double pnu, string sort, int corr_ord, double q2_max);
double Delta_small(double pnu, string sort, string gen, double Q2);

double SigmaLLL(double pnu, double q2, string sort, string gen);
double LowIntLim(double pnu, double q2); /*const double q2_max*/
double ZeroIntLim(double pnu, double q2);
double SubIntLim(double pnu, double q2);


//for testing
double TestStart(double a);
double TestEnd(double b);
double TestFunc(double x);
double TestSimpson(double (*str)(double a), double (*end)(double b), double n, double (*func)(double), double a, double b);

//double R(double pnu, double q2, int n, int m, double(*func)(double pu, double q2, double z), double a, double b);
//double R_Q2(double pnu, string sort, string gen, double q2_max, int n, int m,
//	double(*func)(double pnu, double q2, string sort, string gen), int corr_ord, double a, double b);
double T(double pnu, string sort, string gen,
	double(*func)(double pnu, double q2, string sort, string gen), double q2_max, int corr_ord, int n);

double x2(double x);
double x5(double x);
double Trapeze(double a, double b, int n);
double CheckAccur();



int main()
{
	double pnu;
	

	//cout << "Enter the type of particle" << endl;
	/*cin >> sort;*/

	cout << "Enter the neutrino pulse (GeV): " << endl;
	cin >> pnu;

	pnu *= GeVtoMeV;

	//double pe;
	int steps = 20;
	double* q2 = new double[steps];


	const double q2_max = (4 * pnu * pnu * Mn) / (2 * pnu + Mn);
	const double q2_min = q2_max * 0.0001;
	double delta = (q2_max -q2_min) / steps;

	double d_sigma; 
	double sigma1, sigma2, sigmalll;
	//double totalEnergy;

	//double protonImpulse;
	//double protonAngle;

	string paths[9] = { "DATA.txt", "PULSE.txt", "PROTON.txt", "ENERGY.txt", "DIFFERENCE.txt", "RAD_CORR.txt", "DELTA.txt", "DSMALLe.txt", "DSMALLm.txt" };
	ofstream files[9];

	//cout << "x2" << "    " << CheckAccur() << endl;
	//cout << "x5" << "    " << Trapeze(0, 1, 1000) << endl;
	for (int i = 0; i < 9; i++)
		files[i].open(paths[i]);


	double theta_min = 0;

	int i = 0;
	int j = 0;
	/*files[0] << eps2 << endl << endl;
	files[0] << "Test function: " << "1 / z" << " ; " << "integration boundaries: " << "[eps2, 1]" << endl << endl;
	files[0] << "-ln(eps2)=" << -log(eps2) << endl << endl;*/



	//while (i<steps)
	//{
	//	q2[i] = i * delta;
	//	if (i == 0)
	//		q2[i] = delta / 2;
	//	
	//	sigma2 = ScatteringCrossSection(pnu, q2[i], sort);

	//	
	//	sigma1 = ScatteringCrossSectionMasslessCase(pnu, q2[i], sort);
	//	double born = sigma1;
	//	sigmalll = SigmaLLL(pnu, q2[i], sort, "electron");

	//	d_sigma = (sigma2 - sigma1) / sigma1;

	//	/*data*/ files[0] << q2[i] << " " << sigma2 << " " << sigmalll << endl;
	//	/*difference*/ files[4] << q2[i] << " " << sigmalll << endl;
	//	/*rad_corr */files[5] << q2[i] << " " << sigmalll << endl;
	//	i++;
	//}



	double pnu_min = 0.18 * GeVtoMeV;
	double pnu_max = 2 * GeVtoMeV;
	const int N = 25;
	const double pnu_step = (pnu_max - pnu_min) / N;
	double pnu_curr = pnu_min;
	double diff;
	double Q2_max;

	/*while (j < steps)
	{
		diff = Delta_diff(pnu, q2, sort, 1, q2[j]);
		files[6] << q2[j] << " " << diff << endl;
		j++;
	}*/

	while (pnu_curr <= pnu_max)
	{
		Q2_max = (4 * pnu_curr * pnu_curr * Mn) / (2 * pnu_curr + Mn);
		diff = Delta_diff(pnu_curr, sort, 1, Q2_max);
		files[6] << pnu_curr << " " << diff << endl;
		pnu_curr += pnu_step;
	}

	const int NN = 25;
	const double q2_step = (q2_max - q2_min) / NN;
	double q2_curr = q2_min;
	double delta_small;

	/*while (q2_curr <= q2_max)
	{
		delta_small = Delta_small(pnu, "anti", "electron", q2_curr);
		files[7] << q2_curr << " " << delta_small << endl;

		delta_small = Delta_small(pnu, "anti", "muon", q2_curr);
		files[8] << q2_curr << " " << delta_small << endl;

		q2_curr+=q2_step;
	}*/


	files[1]<< endl << endl;



	for (int i = 0; i < 9; i++)
		files[i].close();

	delete[]q2;

	return 0;
}

void Pulse(double pnu, double q2, double* pe)
{
	(*pe) = (-q2) / (2 * Mn) + pnu;
}


//double ScatteringCrossSection(double pnu, double q2) //  ds/d(cos[theta])
//{
//	const double Fa = 1.25, Fv = 1, Fm = 3.71, hc2 = 0.389*10E3;
//	const double Gf = 1.166 * 1E-11;
//
//	double theta_C = 13 * radian; // the Cabbibo angle
//
//	/*double pe = q2 / (2 * Mn) - pnu;
//	double Ee = sqrt(pe * pe + me * me);*/
//	double y = q2 / (2 * Mn * pnu);
//
//	double sigma = (Gf * Gf / M_PI * cos(theta_C) * cos(theta_C) *
//		(pow((Fv + Fa) / 2, 2) + pow((1 - y), 2) * pow((Fv - Fa) / 2, 2)
//			+ Mn * y / (4 * pnu) * (Fa * Fa - Fv * Fv)
//			+ 0.5 * y * Fm * ((1 - y) * pnu / (2 * Mn) * Fm + y * (Fv + 0.25 * Fm - Fa) + 2 * Fa)));
//
//	return sigma * hc2;
//}

double ScatteringCrossSection(double pnu, double q2, string sort) //for testing
{
	const double Fa = 1.25, Fv = 1, Fm = 3.71, hc2 = 0.389 * 10E18;
	const double Gf = 1.166 * 1E-11;

	double theta_C = 13 * radian; // the Cabbibo angle

	/*double Ee = sqrt(pe * pe + me * me);
	double q2 = -me * me + 2 * pnu * Ee - 2 * pnu * pe * cos(theta);*/
	double y = q2 / (2 * Mn * pnu);


	double sigma = 0;

	/*sort = "norm"; */
	if (sort == "norm")
	{
		sigma = (Gf * Gf / M_PI * cos(theta_C) * cos(theta_C) *
			(pow((Fv + Fa) / 2, 2) + pow((1 - y), 2) * pow((Fv - Fa) / 2, 2)
				+ Mn * y / (4 * pnu) * (Fa * Fa - Fv * Fv)
				+ 0.5 * y * Fm * ((1 - y) * pnu / (2 * Mn) * Fm + y * (Fv + 0.25 * Fm - Fa) + 2 * Fa)));
	}
	if (sort == "anti")
	{
		sigma = (Gf * Gf / M_PI * cos(theta_C) * cos(theta_C) *
			(pow((Fv - Fa) / 2, 2) + pow((1 - y), 2) * pow((Fv + Fa) / 2, 2)
				+ Mn * y / (4 * pnu) * (Fa * Fa - Fv * Fv)
				+ 0.5 * y * Fm * ((1 - y) * pnu / (2 * Mn) * Fm + y * (Fv + 0.25 * Fm + Fa) - 2 * Fa)));
	}
	//return sigma * hc2;
	//return 1.0;
	//cout << q2 << " " << pnu << " " << hc2 << endl;
	return sigma * hc2;
}


double ScatteringCrossSectionMasslessCase(double pnu, double q2, string sort)
{
	const double Fa = 1.25, Fv = 1, Fm = 3.71, hc2 = 0.389 * 10E3; 
	const double Gf = 1.166 * 1E-11; //MeV-2

	/*double q2 = 2 * pnu * pe0 - 2 * pnu * pe0 * cos(theta);*/
	double y = q2 / (2 * Mn * pnu);

	double theta_C = M_PI / 180 * 13; // the Cabbibo angle


	double sigma = (Gf * Gf / M_PI * cos(theta_C) * cos(theta_C) *
		(pow((Fv + Fa) / 2, 2) + pow((1 - y), 2) * pow((Fv - Fa) / 2, 2)
			+ Mn * y / (4 * pnu) * (Fa * Fa - Fv * Fv)
			+ 0.5 * y * Fm * ((1 - y) * pnu / (2 * Mn) * Fm + y * (Fv + 0.25 * Fm - Fa) + 2 * Fa)));

	return sigma * hc2;
}

double func1(double pnu, double q2, double z)
{
	double q2z = q2 / z;
	double f1 =  (ScatteringCrossSection(pnu, q2z,sort) / z) * (1 + z * z) / (1 - z); 
	//double f1 =  1 / (1 - z); 
	return f1;
}

double func2(double pnu, double q2, double z)
{
	double f2 = ScatteringCrossSection(pnu, q2, sort) * (1 + z * z) / (1 - z);
	/*double f2 =  1 / z;*/
	/*double f2 =  1 / (1 - z);*/
	/*double f2 =  1;*/
	return f2;
}

//double R(double pnu, double q2, int n, int m, double(*func)(double pu, double q2, double z), double a, double b)
//{
//	double sum = 0;
//	double h;
//	if (n == 0 && m == 0)
//	{
//		return (b - a) / 2;
//	}
//	else if (m == 0 && n != 0)
//	{
//		h = (b - a) / pow(2, n);
//		for (int i = 1; i < pow(2, n - 1); i++)
//		{
//			sum += func(pnu, q2, a + (2 * i - 1) * h);
//		}
//		return R(pnu, q2, n - 1, 0, func, a, b) / 2 + sum * h;
//	}
//	else
//	{
//		return 1 / (pow(4, m) - 1) * (pow(4, m) * R(pnu, q2, n, m - 1, func, a, b) - R(pnu, q2, n - 1, m - 1, func, a, b));
//	}
//}

//double R_Q2(double pnu, string sort, string gen, double q2_max, int n, int m,
//	double(*func)(double pnu, double q2, string sort, string gen), int corr_ord, double a, double b)
//{
//	double sum = 0;
//	double h;
//
//
//	if (corr_ord == 1)
//	{
//		if (n == 0 && m == 0)
//		{
//			return (b - a) / 2;
//		}
//		else if (m == 0 && n != 0)
//		{
//			h = (b - a) / pow(2, n);
//			for (int i = 1; i < pow(2, n - 1); i++)
//			{
//				sum += func(pnu, q2, a + (2 * i - 1) * h);
//			}
//			return R(pnu, q2, n - 1, 0, func, a, b) / 2 + sum * h;
//		}
//		else
//		{
//			return 1 / (pow(4, m) - 1) * (pow(4, m) * R(pnu, q2, n, m - 1, func, a, b) - R(pnu, q2, n - 1, m - 1, func, a, b));
//		}
//	}
//	
//}

double T(double pnu, string sort, string gen, double(*func)(double pnu, double q2, string sort, string gen), 
	double q2_max, int corr_ord, int n)
{
	double q2_min = 5000;
	double h = (q2_max - q2_min) / n;
	double sum = 0.0;

	if (corr_ord == 1)
	{
		//sum = SigmaLLL(pnu, q2_min, sort, gen) + SigmaLLL(pnu, q2_max, sort, gen);
		sum = func(pnu, q2_min, sort, gen) + func(pnu, q2_max, sort, gen);
		for (int i = 1; i <= n - 1; i++)
			//sum += 2 * SigmaLLL(pnu, q2_min + i * h, sort, gen);
			sum += 2 * func(pnu, q2_min + i * h, sort, gen);
		sum *= h / 2;
	}
	return sum;
}

double x2(double x)
{
	return x*x;
}

double x5(double x)
{
	return pow(x,5);
}

double Trapeze(double a, double b, int n)
{
	double h = (b - a) / n;
	double sum = x2(a) + x2(b);

	for (int i = 0; i <= n - 1; i++) {
		sum += 2 * x2(a + i * h);
	}
	sum *= h / 2;

	return sum;
}

double CheckAccur()
{
	double precise_val = 0.3333333333;
	//double num_val = Trapeze(0, 1, 10);


	double err = 1e-8;
	double res, res1;
	int i = 1;
	double d = 0.0;

	res = Trapeze(0, 1, 5);
	do
	{
		//res1 = res;
		res = Trapeze(0, 1, 5 * (i + 1));
		d = res - precise_val;
		i++;
	} while (abs(d) > err);

	return res;

}

//double Romberg(double pnu, double q2, double (*str)(double pnu, double q2), double end, 
//	double (*func)(double pnu, double q2, double z))
//{
//	double err = 1e-4;
//
//	double d, res, res1;
//	int i = 1;
//
//	if (str(pnu, q2) < end)
//	{
//		res = R(pnu, q2, 0, 0, func, str(pnu, q2), end);
//		int i = 1;
//		do {
//			res1 = res;
//			res = R(pnu, q2, i, 0, func, str(pnu, q2), end);
//			d = res - res1;
//			i++;
//		} while (abs(d) > err);
//		return res;
//	}
//}

//double Romberg_Q2(double pnu, double q2_max, string sort, string gen, double str, double end, double (*func)(double pnu, double q2, string sort, string gen))
//{
//	double err = 1e-4;
//
//	double d, res, res1;
//	int i = 1;
//
//	if (str < end)
//	{
//		res = R_Q2(pnu, q2_max, 0, 0, func, str, end);
//		int i = 1;
//		do {
//			res1 = res;
//			res = R_Q2(pnu, q2_max, i, 0, func, str, end);
//			d = res - res1;
//			i++;
//		} while (abs(d) > err);
//		return res;
//	}
//}

double Simpson(double pnu, double q2, double(*str)(double pnu, double q2), double end, double splits, double(*func)(double pnu, double q2, double z))
{
	//double sum1 = 0;
	double sum = 0;

	if (str(pnu, q2) < end)
	{
		double h = (end - str(pnu, q2)) / splits;

		sum = func(pnu, q2, str(pnu, q2)) + func(pnu, q2, end);

		for (int i = 1; i <= splits - 1; i++)
		{
			if (i % 3 == 0)
				sum += 2 * func(pnu, q2, str(pnu, q2) + i * h);
			else
				sum += 3 * func(pnu, q2, str(pnu, q2) + i * h);
		}
		sum *= 3 * h / 8;


		/*double h2 = (end - str(pnu, q2)) / splits2;

		sum2 = func(pnu, q2, str(pnu, q2)) + func(pnu, q2, end);

		for (int i = 1; i <= splits2 - 1; i++)
		{
			if (i % 3 == 0)
				sum2 += 2 * func(pnu, q2, str(pnu, q2) + i * h2);
			else
				sum2 += 3 * func(pnu, q2, str(pnu, q2) + i * h2);
		}
		sum2 *= 3 * h2 / 8;*/


	/*	do {
			sum2 = sum1;

			splits1 *= 2;


			diff = sum2 - sum1;
		} while (abs(diff) > err);*/

	}
	return sum;
}

//double Trapezium(double pnu, double q2[], double q2_max_curr, string sort, string gen, int corr_ord, int n)
double Trapezium(double pnu, double q2_max, string sort, string gen, int corr_ord, int n)
{
	//double h = (q2_max_curr - q2[0]) / n;
	//double sum = 0.0;


	//if (corr_ord == 0)
	//{
	//	sum = ScatteringCrossSection(pnu, q2[0], sort) + ScatteringCrossSection(pnu, q2_max_curr, sort);
	//	for (int i = 1; i <= n - 1; i++)
	//	{
	//		//sum += 2 * ScatteringCrossSection(pnu, q2[i], sort);
	//		sum += 2 * ScatteringCrossSection(pnu, q2[0] + i * h, sort);
	//	}
	//	sum *= h / 2;
	//}

	//if (corr_ord == 1)
	//{
	//	sum = SigmaLLL(pnu, q2[0], sort, gen) + SigmaLLL(pnu, q2_max_curr, sort, gen);
	//	for (int i = 1; i <= n - 1; i++)
	//	{
	//		//sum += 2 * SigmaLLL(pnu, q2[i], sort, gen);
	//		sum += 2 * SigmaLLL(pnu, q2[0] + i * h, sort, gen);
	//	}
	//	sum *= h / 2;
	//}

	//return sum;

	double err = 1e-3;
	double res, res1;
	int i = 1;
	double d = 0.0;

	res = T(pnu, sort, gen, SigmaLLL, q2_max, corr_ord, n);
	do
	{
		res1 = res;
		res = T(pnu, sort, gen, SigmaLLL, q2_max, corr_ord, n * (i + 1));
		d = res - res1;
		i++;
	} while (abs(d) > err);

	return res;
}

double Delta_diff(double pnu, string sort, int corr_ord, double q2_max)
{
	double diff = 0.0;  

	if (corr_ord == 0)
	{
		diff = 0;
	}

	if (corr_ord == 1)
		diff = Trapezium(pnu, q2_max, sort, "muon", 1, 15) / Trapezium(pnu, q2_max, sort, "electron", 1, 15) - 1;
		//diff = Romberg_Q2(pnu, q2_max, sort, "muon", 1, 25) / Romberg_Q2(pnu, q2_max, sort, "electron", 1, 25) - 1;

	return diff;
}

double Delta_small(double pnu, string sort, string gen, double Q2)
{
	double born = ScatteringCrossSection(pnu, Q2, sort);
	double corr = SigmaLLL(pnu, Q2, sort, gen);

	double delta = (corr / born - 1) * 100;
	return delta;
}


double SigmaLLL(double pnu, double q2, string sort, string gen)
{
	const double alpha = 0.0073; //fine structure constant
	double E_astr = 0.0;
	double s = 0.0;
	double sigmalll = 0.0;

	if (gen == "electron")
	{
		s = me * me + 2 * pnu * Mn + Mn * Mn; //total energy of the system in the SCM
		E_astr = (s + me * me - Mn * Mn) / (2 * sqrt(s)); //electron energy in the SCM


		//sigmalll = ScatteringCrossSection(pnu, q2, sort) + alpha / (2 * M_PI) * log((4 * E_astr * E_astr) / (me * me)) *
		//	(Simpson(pnu, q2, LowIntLim, 1 - eps2, 1E3 - 1,
		//		func1) + Simpson(pnu, q2, SubIntLim, 1 - eps1, 1E6 - 1,
		//			func1) - Simpson(pnu, q2, ZeroIntLim, 1 - eps2, 1E3 - 1,
		//				func2) - Simpson(pnu, q2, SubIntLim, 1 - eps1, 1E6 - 1,
		//					func2));

		sigmalll = ScatteringCrossSection(pnu, q2, sort) + alpha / (2 * M_PI) * log((4 * E_astr * E_astr) / (me * me)) *
			(Simpson(pnu, q2, LowIntLim, 1 - eps2, 666,
				func1) + Simpson(pnu, q2, SubIntLim, 1 - eps1, 199998,
					func1) - Simpson(pnu, q2, ZeroIntLim, 1 - eps2, 666,
						func2) - Simpson(pnu, q2, SubIntLim, 1 - eps1, 199998,
							func2));
	}

	if (gen == "muon")
	{
		s = mmu * mmu + 2 * pnu * Mn + Mn * Mn; //total energy of the system in the SCM
		E_astr = (s + mmu * mmu - Mn * Mn) / (2 * sqrt(s)); //muon energy in the SCM


		/*sigmalll = ScatteringCrossSection(pnu, q2, sort) + alpha / (2 * M_PI) * log((4 * E_astr * E_astr) / (mmu * mmu)) * 
			(Simpson(pnu, q2, LowIntLim, 1 - eps2, 1E3 - 1,
				func1) + Simpson(pnu, q2, SubIntLim, 1 - eps1, 1E6 - 1,
					func1) - Simpson(pnu, q2, ZeroIntLim, 1 - eps2, 1E3 - 1,
						func2) - Simpson(pnu, q2, SubIntLim, 1 - eps1, 1E6 - 1,
							func2));*/

		sigmalll = ScatteringCrossSection(pnu, q2, sort) + alpha / (2 * M_PI) * log((4 * E_astr * E_astr) / (mmu * mmu)) *
			(Simpson(pnu, q2, LowIntLim, 1 - eps2, 666,
				func1) + Simpson(pnu, q2, SubIntLim, 1 - eps1, 199998,
					func1) - Simpson(pnu, q2, ZeroIntLim, 1 - eps2, 666,
						func2) - Simpson(pnu, q2, SubIntLim, 1 - eps1, 199998,
							func2));
	}

	return sigmalll;
}

double LowIntLim(double pnu, double q2) 
{
	double q2_max = (4 * pnu * pnu * Mn) / (2 * pnu + Mn);
	double low = q2 / q2_max;
	if (low > 1 - eps2)
		low = 1 - eps2;
	//double a, b, c;
	//CalcQuotients(pnu, M_PI, &a, &b, &c); //max
	//double pe_max = MaxImpulse(pnu, a, b, c);

	//double low = (2 * pnu * pe * (1 - cos(theta)) - me * me) / (4 * pnu * pe_max - me * me);
	////return eps2;
	return low;
}

double ZeroIntLim(double pnu, double q2)
{
	/*return eps2;*/
	return 0;
}

double SubIntLim(double pnu, double q2)
{
	return 1 - eps2;
}

//double MaxImpulse(double pnu, double a, double b, double c)
//{
//	double pe_max = (-b + sqrt(Discriminant(a, b, c))) / (2 * a);
//	return pe_max;
//}

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


