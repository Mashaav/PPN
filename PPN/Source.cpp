/*
Zaimplementuj velocity Verlet dla wody. oprócz oddziaływań niewiążących
liczone mają być wiązania kowalencyjne oraz kąty między atomami. 
Za PBC przyjmujemy wartość minimalną i maksymalną danej składowej w pliku PDB w każdym kierunku ±0.5 Å. x y z

*/
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <iterator>
using namespace std;
#define N 30 // liczba atomów

//stałe
float MasO = 16, MasH = 1, MasH2O = 18; //h2o

/*
Vbond = sum Kb*(b - b0) ^ 2
Fbond = sum Kb* 2(b - b0)

Kb - stretching force constant, -> 450
b - bond length,
b0 - equilibrium distance parameter -> 0.9572

	
VbondAngle = sum ko*(o - o0) ^ 2
Fbond = sum Ko*2(o - o0)

Ko - bending angle constant, -> 55
o - tetha - bond angle, 
o0 - equilibrium value parameter -> 105°(equilibrium)

*/

//Fbond
double Kb = 450, //stała
b0 = 0.9572, //equilibrium

//FbondAngle
Ktheta = 55, //stała
theta0 = 105; // equilibrium

float sigmaO = 1.7682,
	  sigmaH = 0.2245,
	  epsO = 0.1521,
	  epsH = 0.0460;

void writePDB(double r[][3], char nazwaPDB[])
{
	int i;
	FILE *plik_pdb;
	plik_pdb = fopen(nazwaPDB, "a");
	if (plik_pdb == NULL)
	{
		perror("NIE MOGE OTWORZYC PLIKU PDB! ABY DOPISAC DO NIEGO DANE");
	}
	for (i = 0; i<N; i++)
	{
		if (i % 3 == 0)//tlen
			fprintf(plik_pdb, "ATOM   %4.d  OH2  TIP3W%7.1d  %8.3lf%8.3lf%8.3lf  1.00  0.00      WT1 O\n", i + 1, i, r[i][0], r[i][1], r[i][2]);
		if (i % 3 == 1)//wodór 1
			fprintf(plik_pdb, "ATOM   %4.d  H1   TIP3W%7.1d  %8.3lf%8.3lf%8.3lf  1.00  0.00      WT1 H\n", i + 1, i-1, r[i][0], r[i][1], r[i][2]);
		if (i % 3 == 2)//wodór 2
			fprintf(plik_pdb, "ATOM   %4.d  H2   TIP3W%7.1d  %8.3lf%8.3lf%8.3lf  1.00  0.00      WT1 H\n", i + 1, i-2 , r[i][0], r[i][1], r[i][2]);

	}
	fprintf(plik_pdb, "END\n");
	fclose(plik_pdb);
}

double Ekin(double v[N][3], float Mas) //funkcja do mierzenia średniej Energii kin
{
	int i;
	double Ek = 0; //energia kinetyczna 
	for (i = 0; i<N; i++) {
		Ek += Mas*0.5*(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
	}
	return Ek;
}

double GetAngleABC(double a[3], double b[3], double c[3])
{
	//odległości a do b oraz b do c
	double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
	double bc[3] = { c[0] - b[0], c[1] - b[1], c[2] - b[2] };

	//dVec = sqrt(x^2 + y^2 + z^2)
	double abVec = sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]);
	double abNorm[3] = { ab[0] / abVec, ab[1] / abVec, ab[2] / abVec };

	double bcVec = sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]);
	double bcNorm[3] = { bc[0] / bcVec, bc[1] / bcVec, bc[2] / bcVec };

	double res = abNorm[0] * bcNorm[0] + abNorm[1] * bcNorm[1] + abNorm[2] * bcNorm[2];

	return (3.141592653589793 - acos(res)) * 180.0 / 3.141592653589793; //stopnie. Radiany -> zakomentować przed mnożeniem
	}

double Distance(double rO[3], double rH[3])
{
	// d = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
	return sqrt(pow(rO[0] - rH[0], 2) + pow(rO[1] - rH[1], 2) + pow(rO[2] - rH[2], 2));
}

void Acc(const double r[N][3], double a[N][3]) // obliczanie przyspieszeń / sił działających na atomy
{
	int i, j, k;
	double dx, dy, dz, dr, FnonBounded, F, Fbounded[N/3] = {0};
	float Mas = 0, eps = 0, sigma = 0;
	
	//podział tablicy głównej na O, H1, H2
	double rO[N / 3][3] = { 0 };
	double rH1[N / 3][3] = { 0 };
	double rH2[N / 3][3] = { 0 };

	//tablice sił
	double FBond[N / 3] = { 0 }, //liczona pomiędzy 'O-H', suma dla obu atomów wodoru
		FBondAngle[N / 3] = { 0 };

	//tablice zmiennych, dla każdego 'O' wyznaczane na podstawie pozycji atomów: b -> długość wiązania'O-H', 2. theta -> kąt pomiędzy wiązaniami 'H-O-H'
	double b[(2*N / 3)] = { 0 },//bond length -> dla N 30, ilość wiązań to 20
		theta[N / 3] = { 0 }; // bond angle -> dla N 30, ilość kątów to 10

	//ustawienie wartości z tablicy głównej zawierającej wszystkie pozycje, na poszczególne O, H1 i H2
	for (i = 0; i < N; i++)
	{
		if (i % 3 == 0)
		{
			rO[i / 3][0] = r[i][0];
			rO[i / 3][1] = r[i][1];
			rO[i / 3][2] = r[i][2];
			rH1[i / 3][0] = r[i + 1][0];
			rH1[i / 3][1] = r[i + 1][1];
			rH1[i / 3][2] = r[i + 1][2];
			rH2[i / 3][0] = r[i + 2][0];
			rH2[i / 3][1] = r[i + 2][1];
			rH2[i / 3][2] = r[i + 2][2];
		}
	}

	//1. Obliczenie b oraz theta
	for (i = 0; i < N; i++)
	{
		if (i % 3 == 0) // dla każdego tlenu
		{
			//bond O-H1
			b[(2 * i / 3)] = Distance(
				new double[3]{ rO[i / 3][0], rO[i / 3][1], rO[i / 3][2] },
				new double[3]{ rH1[i / 3][0], rH1[i / 3][1], rH1[i / 3][2] }
			);

			//bond O-H2
			b[(2 * i / 3) + 1] = Distance(
				new double[3]{ rO[i / 3][0], rO[i / 3][1], rO[i / 3][2] },
				new double[3]{ rH2[i / 3][0], rH2[i / 3][1], rH2[i / 3][2] }
			);

			//kąt H-O-H w stopniach, -> w razie potrzeby można zmień w funkcji na radiany
			theta[i/3] = GetAngleABC(
				new double[3]{ rH1[i / 3][0], rH1[i / 3][1], rH1[i / 3][2] },	//xyz rH1
				new double[3]{ rO[i / 3][0], rO[i / 3][1], rO[i / 3][2] },		//xyz rO
				new double[3]{ rH2[i / 3][0], rH2[i / 3][1], rH2[i / 3][2] }	//xyz rH2
			);
		}
	}

	//2. Obliczanie sił
	for (i = 0; i < N; i++)
	{
		//dla każdego tlenu
		if(i % 3 == 0)
		{
			//suma sił dla obu łączeń
			FBond[i/3] = (Kb * 2 * (b[2*i/3] - b0)) + (Kb * 2 * (b[(2*i/3) + 1] - b0));
			
			//siła kąta
			FBondAngle[i/3] = Ktheta * 2 * (theta[i/3] - theta0);
			
			//siła wiążąca
			Fbounded[i/3] = FBond[i/3] + FBondAngle[i/3];
		}
	}

	//3.zerowanie przyspieszeń
	for (i = 0; i < N; i++)
		for (j = 0; j < 3; j++) 
		{
			a[i][j] = 0; 
		}

	//4. Obliczanie przyspieszeń
	for (i = 0; i < N; i++) 
	{
		for (j = i + 1; j < N; j++) 
		{
			k = 0;
			eps = 0;
			sigma = 0;
			F = 0;
			dx = 0;
			dy = 0;
			dz = 0;
			dx = r[i][0] - r[j][0];
			dy = r[i][1] - r[j][1];
			dz = r[i][2] - r[j][2];
			dr = sqrt(dx*dx + dy*dy + dz*dz);
			

			//w zależności od tego z jakim aktualny atom (i) oddziaływuje (j), ustawiam odpowienio epsilon i sigma, oraz jaka siła jest liczona. W przypadku, gdy nie są to 'sąsiedzi cząsteczkowi' to nie uwzględniam sił wiążących
			if (i % 3 == 0)//i - tlen
			{
				if (j % 3 == 0)//j - tlen
				{
					eps = epsO;
					sigma = sigmaO;
				}
				else// j - wodór
				{
					eps = sqrt(epsO * epsH);
					sigma = 0.5 * sigmaO + 0.5 * sigmaH;
				}
				
				Mas = 16;
				FnonBounded = eps*(48 * pow(dr / sigma, -13.0) - 24 * pow(dr / sigma, -7.0));
				F = FnonBounded;
				
				k = i / 3;
				if (j == i + 1 || j == i + 2)
				{
					F = Fbounded[k] + FnonBounded;
				}
			}
			else if (i % 3 == 1)//i - wodór 1
			{
				if (j % 3 != 0) //j - wodór
				{
					eps = epsH;
					sigma = sigmaH;
				}
				else // j - tlen
				{
					eps = sqrt(epsO * epsH);
					sigma = 0.5 * sigmaO + 0.5 * sigmaH;
				}

				Mas = 1;
				FnonBounded = eps*(48 * pow(dr / sigma, -13.0) - 24 * pow(dr / sigma, -7.0));
				F = FnonBounded;
				
				k = (i - 1) / 3;
				if (/*j == i - 1 ||*/ j == i + 1) // warunek zakomentowany - nie sprawdzany przez budowę pętli
				{
					F = FnonBounded + Fbounded[k];
				}
			}
			else if (i % 3 == 2)//i - wodór 2
			{
				if (j % 3 != 0)//j - wodór
				{
					eps = epsH;
					sigma = sigmaH;
				}
				else// j - tlen
				{
					eps = sqrt(epsO * epsH);
					sigma = 0.5 * sigmaO + 0.5 * sigmaH;
				}

				Mas = 1;
				FnonBounded = eps*(48 * pow(dr / sigma, -13.0) - 24 * pow(dr / sigma, -7.0));
				F = FnonBounded;
				
				//k = (i - 2) / 3;
				//if (j == i - 1 || j == i - 2)//nie sprawdzam sił z poprzednikami - sąsiadami cząsteczkowymi - (budowa pętli) więc to nie zostanie sprawdzone
				//{
				//	F = FnonBounded + Fbounded[k];
				//}
			}
			
			a[i][0] += F*dx / (Mas * dr);
			a[j][0] -= F*dx / (Mas * dr);
			a[i][1] += F*dy / (Mas * dr);
			a[j][1] -= F*dy / (Mas * dr);
			a[i][2] += F*dz / (Mas * dr);
			a[j][2] -= F*dz / (Mas * dr);

		}
	}
}


int main() {
	
	// polozenia N-indeksuje atomy, x->0,y->1,z->2
	double rO[N/3][3] = {
		{ -4.256, 1.799, 2.660 },
		{ 3.668,   1.864, -1.983 },
		{ -5.345, -0.364, -5.026 },
		{ 3.420,   0.559, -4.565 },
		{ 3.920, -1.680, -5.951 },
		{ 1.886,   5.932, -1.648 },
		{ 2.929,   5.831,   3.222 },
		{ 1.817,   3.536,   1.858 },
		{ 2.072,   4.922,   5.765 },
		{ 3.276, -6.405, -4.916 }
	};
	double rH1[N/3][3] = {
		{ -4.896,   1.245,   3.029 },
		{ 3.992,   1.826 ,-1.037 },
		{ -5.063, -1.306, -4.988 },
		{ 3.755,   1.356 ,-4.975 },
		{ 4.132, -2.276 ,-5.215 },
		{ 1.113,   5.433 ,-1.526 },
		{ 2.703,   5.053 ,  2.665 },
		{ 1.177,   3.708 ,  1.179 },
		{ 1.222,   4.516 ,  6.024 },
		{ 3.361, -6.086 ,-3.953 }
	};
	double rH2[N/3][3] = {
		{ -4.842,   2.200,   1.993 },
		{ 2.742,   2.225, -1.894 },
		{ -6.228, -0.446, -5.486 },
		{ 3.242,   0.844, -3.661 },
		{ 3.651, -0.833, -5.476 },
		{ 2.336,   5.512, -2.433 },
		{ 2.466,   5.721,   4.049 },
		{ 1.383,   3.016,   2.565 },
		{ 2.506,   5.149,   6.591 },
		{ 2.737, -5.722, -5.302 }
	};
	double r[N][3] = { 0 };
	
	//wypełnienie tablicy głównej z danych cząstkowych (osobno tleny, osobno wodory)
	for (int i = 0; i < N; i++)
	{
		if (i == 0 || i % 3 == 0)
		{
			r[i][0] = rO[i / 3][0];
			r[i][1] = rO[i / 3][1];
			r[i][2] = rO[i / 3][2];
			r[i + 1][0] = rH1[i / 3][0];
			r[i + 1][1] = rH1[i / 3][1];
			r[i + 1][2] = rH1[i / 3][2];
			r[i + 2][0] = rH2[i / 3][0];
			r[i + 2][1] = rH2[i / 3][1];
			r[i + 2][2] = rH2[i / 3][2];
		}
	}
	
	// predkosci N-indeksuje atomy, x->0,y->1,z->2
	double vO[N/3][3] = {
		{ 2.283, -0.112, -1.464 },
		{ 11.920, 2.644, 1.029 },
		{ -5.456, 4.579, 6.238 },
		{ 2.340, -1.936, 0.826 },
		{ -2.165, 2.500, -5.092 },
		{ -3.376 , 1.966, -4.115 },
		{ -1.358, -2.464, 4.517 },
		{ -9.440, 3.746, 2.581 },
		{ -4.861, -1.532, -1.235 },
		{ 7.517, -0.267, 1.367 }
	};
	double vH1[N/3][3] = {
		{ -13.660, -10.622, 17.822 },
		{ -16.554, 16.892, -22.647 },
		{ -6.530, -0.957, 6.291 },
		{ 12.040, 8.431, 0.436 },
		{ -0.835, 3.812, 16.868 },
		{ -8.404, -12.044, -6.300 },
		{ -0.263, -25.157, 14.399 },
		{ -23.196, -19.135, -4.375 },
		{ 6.942, 17.259, 0.043 },
		{ -2.662, 17.386, -26.157 }
	};
	double vH2[N/3][3] = {
		{ -5.968, 6.647 -0.448 },
		{ 21.376, 3.860, 0.046 },
		{ 19.546, 24.586, -2.407 },
		{ -8.141, 10.003, 11.059 },
		{ 20.312, -0.962, -9.369 },
		{ -27.244, 2.290, 7.731 },
		{ 6.389, 18.205, -4.177 },
		{ 6.959, -3.642, -11.562 },
		{ -9.771, 12.903, -7.065 },
		{ 6.155, 14.977, -4.926 }
	};
	double v[N][3] = { 0 };
	
	//wypełenienie v z danych cząstkowych
	for (int i = 0; i < N; i++)
	{
		if (i == 0 || i % 3 == 0)
		{
			v[i][0] = vO[i / 3][0];
			v[i][1] = vO[i / 3][1];
			v[i][2] = vO[i / 3][2];
			v[i + 1][0] = vH1[i / 3][0];
			v[i + 1][1] = vH1[i / 3][1];
			v[i + 1][2] = vH1[i / 3][2];
			v[i + 2][0] = vH2[i / 3][0];
			v[i + 2][1] = vH2[i / 3][1];
			v[i + 2][2] = vH2[i / 3][2];
		}
	}

	// przyspieszenia N-indeksuje atomy, x->0,y->1,z->2
	double a[N][3] = { 0 };
	
	//zmienne pomocnicze
	int i, j, k;

	double vMax = 5.775;   // A/ps
	double dt = 0.001;  // timestep krok czasowy w ps
	int TIME = 10000; 

	//Periodic Boundary Conditions - wczytywane najwyższe wartości z pliku dla x y z +- 0.5
	float L = 14; // orginał(argon)
	float Lx = 6.728, Ly = 6.905, Lz = 7.091; // verlet dla wody

	char nPDB[20] = "H2O_MD.pdb"; //nazwa pliku PDB
	double kb = 1.38064852;  //Stała boltzmanna w J/K
	double Ek; //energia kinetyczna

	FILE *plik_pdb;
	FILE *data;

	srand(time(NULL));

	plik_pdb = fopen(nPDB, "w");
	if (plik_pdb == NULL)
	{
		perror("NIE MOGE OTWORZYC PLIKU PDB!");
		return 1;
	}
	fclose(plik_pdb);

	data = fopen("data.dat", "w");
	if (plik_pdb == NULL)
	{
		perror("NIE MOGE OTWORZYC PLIKU PDB!");
		return 1;
	}

	//writePDB(r, nPDB);

	for (k = 1; k <= TIME; k++) 
	{
		//odświeżam przyspieszenia
		Acc(r, a);
		
		//odświeżam pozycję i na wpół prędkości
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < 3; j++) 
			{
				//pozycja
				r[i][j] += v[i][j] * dt + 0.5 * a[i][j] * dt * dt; 
				
				// PBC - periodic boundary conditions
				
				//przerzucanie cząsteczki
				if (i % 3 == 0)//tlen
				{

					if (j == 0 && r[i][j] < -Lx)//>x
					{
						try
						{
							r[i][j] += 2 * Lx;
							r[i + 1][j] += 2 * Lx;
							r[i + 2][j] += 2 * Lx;
						}
						catch(exception e)
						{ }
					}
					if (j == 1 && r[i][j] < -Ly)//>y
					{
						try
						{
							r[i][j] += 2 * Ly;
							r[i + 1][j] += 2 * Ly;
							r[i + 2][j] += 2 * Ly;
						}
						catch (exception e)
						{
						}
					}
					if (j == 2 && r[i][j] < -Lz)//>z
					{
						try
						{
							r[i][j] += 2 * Lz;
							r[i + 1][j] += 2 * Lz;
							r[i + 2][j] += 2 * Lz;
						}
						catch (exception e)
						{
						}
					}
				
					if (j == 0 && r[i][j] > Lx)//>x
					{
						try
						{
							r[i][j] -= 2 * Lx;
							r[i+1][j] -= 2 * Lx;
							r[i+2][j] -= 2 * Lx;
						}
						catch (exception e)
						{
						}
					}
					if (j == 1 && r[i][j] > Ly)//>y
					{
						try
						{
							r[i][j] -= 2 * Ly;
							r[i + 1][j] -= 2 * Ly;
							r[i + 2][j] -= 2 * Ly;
						}
						catch (exception e)
						{
						}
					}
					if (j == 2 && r[i][j] > Lz)//>z
					{
						try
						{
							r[i][j] -= 2 * Lz;
							r[i + 1][j] -= 2 * Lz;
							r[i + 2][j] -= 2 * Lz;
						}
						catch (exception e)
						{
						}
					}
				}

				if (i % 3 == 1)//wodór1
				{

					if (j == 0 && r[i][j] < -Lx)//<x
					{
						try
						{
							r[i - 1][j] += 2 * Lx;
							r[i][j] += 2 * Lx;
							r[i + 1][j] += 2 * Lx;
						}
						catch (exception e)
						{
						}
					}
					if (j == 1 && r[i][j] < -Ly)//<y
					{
						try
						{
							r[i - 1][j] += 2 * Ly;
							r[i][j] += 2 * Ly;
							r[i + 1][j] += 2 * Ly;
						}
						catch (exception e)
						{
						}
					}
					if (j == 2 && r[i][j] < -Lz)//<z
					{
						try
						{
							r[i - 1][j] += 2 * Lz;
							r[i][j] += 2 * Lz;
							r[i + 1][j] += 2 * Lz;
						}
						catch (exception e)
						{
						}
					}

					if (j == 0 && r[i][j] > Lx)//>x
					{
						try
						{
							r[i - 1][j] -= 2 * Lx;
							r[i][j] -= 2 * Lx;
							r[i + 1][j] -= 2 * Lx;
						}
						catch (exception e)
						{
						}
					}
					if (j == 1 && r[i][j] > Ly)//>y
					{
						try
						{
							r[i - 1][j] -= 2 * Ly;
							r[i][j] -= 2 * Ly;
							r[i + 1][j] -= 2 * Ly;
						}
						catch (exception e)
						{
						}
					}
					if (j == 2 && r[i][j] > Lz)//>z
					{
						try
						{
							r[i - 1][j] -= 2 * Lz;
							r[i][j] -= 2 * Lz;
							r[i + 1][j] -= 2 * Lz;
						}
						catch (exception e)
						{
						}
					}
				}

				if (i % 3 == 0)//wodór2
				{

					if (j == 0 && r[i][j] < -Lx)//<x
					{
						try
						{
							r[i - 2][j] += 2 * Lx;
							r[i - 1][j] += 2 * Lx;
							r[i][j] += 2 * Lx;
						}
						catch (exception e)
						{
						}
					}
					if (j == 1 && r[i][j] < -Ly)//<y
					{
						try
						{
							r[i - 2][j] += 2 * Ly;
							r[i - 1][j] += 2 * Ly;
							r[i][j] += 2 * Ly;
						}
						catch (exception e)
						{
						}
					}
					if (j == 2 && r[i][j] < -Lz)//<z
					{
						try
						{
							r[i - 2][j] += 2 * Lz;
							r[i - 1][j] += 2 * Lz;
							r[i][j] += 2 * Lz;
						}
						catch (exception e)
						{
						}
					}

					if (j == 0 && r[i][j] > Lx)//>x
					{
						try
						{
							r[i - 2][j] -= 2 * Lx;
							r[i - 1][j] -= 2 * Lx;
							r[i][j] -= 2 * Lx;
						}
						catch (exception e)
						{
						}
					}
					if (j == 1 && r[i][j] > Ly)//>y
					{
						try
						{
							r[i - 2][j] -= 2 * Ly;
							r[i - 1][j] -= 2 * Ly;
							r[i][j] -= 2 * Ly;
						}
						catch (exception e)
						{
						}
					}
					if (j == 2 && r[i][j] > Lz)//>z
					{
						try
						{
							r[i - 2][j] -= 2 * Lz;
							r[i - 1][j] -= 2 * Lz;
							r[i][j] -= 2 * Lz;
						}
						catch (exception e)
						{
						}
					}
				}
				
				//0.5 * prędkość
				v[i][j] += 0.5 * a[i][j] * dt; 
			}
		}

		// odświeżam przyspieszenia
		Acc(r, a); 
		
		// drugie 0.5 * prędkość
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < 3; j++) 
			{
				v[i][j] += 0.5 * a[i][j] * dt; 
			}
		}

		if (k % 10 == 0) 
		{
			writePDB(r, nPDB);
			Ek = Ekin(v, MasH2O);
			fprintf(data, "%d\t%lf\t%lf\n", k, Ek, 2.0*1.66054*Ek / (3 * N*kb)); //2.0*1.66054*Ek/(3*N*kb) bierze się z przeliczenia EK na J (1.66054) oraz podstawienia kb w J/K (1.38064852) 10^34 w liczniku i mianowniku skraca sie
		}
	}

	fclose(data);
	printf("GOTOWE!\nPlik znajdują się w folderze projektu. Plik wynikowy to H2O_MD.pdb\n");
	system("Pause");
	return 0;
}