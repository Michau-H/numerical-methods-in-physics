#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>

const double mu = 1.0, ro = 1.0, dz = 0.01;
const double Q = 0;
const int i1 = -100, i2 = 150, jj1 = -40, jj2 = 40;
const int i12 = i2 - i1, j12 = jj2 - jj1;
const int ik = 5, jk = 10;
const double Y1 = jj1 * dz, Y2 = jj2 * dz;

double phi0(int, int);

double relaksacjaPhi(int, int, double**);

void zapiszU(double**, std::ofstream&);

int main() {
    std::string folder = "/Users/michau/Documents/MOFIT1/results_p03/results_3/";

    

    for(int iii=0; iii<1; iii++){

        std::ofstream outPhi(folder + "Phi_"+ std::to_string(iii) + ".txt");

        double** phi = new double*[i12 + 1];
        for (int i = 0; i <= i12; i++) {
            phi[i] = new double[j12 + 1]();
        }


        for (int i = 0; i <= i12; i++) {
            for (int j = 0; j <= j12; j++) {
                phi[i][j] = phi0(i, j);
            }
        }

        for (int iter = 0; iter < 10000; iter++) {
            for (int i = 1; i < i12; i++) {
                for (int j = 1; j < j12; j++) {

                    bool zastawka = !((i >= -ik - i1 && i <= ik - i1) && j <= jk - jj1);
                    if (zastawka) {
                        phi[i][j] = relaksacjaPhi(i, j, phi);
                    }
                }
            }
        }



        zapiszU(phi, outPhi);

        outPhi.close();

        for (int i = 0; i <= i12; i++) {
            delete[] phi[i];
        }
        delete[] phi;
    }

    return 0;
}

double phi0(int i, int j) {
    double y = (j + jj1) * dz;
    double A=1;
    if (j == 0 || j == j12 || i == 0 || i == i12)
        return A*y;

    if ((i == ik - i1 && j <= jk - jj1) || (i == -ik - i1 && j <= jk - jj1))
        return A*Y1;

    if ((j == jk - jj1) && (i > -ik - i1 && i < ik - i1))
        return A*Y1;

    return 0.0;
}

double relaksacjaPhi(int i, int j, double** phi) {
    return (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]) / 4.0;
}

void zapiszU(double** u, std::ofstream& out) {
    for (int i = 0; i <= i12; i++) {
        for (int j = 0; j <= j12; j++) {
            out << u[i][j] << " ";
        }
        out << "\n";
    }
}

