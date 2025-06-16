#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>

const double mu = 1.0, ro = 1.0, dz = 0.01;
const double Q[5] = {-1,-10,-100,-200,-400};
const int i1 = -100, i2 = 150, jj1 = -40, jj2 = 40;
const int i12 = i2 - i1, j12 = jj2 - jj1;
const int ik = 5, jk = 10;
const double Y1 = jj1 * dz, Y2 = jj2 * dz;

double phi0(int, int, int);

void ksi0(double**, double**);

double relaksacjaPhi(int, int, double**, double);

double relaksacjaKsi(int, int, double**, double**);

void zapiszU(double**, std::ofstream&);

int main() {
    std::string folder = "/Users/michau/Documents/MOFIT1/results_p03/results_2/";

    

    for(int iii=0; iii<5; iii++){

        std::ofstream outPhi(folder + "Phi_"+ std::to_string(iii) + ".txt");
        std::ofstream outKsi(folder + "Ksi_"+ std::to_string(iii) + ".txt");

        double** phi = new double*[i12 + 1];
        double** ksi = new double*[i12 + 1];
        for (int i = 0; i <= i12; i++) {
            phi[i] = new double[j12 + 1]();
            ksi[i] = new double[j12 + 1]();
        }


        for (int i = 0; i <= i12; i++) {
            for (int j = 0; j <= j12; j++) {
                phi[i][j] = phi0(i, j, iii);
                ksi[i][j] = 0.0;
            }
        }
        ksi0(ksi, phi);

        for (int iter = 0; iter < 10000; iter++) {
            for (int i = 1; i < i12; i++) {
                for (int j = 1; j < j12; j++) {

                    bool zastawka = !((i >= -ik - i1 && i <= ik - i1) && j <= jk - jj1);
                    if (zastawka) {
                        phi[i][j] = relaksacjaPhi(i, j, phi, ksi[i][j]);
                        ksi[i][j] = relaksacjaKsi(i, j, phi, ksi);
                    }
                }
            }
            ksi0(ksi, phi);
        }



        zapiszU(phi, outPhi);
        zapiszU(ksi, outKsi);

        outPhi.close();
        outKsi.close();

        for (int i = 0; i <= i12; i++) {
            delete[] phi[i];
            delete[] ksi[i];
        }
        delete[] phi;
        delete[] ksi;
    }

    return 0;
}

double phi0(int i, int j, int iii) {
    double y = (j + jj1) * dz;

    if (j == 0 || j == j12 || i == 0 || i == i12)
        return 0.5 * Q[iii] / mu * (pow(y, 3) / 3.0 - pow(y, 2) / 2.0 * (Y1 + Y2) + Y1*Y2*y);

    if ((i == ik - i1 && j <= jk - jj1) || (i == -ik - i1 && j <= jk - jj1))
        return 0.5 * Q[iii] / mu * (pow(Y1, 3) / 3.0 - pow(Y1, 2) / 2.0 * (Y1 + Y2) + Y1*Y2*Y1);

    if ((j == jk - jj1) && (i > -ik - i1 && i < ik - i1))
        return 0.5 * Q[iii] / mu * (pow(Y1, 3) / 3.0 - pow(Y1, 2) / 2.0 * (Y1 + Y2) + Y1*Y2*Y1);

    return 0.0;
}

void ksi0(double** ksi, double** phi) {
    for (int i = 0; i <= i12; i++) {
        ksi[i][j12] = 2 * (phi[i][j12 - 1] - phi[i][j12]) / (dz * dz);
        ksi[i][0] = 2 * (phi[i][1] - phi[i][0]) / (dz * dz);
    }
    for (int j = 0; j <= jk - jj1; j++) {
        ksi[ik - i1][j] = 2 * (phi[ik - i1 + 1][j] - phi[ik - i1][j]) / (dz * dz);
        ksi[-ik - i1][j] = 2 * (phi[-ik - i1 - 1][j] - phi[-ik - i1][j]) / (dz * dz);
    }
    for (int i = -ik - i1 + 1; i < ik - i1; i++) {
        ksi[i][jk - jj1] = 2 * (phi[i][jk - jj1 + 1] - phi[i][jk - jj1]) / (dz * dz);
    }
    ksi[-ik - i1][jk - jj1] = 0.5 * (ksi[-ik - i1][jk - jj1 - 1] + ksi[-ik - i1 + 1][jk - jj1]);
    ksi[ik - i1][jk - jj1] = 0.5 * (ksi[ik - i1][jk - jj1 - 1] + ksi[ik - i1 - 1][jk - jj1]);
}

double relaksacjaPhi(int i, int j, double** phi, double ksi) {
    return (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1] - ksi * dz * dz) / 4.0;
}

double relaksacjaKsi(int i, int j, double** phi, double** ksi) {
    return (ksi[i + 1][j] + ksi[i - 1][j] + ksi[i][j + 1] + ksi[i][j - 1]) / 4.0
        - ((phi[i][j + 1] - phi[i][j - 1]) * (ksi[i + 1][j] - ksi[i - 1][j])
         - (phi[i + 1][j] - phi[i - 1][j]) * (ksi[i][j + 1] - ksi[i][j - 1])) / 16.0;
}

void zapiszU(double** u, std::ofstream& out) {
    for (int i = 0; i <= i12; i++) {
        for (int j = 0; j <= j12; j++) {
            out << u[i][j] << " ";
        }
        out << "\n";
    }
}

