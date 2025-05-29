#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

const double d=4, x0=4, dx=1, w=1.9;
const int N=31;

void kopiuj(double**, double**, int);

double ro(int, int);

double roP(int, int, double**);

double S(double**, double**);

double S_loc (double**, int, int, double**);

double S_d (double**, double, int, int, double, double**);

double update_point(int, int, double**, double**);

void zapiszS(int, double, std::ofstream *);

void zapiszU(double**, std::ofstream *);

int main(){


    double** u = new double*[2*N+1];
    for (int i = 0; i < 2*N+1; i++) {
        u[i] = new double[2*N+1]();
    }

    double** ro_m = new double*[2*N+1];
    for (int i = 0; i < 2*N+1; i++) {
        ro_m[i] = new double[2*N+1]();
    }

    double** d_m = new double*[2*N+1];
    for (int i = 0; i < 2*N+1; i++) {
        d_m[i] = new double[2*N+1]();
    }

    // sciezka folderu do zapisania
    std::string folder = "/Users/michau/Documents/MOFIT1/results_p02/results_3/";


    std::string path500 = folder + "u_at_500.txt"; 
    std::ofstream outW500(path500);

    std::string pathS = folder + "wart_S.txt"; 
    std::ofstream outWS(pathS);

    for(int i=0; i<2*N+1; i++)
        for(int j=0; j<2*N+1; j++)
            ro_m[i][j] = ro(i,j);
    
    
    double wartS=0;
    for(int iter=0; iter<500; iter++){
        for(int i=1; i<2*N; i++){
            for(int j=1; j<2*N; j++){
                u[i][j] = update_point(i,j,u, ro_m);
            } 
        }

        wartS = S(u, ro_m);
        zapiszS(iter, wartS, &outWS);
        std::cout<< wartS << std:: endl;
    }

    zapiszU(u, &outW500);
    outW500.close();

    outWS.close();


    for (int i = 0; i < 2*N+1; i++) {
        delete[] u[i];
        delete[] ro_m[i];
        delete[] d_m[i];
    }
    delete[] u;
    delete[] ro_m;
    delete[] d_m;
    

    return 0;
}

void kopiuj(double** dest, double** src, int rozmiar) {
    for (int i = 0; i < rozmiar; i++) {
        for (int j = 0; j < rozmiar; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

double ro(int x, int y){
    return exp(-(pow((x-x0-N),2)+pow(y-N, 2))/(d*d)) - exp(-(pow((x+x0-N),2)+pow(y-N, 2))/(d*d));
}

double roP(int i, int j, double** u){
    return -(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] -4*u[i][j])/pow(dx,2);
}

double S(double ** u, double** ro){
    double suma=0;
    for(int i=1; i<2*N; i++){
        for(int j=1; j<2*N; j++){
            suma += u[i][j]*(
                0.5*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j])
                + ro[i][j]*pow(dx,2));
        }
    }
    return -suma;
}

double S_loc (double** u, int ip, int jp, double** ro){
    double suma=0;
    // std::cout<<ip<<" "<<jp<<std::endl;
    if(ip!=1 && jp!=1 && ip!=2*N-1 && jp!=2*N-1)
        for(int i=ip-1; i<ip+2; i++){
            for(int j=jp-1; j<jp+2; j++){
                suma += u[i][j]*(
                    0.5*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j])
                    + ro[i][j]*pow(dx,2));
            }
        }
    else
        suma += u[ip][jp]*(
            0.5*(u[ip+1][jp] + u[ip-1][jp] + u[ip][jp+1] + u[ip][jp-1] - 4*u[ip][jp])
            + ro[ip][jp]*pow(dx,2));
    return -suma;
}

double S_d (double** u, double d, int ip, int jp, double S0, double** ro){
    double S = S_loc(u, ip, jp, ro);
    u[ip][jp] += d;

    double Sd = S_loc(u, ip, jp, ro);
    u[ip][jp] -= d;

    return S0 - S + Sd;
}

double update_point(int ip, int jp, double** u, double** ro){
    double S1 = S(u, ro);
    double S2 = S_d(u, 0.5, ip, jp, S1, ro);
    double S3 = S_d(u, 1.0, ip, jp, S1, ro);
    double d4 = 0.25* (3*S1 - 4*S2 + S3)/(S1 - 2*S2 + S3);
    double S4 = S_d(u, d4, ip, jp, S1, ro);

    std::vector<double> d = {0.0, 0.5, 1.0, d4};
    std::vector<double> SS = {S1, S2, S3, S4};
    int imin = min_element(SS.begin(), SS.end()) - SS.begin();
    return u[ip][jp] + d[imin];
}


void zapiszS(int i, double s, std::ofstream * out){
    *out << i << " " << s << std::endl;
}

void zapiszU(double** u, std::ofstream * out){
    for(int i=0; i<2*N+1; i++){
        for(int j=0; j<2*N+1; j++){
            *out << u[i][j] << " ";
        }
        *out << std::endl;
    }
}


