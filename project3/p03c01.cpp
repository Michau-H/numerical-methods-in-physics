#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

const double Q=-1., mu=1, ro=1, dz=0.01;
const int i1=-100, i2=150, jj1=-40, jj2=40;
const int i12 = i2-i1, j12=jj2-jj1;
const double X1=i1*dz, X2=i2*dz, Y1=jj1*dz, Y2=jj2*dz;

void kopiuj(double**, double**, int);

double phi0(int, int);

double ksi0(int, int);

double relaksacjaPhi(int, int, double**, double);

double relaksacjaKsi(int, int, double**, double**);

double S(double**);

void zapiszS(int, double, std::ofstream *);

void zapiszU(double**, std::ofstream *);

int main(){


    double** ksi = new double*[i12+1];
    for (int i = 0; i < i12+1; i++) {
        ksi[i] = new double[j12+1]();
    }

    double** phi = new double*[i12+1];
    for (int i = 0; i < i12+1; i++) {
        phi[i] = new double[j12+1]();
    }

    // sciezka folderu do zapisania
    std::string folder = "/Users/michau/Documents/MOFIT1/results_p03/results_1/";

    std::string pathPhi = folder + "Phi.txt"; 
    std::ofstream outPhi(pathPhi);

    std::string pathKsi = folder + "Ksi.txt"; 
    std::ofstream outKsi(pathKsi);


    for(int i=0; i<i12+1; i++){
        for(int j=0; j<j12+1; j++){
            phi[i][j] = phi0(i,j);
            ksi[i][j] = ksi0(i,j);
        } 
    }
    
    
    for(int iter=0; iter<5000; iter++){
        for(int i=1; i<i12; i++){
            for(int j=1; j<j12; j++){
                phi[i][j] = relaksacjaPhi(i,j,phi, ksi[i][j]);
                ksi[i][j] = relaksacjaKsi(i,j,phi, ksi);
            } 
        }
    }
    zapiszU(phi, &outPhi);
    outPhi.close();
        

    zapiszU(ksi, &outKsi);
    outKsi.close();


    for (int i = 0; i < i12+1; i++) {
        delete[] phi[i];
        delete[] ksi[i];
        // delete[] d_m[i];
    }
    delete[] phi;
    delete[] ksi;
    // delete[] d_m;
    

    return 0;
}

void kopiuj(double** dest, double** src, int rozmiar) {
    for (int i = 0; i < rozmiar; i++) {
        for (int j = 0; j < rozmiar; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

double phi0(int i, int j){
    if(i!=0 && i!=i12 && j!=0 && j!=j12)
        return 0;
    else{
        double y = (j+jj1)*dz;
        return 0.5*Q/mu*(pow(y,3)/3.-pow(y,2)/2.*(Y1+Y2)+Y1*Y2*y);
    }
}

double ksi0(int i, int j){
    if(i!=0 && i!=i12 && j!=0 && j!=j12)
        return 0;
    else{
        double y = (j+jj1)*dz;
        return 0.5*Q/mu*(2*y-Y1-Y2);
    }
}

double relaksacjaPhi(int i, int j, double** phi, double ksi){
    return (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - ksi*pow(dz,2))/4;
}

double relaksacjaKsi(int i, int j, double** phi, double** ksi){
    return (ksi[i+1][j] + ksi[i-1][j] + ksi[i][j+1] + ksi[i][j-1])/4
            -((phi[i][j+1]-phi[i][j-1])*(ksi[i+1][j]-ksi[i-1][j])
            -(phi[i+1][j]-phi[i-1][j])*(ksi[i][j+1]-ksi[i][j-1]))/16;
}


void zapiszU(double** u, std::ofstream * out){
    for(int i=0; i<i12+1; i++){
        for(int j=0; j<j12+1; j++){
            *out << u[i][j] << " ";
        }
        *out << std::endl;
    }
}


