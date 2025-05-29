#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

const double Q=-1., mu=1, ro=1, dz=0.01;
const int i1=-100, i2=150, j1=-40, j2=40;
const double X1=i1*dz, X2=i2*dz, Y1=-j1*dz, Y2=j2*dz;

void kopiuj(double**, double**, int);

double phi0(int, int);

double ksi0(int, int);

double relaksacjaPhi(int, int, double**, double);

double relaksacjaKsi(int, int, double**, double**);

double S(double**);

void zapiszS(int, double, std::ofstream *);

void zapiszU(double**, std::ofstream *);

int main(){


    double** ksi = new double*[i2-i1+1];
    for (int i = 0; i < i2-i1+1; i++) {
        ksi[i] = new double[j2-j1+1]();
    }

    double** phi = new double*[i2-i1+1];
    for (int i = 0; i < i2-i1+1; i++) {
        phi[i] = new double[j2-j1+1]();
    }

    // double** ro_m = new double*[2*N+1];
    // for (int i = 0; i < 2*N+1; i++) {
    //     ro_m[i] = new double[2*N+1]();
    // }

    // double** d_m = new double*[2*N+1];
    // for (int i = 0; i < 2*N+1; i++) {
    //     d_m[i] = new double[2*N+1]();
    // }

    // sciezka folderu do zapisania
    std::string folder = "/Users/michau/Documents/MOFIT1/results_p02/results_1/";

    std::string path100 = folder + "u_at_100.txt"; 
    std::ofstream outW100(path100);

    std::string path500 = folder + "u_at_500.txt"; 
    std::ofstream outW500(path500);

    std::string pathS = folder + "wart_S.txt"; 
    std::ofstream outWS(pathS);

    std::string pathRoP100 = folder + "roP_at_100.txt"; 
    std::ofstream outWro100(pathRoP100);

    std::string pathRoP500 = folder + "roP_at_500.txt"; 
    std::ofstream outWro500(pathRoP500);

    std::string pathD100 = folder + "d_at_100.txt"; 
    std::ofstream outWd100(pathD100);

    std::string pathD500 = folder + "d_at_500.txt"; 
    std::ofstream outWd500(pathD500);

    
    for(int iter=0; iter<500; iter++){
        for(int i=1; i<2*N; i++){
            for(int j=1; j<2*N; j++){
                u[i][j] = relaksacja(i,j,u);
            } 
        }
        for(int i=1; i<2*N; i++){
            for(int j=1; j<2*N; j++){
                ro_m[i][j] = roP(i,j,u);
                d_m[i][j] = roP(i,j,u) - ro(i,j);
            } 
        }
        // kopiuj(u, v, 2*N+1);

        if(iter==99){
            zapiszU(u, &outW100);
            outW100.close();
            
            zapiszU(ro_m, &outWro100);
            outWro100.close();

            zapiszU(d_m, &outWd100);
            outWd100.close();

        }
        if(iter==499){
            zapiszU(u, &outW500);
            outW500.close();

            zapiszU(ro_m, &outWro500);
            outWro500.close();

            zapiszU(d_m, &outWd500);
            outWd500.close();

        }
        
        zapiszS(iter, S(u), &outWS);
        std::cout<< S(u) << std:: endl;
    }

    outWS.close();


    for (int i = 0; i < j2-j1+1; i++) {
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
    if(i!=0 && i!=i1-i2 && j!=0 && j!=j1-j2)
        return 0;
    else{
        double x = i*dz;
        double y = j*dz;
        return 0.5*Q/mu*(pow(y,3)/3-pow(y,2)/2*(Y1+Y2)+Y1*Y2*y);
    }
}

double ksi0(int i, int j){
    if(i!=0 && i!=i1-i2 && j!=0 && j!=j1-j2)
        return 0;
    else{
        double x = i*dz;
        double y = j*dz;
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

double S(double ** u){
    double suma=0;
    for(int i=1; i<2*N; i++){
        for(int j=1; j<2*N; j++){
            suma += u[i][j]*(
                0.5*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j])
                + ro(i,j)*pow(dx,2));
        }
    }
    return -suma;
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


