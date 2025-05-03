#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

const double d=4, x0=4, dx=1, w=1.9;
const int N=31;


double ro(int, int);

double roP(int, int, double**);

double nad_relaksacja(int, int, double**);

double S(double**);

void zapiszS(int, double, std::ofstream *);

void zapiszU(double**, std::ofstream *);

int main(){


    double** u = new double*[2*N+1];
    for (int i = 0; i < 2*N+1; i++) {
        u[i] = new double[2*N+1]();
    }

    double** v = new double*[2*N+1];
    for (int i = 0; i < 2*N+1; i++) {
        v[i] = new double[2*N+1]();
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
    std::string folder = "/Users/michau/Documents/MOFIT1/results_p02/results_2/";

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
                v[i][j] = nad_relaksacja(i,j,u);
                ro_m[i][j] = roP(i,j,v);
                d_m[i][j] = roP(i,j,v) - ro(i,j);
            } 
        }
        u = v;

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


    for (int i = 0; i < 2*N+1; i++) {
        delete[] v[i];
        delete[] ro_m[i];
        delete[] d_m[i];
    }
    delete[] v;
    delete[] ro_m;
    delete[] d_m;
    

    return 0;
}

double ro(int x, int y){
    return exp(-(pow((x-x0-N),2)+pow(y-N, 2))/(d*d)) - exp(-(pow((x+x0-N),2)+pow(y-N, 2))/(d*d));
}

double roP(int i, int j, double** u){
    return -(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] -4*u[i][j])/pow(dx,2);
}

double nad_relaksacja(int i, int j, double** u){
    return (1-w)*u[i][j] + w*((u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] + ro(i,j)*pow(dx,2))/4);
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


