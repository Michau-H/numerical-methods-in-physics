#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

const double d=4, x0=4, dx=1;

const int N=31;

double u[2*N+1][2*N+1]= {};

double ro(int, int);

double relaksacja(int, int);

double S();

int main(){


    return 0;
}

double ro(int x, int y){
    return exp(-(pow((x-x0),2)+pow(y, 2))/(d*d)) - exp(-(pow((x+x0),2)+pow(y, 2))/(d*d));
}

double relakscaja(int i, int j){
    return (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] + ro(i,j)*pow(dx,2))/4;
}

double S(){
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
