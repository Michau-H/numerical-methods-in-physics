#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

const double l1=1, l2=1/sqrt(8);
const double m0 = 1.;

struct Stan {
    long double t;
    double x,v,e;
};


// funkcje
double phi(double);

double a(double);

Stan war_pocz();

Stan Euler(Stan, double);

Stan Verlet(Stan, double);

Stan RK4(Stan, double);

void zapisz(Stan, std::ofstream *);


// main
int main() {
    
    /*
    long double DeltaTE[5] = {0.001,0.0005,0.0001,0.00005,0.00001};
    long double DeltaTV[5] = {0.1,0.05,0.01,0.005,0.001};
    long double DeltaTR[5] = {0.1,0.05,0.01,0.005,0.001};
    std::string name[3] = {"Euler", "Verlet", "RK4"};
    for(int k=0; k<3; k++){
        for(int j=0; j<5; j++){

            long double dt;
            if(k==0){
                dt=DeltaTE[j];
            }
            else if(k==1){
                dt=DeltaTV[j];
            }
            else if(k==2){
                dt=DeltaTR[j];
            }

            std::string n = name[k];
            n += ("_" + std::to_string(j) + ".txt");
            std::string folder = "/Users/michau/Documents/MOFIT1/results_p01/results_1a/";
            std::string path = folder + n; 
            std::ofstream outW(path);
            
            Stan niestan = {dt,dt,dt,dt};
            zapisz(niestan, &outW);
            int rozm = int(100./dt)+1;
            Stan * test = new Stan[rozm];

            test[0] = war_pocz();
            zapisz(test[0], &outW);
            
            if(k==0){
                
                for(int i=0; i <rozm; i++){
                    test[i+1] = Euler(test[i], dt);
                    zapisz(test[i+1], &outW);
                }
            }
            else if(k==1){
                for(int i=0; i <rozm; i++){
                    test[i+1] = Verlet(test[i], dt);
                    zapisz(test[i+1], &outW);
                }
            }
            else if(k==2){
                for(int i=0; i <rozm; i++){
                    test[i+1] = RK4(test[i], dt);
                    zapisz(test[i+1], &outW);
                }
            }
            delete[] test;
            outW.close();
            std::cout << "zapisane " << n << std::endl;
        }
    }
    */

    ////////////////////////////////////////////////////////////////////////////////////////
    // do badania zbierzności
    long double DeltaTE[7] = {0.0005,0.00025,0.0002,0.0001,0.00005,0.00002,0.00001};
    long double DeltaTV[7] = {0.02,0.016,0.0125,0.01,0.005,0.002,0.001};
    long double DeltaTR[7] = {0.02,0.016,0.0125,0.01,0.005,0.002,0.001};
    std::string name[3] = {"Euler", "Verlet", "RK4"};
    for(int k=0; k<3; k++){
        for(int j=0; j<7; j++){

            long double dt;
            if(k==0){
                dt=DeltaTE[j];
            }
            else if(k==1){
                dt=DeltaTV[j];
            }
            else if(k==2){
                dt=DeltaTR[j];
            }

            std::string n = name[k];
            n += ("_" + std::to_string(j) + ".txt");
            std::string folder = "/Users/michau/Documents/MOFIT1/results_p01/results_1b/";
            std::string path = folder + n; 
            std::ofstream outW(path);
            
            Stan niestan = {dt,0,0,0};
            zapisz(niestan, &outW);
            int rozm = int(10./dt)+1;
            Stan * test = new Stan[rozm];

            test[0] = war_pocz();
            zapisz(test[0], &outW);
            
            if(k==0){
                for(int i=0; i <rozm; i++){
                    test[i+1] = Euler(test[i], dt);
                    zapisz(test[i+1], &outW);
                }
            }
            else if(k==1){
                for(int i=0; i <rozm; i++){
                    test[i+1] = Verlet(test[i], dt);
                    zapisz(test[i+1], &outW);
                }
            }
            else if(k==2){
                for(int i=0; i <rozm; i++){
                    test[i+1] = RK4(test[i], dt);
                    zapisz(test[i+1], &outW);
                }
            }
            delete[] test;
            outW.close();
            std::cout << "zapisane " << n << std::endl;
        }
    }


    return 0;
}

// definicje
double phi(double x){
    return -pow(M_E,(-x*x/l1*l1))-8*pow(M_E,(-pow((x-2),2)/(l2*l2)));
}

double a(double x){
    return -1/m0*(phi(x + 0.001) - phi(x - 0.001))/(2*0.001);
}

Stan war_pocz(){
    Stan s = {0,2.8,0,0};
    s.e = m0*(phi(s.x)+0.5*pow(s.v,2));
    return s;
}

Stan Euler(Stan s1, double dt){
    Stan s2 = {(s1.t/dt+1)*dt, s1.x+dt*s1.v, s1.v+dt*a(s1.x), 0};
    s2.e = m0*(phi(s2.x)+0.5*pow(s2.v,2));
    return s2;
}

Stan Verlet(Stan s1, double dt){
    Stan s2 = {(s1.t/dt+1)*dt, s1.x+dt*s1.v+0.5*a(s1.x)*dt*dt, 0, 0};
    s2.v = s1.v+0.5*dt*(a(s1.x)+a(s2.x));
    s2.e = m0*(phi(s2.x)+0.5*pow(s2.v,2));
    return s2;
}

Stan RK4(Stan s1, double dt){
    double x = s1.x, v = s1.v;
    double k1v = v, k1a=a(x);
    double k2v = v + 0.5*dt*k1a, k2a = a(x + 0.5*dt*k1v);
    double k3v = v + 0.5*dt*k2a, k3a = a(x + 0.5*dt*k2v);
    double k4v = v + dt*k3a, k4a = a(x + dt*k3v);
    Stan s2 = {(s1.t/dt+1)*dt, 
        x+dt/6*(k1v + 2*k2v + 2*k3v + k4v), 
        v+dt/6*(k1a + 2*k2a + 2*k3a + k4a),
        0};
    s2.e = m0*(phi(s2.x)+0.5*pow(s2.v,2));
    return s2;
}

void zapisz(Stan s, std::ofstream * out){
    *out << s.t << " " << s.x << " " << s.v << " " << s.e << std::endl;
}
