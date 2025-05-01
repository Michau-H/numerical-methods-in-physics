#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

const double l1=1, l2=1/sqrt(8);
const double m0 = 1.;

struct Stan {
    long double t;
    double x,v,e;
};


double phi(double);

double a(double);

Stan war_pocz();

Stan Euler(Stan, double);

Stan Verlet(Stan, double);

Stan RK4(Stan, double);

void zapisz(Stan, std::ofstream *);


int main() {
    
    double tol = 1e-8;
    std::string name[3] = {"Euler", "Verlet", "RK4"};
    int d[3]={1,2,4};

    for(int k=0; k<3; k++){
        std::string n = name[k];
        n += ("_" + std::to_string(0) + ".txt");
        // sciezka folderu do zapisania
        std::string folder = "/Users/michau/Documents/MOFIT1/results_p01/results_2/";
        std::string path = folder + n; 
        std::ofstream outW(path);
        double dt = 0.001;
        int rozm = int(100./dt)+1;

        std::vector < Stan > test;

        test.push_back(war_pocz());
        zapisz(test[0], &outW);
            
        if(k==0){  
            while(test.back().t <100){
                Stan u2, u12;
                u2 = Euler(test.back(), 2*dt);
                u12 = Euler(Euler(test.back(), dt), dt);
                double eps_x = abs((u2.x-u12.x)/(pow(2,d[k])-1));
                double eps_v = abs((u2.v-u12.v)/(pow(2,d[k])-1));
                double eps = (eps_x>eps_v ? eps_x : eps_v);
                if(eps <= tol){
                    test.push_back(u12);
                    zapisz(test.back(), &outW);
                }
                dt = 0.9*dt*pow((tol/eps), 1./(d[k]+1));
            }
        }
        else if(k==1){
            while(test.back().t <100){
                Stan u2, u12;
                u2 = Verlet(test.back(), 2*dt);
                u12 = Verlet(Verlet(test.back(), dt), dt);
                double eps_x = abs((u2.x-u12.x)/(pow(2,d[k])-1));
                double eps_v = abs((u2.v-u12.v)/(pow(2,d[k])-1));
                double eps = (eps_x>eps_v ? eps_x : eps_v);
                if(eps <= tol){
                    test.push_back(u12);
                    zapisz(test.back(), &outW);
                }
                dt = 0.9*dt*pow((tol/eps), 1./(d[k]+1));
            }
        }
        else if(k==2){
            while(test.back().t <100){
                Stan u2, u12;
                u2 = RK4(test.back(), 2*dt);
                u12 = RK4(RK4(test.back(), dt), dt);
                double eps_x = abs((u2.x-u12.x)/(pow(2,d[k])-1));
                double eps_v = abs((u2.v-u12.v)/(pow(2,d[k])-1));
                double eps = (eps_x>eps_v ? eps_x : eps_v);
                if(eps <= tol){
                    test.push_back(u12);
                    zapisz(test.back(), &outW);
                }
                dt = 0.9*dt*pow((tol/eps), 1./(d[k]+1));
            }
        }
        test.clear();
        test.shrink_to_fit();

        outW.close();
        std::cout << "zapisane " << n << std::endl;
    }


    return 0;
}


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
