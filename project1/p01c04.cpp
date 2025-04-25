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


// funkcje
double phi(double);

double a(double, double, double);

double ddPhi(double, double);

double F1(Stan, Stan, double);

double F2(Stan, Stan, double, double);

Stan war_pocz();

Stan Euler(Stan, double, double);

// Stan Verlet(Stan, double);

Stan RK4(Stan, double, double);

Stan Trapezy(Stan, double, double);

void zapisz(Stan, std::ofstream *);


// main
int main() {
    
    double tol = 1e-7;
    std::string name = {"Trapezy"};
    int d={2};

    double alpha_tab[2] = {0.5, 5.0};
    for(int i=0; i<2; i++){

        double alpha = alpha_tab[i];
        std::string n = name;
        n += ("_" + std::to_string(i) + ".txt");
        std::string folder = "/Users/michau/Documents/MOFIT1/results_p01/results_4/";
        std::string path = folder + n; 
        std::ofstream outW(path);
        double dt = 0.001;
        int rozm = int(100./dt)+1;
    
        std::vector < Stan > test;
        test.push_back(war_pocz());
        zapisz(test[0], &outW);
        
          
        while(test.back().t <100){
            Stan u2, u12;
            u2 = Trapezy(test.back(), 2*dt, alpha);
            u12 = Trapezy(Trapezy(test.back(), dt, alpha), dt, alpha);
            double eps_x = abs((u2.x-u12.x)/(pow(2,d)-1));
            double eps_v = abs((u2.v-u12.v)/(pow(2,d)-1));
            double eps = (eps_x>eps_v ? eps_x : eps_v);
            if(eps <= tol){
                u12.t = test.back().t+2*dt;
                test.push_back(u12);
                zapisz(test.back(), &outW);
            }
            dt = 0.9*dt*pow((tol/eps), 1./(d+1));
        }
            

            
        test.clear();
        test.shrink_to_fit();

        outW.close();
        std::cout << "zapisane " << n << std::endl;

    
    }


    return 0;
}

// definicje
double phi(double x){
    return -pow(M_E,(-x*x/l1*l1))-8*pow(M_E,(-pow((x-2),2)/(l2*l2)));
}

double a(double x, double v, double alpha){
    return -1/m0*(phi(x + 0.001) - phi(x - 0.001))/(2*0.001) - alpha*v;
}

double ddPhi(double x, double alpha){
    return (phi(x + 0.001) - 2*phi(x) + phi(x - 0.001))/(0.001*0.001);
}

double F1(Stan s0, Stan s1, double dt){
    return s1.x - s0.x - dt/2*(s1.v + s0.v);
}

double F2(Stan s0, Stan s1, double dt, double alpha){
    return s1.v - s0.v - dt/2*(a(s1.x, s1.v, alpha) + a(s0.x, s0.v, alpha));
}

Stan war_pocz(){
    Stan s = {0,2.8,0,0};
    s.e = m0*(phi(s.x)+0.5*pow(s.v,2));
    return s;
}

Stan Euler(Stan s1, double dt, double alpha){
    Stan s2 = {(s1.t/dt+1)*dt, s1.x+dt*s1.v, s1.v+dt*a(s1.x, s1.v, alpha), 0};
    s2.e = m0*(phi(s2.x)+0.5*pow(s2.v,2));
    return s2;
}

Stan Trapezy(Stan s1, double dt, double alpha){
    Stan s2 = s1;
    s2.t = (s1.t/dt+1)*dt;
    for(int i=0; i<5; i++){
        double f1 = F1(s1, s2, dt);
        double f2 = F2(s1, s2, dt, alpha);
        double ww = 1 + dt/2*(alpha + dt/(2*m0)*ddPhi(s2.x, alpha));
        double wx = -f1*(1+dt/2*alpha) - f2*dt/2*alpha;
        double wv = dt/(2*m0)*ddPhi(s2.x, alpha)*f1 - f2;

        s2.x += wx/ww;
        s2.v += wv/ww;
    }
    s2.e = m0*(phi(s2.x)+0.5*pow(s2.v,2));
    return s2;
}

// Stan Verlet(Stan s1, double dt){
//     Stan s2 = {(s1.t/dt+1)*dt, s1.x+dt*s1.v+0.5*a(s1.x)*dt*dt, 0, 0};
//     s2.v = s1.v+0.5*dt*(a(s1.x)+a(s2.x));
//     s2.e = m0*(phi(s2.x)+0.5*pow(s2.v,2));
//     return s2;
// }

Stan RK4(Stan s1, double dt, double alpha){
    double x = s1.x, v = s1.v;
    double k1v = v, k1a=a(x, v, alpha);
    double k2v = v + 0.5*dt*k1a, k2a = a(x + 0.5*dt*k1v, v + 0.5*dt*k1a, alpha);
    double k3v = v + 0.5*dt*k2a, k3a = a(x + 0.5*dt*k2v, v + 0.5*dt*k2a, alpha);
    double k4v = v + dt*k3a, k4a = a(x + dt*k3v, v + dt*k3a, alpha);
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
