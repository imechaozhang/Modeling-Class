#include <iostream>
#include <string>
#include "ctime"
#include "cstdlib"
#include <fstream>
#include <cmath>
#include "random_mars.cpp"

using namespace std;

const int Np = 1000001; //loop
const int SI = 500; //snapshot interval
const int M  = 200;

const int Na = 125; //Number of atoms
const double rho = 0.2; //density
const double T = 1.0; //temperature
const double rc2 = 9.0;

double V = Na/rho; //volume
double a = pow(V,1.0/3);

class particle
{
    public:
        double x;
        double y;
        double z;
};

particle par[Na];

void init();
//void equilibrium();
void loop();
double potential();
double potential(int);
double potential(int, particle);
double pressure();
void test();

int main() {
    init();
    test();
    //equilibrium();
    loop();
    return 0;
}

void init(){
    int i;
    int naa = ceil(pow(Na,1.0/3)); //number of atoms per dimension
    ofstream myfile("init.dat", ios::out);
    for (i=0;i<Na;i++){
        par[i].x = i/naa/naa * a/naa;
        par[i].y = (i/naa)%naa * a/naa;
        par[i].z = i%naa * a/naa;
        myfile << par[i].x << " "<< par[i].y << " "<< par[i].z << endl;
    }
    myfile.close();
}

void loop(){
    int i, k, j1, j2, m, count=0;
    double p, dx, dy, dz, alpha = pow(rho,-1.0/3);
    double x12=0.0, y12=0.0, z12 = 0.0, dr=0.0;
    double u0, u1, u, pr;
    double rdf[M];
    particle temp;

    for(m=0;m<M; m++){
        rdf[m] = 0.0;
    }

    ofstream P("pressure.dat", ios::out);
    ofstream U("potential.dat", ios::out);
    ofstream RDF("RDF.dat", ios::out);
    for (i = 0; i < Np; i++){
        k = floor(RanMars(6*Np+i).uniform() * Na);
        dx = (RanMars(7*Np+i).uniform()-0.5)*alpha;
        dy = (RanMars(8*Np+i).uniform()-0.5)*alpha;
        dz = (RanMars(9*Np+i).uniform()-0.5)*alpha;

        temp = par[k];
        temp.x += dx;
        temp.x -= floor(temp.x/a)*a;
        temp.y += dy;
        temp.y -= floor(temp.y/a)*a;
        temp.z += dz;
        temp.z -= floor(temp.z/a)*a;

        p = RanMars(10*Np+i).uniform();
        u0 = potential(k);
        u1 = potential(k,temp);
        if(p<exp((u0-u1)/T)){
            par[k] = temp;
        }

        if(i%SI == 0){
            count ++;
            u = potential();
            U << i/SI << " " << u << endl;
            pr = pressure();
            P << i/SI << " " << pr << endl;
            for (j1=0;j1<Na;j1++){
                for (j2=0;j2<Na;j2++){
                    if (j1==j2) continue;
                    x12 = par[j1].x - par[j2].x;
                    y12 = par[j1].y - par[j2].y;
                    z12 = par[j1].z - par[j2].z;
                    x12 -= a*floor(x12/a + 0.5);
                    y12 -= a*floor(y12/a + 0.5);
                    z12 -= a*floor(z12/a + 0.5);
                    dr = sqrt(x12*x12+y12*y12+z12*z12);
                    m = (int)(2*M*dr/a);
                    if(m<M)  {
                        rdf[m]++;
                    }
                }
            }
        }
    }
    for(m=0;m<M; m++){
        rdf[m] /= count*Na*3.14/2*pow(m,2) *pow(a/M,3)*rho;
        RDF << a/2.0*m/M << " " << rdf[m] << endl;
    }

    P.close();
    U.close();
    RDF.close();

    ofstream final("final.dat", ios::out);
    for (i=0;i<Na;i++){
        final << par[i].x << " "<< par[i].y << " "<< par[i].z << endl;
    }
    final.close();
}

double potential(){
    int k;
    double u = 0.0;
    for (k = 0; k < Na; k++) {
        u += potential(k);
    }
    u = u/2 + 8.0/3*3.14*rho*Na*(1.0/3*pow(rc2,-9.0/2)-pow(rc2,-3.0/2));
    return u;
}
double potential(int k){
    int i;
    double dx, dy, dz, dr2, u = 0;
    for (i = 0; i < Na; i++){
        if (i == k) continue;
        dx = par[i].x - par[k].x;
        dy = par[i].y - par[k].y;
        dz = par[i].z - par[k].z;
        dx -= floor(dx/a + 0.5)*a;
        dy -= floor(dy/a + 0.5)*a;
        dz -= floor(dz/a + 0.5)*a;

        dr2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
        if (dr2<rc2) {
            u += 4 * (pow(dr2, -6.0) - pow(dr2, -3.0));
        }
    }
    return u;
}
double potential(int k, particle temp){
    int i;
    double dx, dy, dz, dr2, u = 0;
    for (i = 0; i < Na; i++){
        if (i == k) continue;
        dx = par[i].x - temp.x;
        dy = par[i].y - temp.y;
        dz = par[i].z - temp.z;
        dx -= floor(dx/a + 0.5)*a;
        dy -= floor(dy/a + 0.5)*a;
        dz -= floor(dz/a + 0.5)*a;

        dr2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
        if (dr2<rc2) {
            u += 4 * (pow(dr2, -6.0) - pow(dr2, -3.0));
        }
    }
    return u;
}

double pressure(){
    int i,j;
    double dx, dy, dz, dr, f, pr = 0;

    for (i =0;i<Na;i++){
        for (j = i+1;j<Na;j++){
            dx = par[i].x - par[j].x;
            dy = par[i].y - par[j].y;
            dz = par[i].z - par[j].z;
            dx -= floor(dx/a + 0.5)*a;
            dy -= floor(dy/a + 0.5)*a;
            dz -= floor(dz/a + 0.5)*a;

            dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
            f = 24.0*(pow(dr,-7.0)-2.0*pow(dr, -13.0));
            pr += f*dr;
        }
    }
    pr = pr/3.0/V + rho*T + 16.0/3*pow(rho,2)*(2.0/3*pow(dr,-9)-pow(dr,-3));

    return pr;
}

void test(){
    int i, j;
    double k;
    particle temp = par[1];
    for (i =0;i<Na;i++){
        for (j = i+1;j<Na;j++){
            //cout << i << " " << j << endl;
        }
    }
}
