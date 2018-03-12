#include <iostream>
#include <fstream>
#include <cmath>
#include "random_mars.cpp"

//Path Integral Molecular Dynamics

using namespace std;

const int P = 32; //Number of beads
const double T = 0.1; //Temperature
const int Nt = 500000; //Time Integration
const int SI = 300;//snapshot interval
const double dt = 0.01; //Time Intevals
const double b = 1.0/64; //parameter for double well

double beta = 1.0/T; //inverse temperature, kB = 1
double w = P/beta/beta; //stiffness of the spring

//particle class, position, velocity and force.
class particle {
public:
    double x;
    double v;
    double f;
};

particle par[P]; //P-bead ring polymer

void init(); //Initialization
void force(); //calculating the forces
double potential(); //calculating the potential energies
double kinetic(); //calculating the kinetic energies
void renorm(double); //renormalization velocities

int main() {
    double u, ke, tot; //potential, kinetic, and total energy
    double hist[100]; //position probability histogram
    int count;

    for (auto &i : hist){
        i = 0;
    }

    init();
    ofstream energy("energy.dat", ios::out);
    //Time Integration
    for (int t = 0; t < Nt; t++){
        //velocity verlet algorithm
        for (auto &i : par) {
            particle temp = i;
            i.x += i.v*dt + 0.5* i.f*dt*dt;
            force();
            i.v += 0.5*(i.f+temp.f)*dt;
        }
        if (0 == t*100 % Nt) cout << 100*t/Nt + 1 << "%" << endl; //monitoring the process
        //Information collection
        if (0 == t % SI){
            u = potential();
            ke  = kinetic();
            renorm(ke);
            tot = u+ke;
            energy << dt*t << " " << u << " " << ke << " " << tot << endl;
            if (t >= Nt/2){ //histogram
                for (auto &i : par){
                    count = static_cast<int>(floor((i.x + 8) / 16 * 100));
                    hist[count] += 1.0/P/Nt*SI*2*100/16;
                }
            }
        }
    }

    ofstream position("position.dat", ios::out);
    for (int i = 0; i< 100; i++){
        position << i/100.0*16-8 << ' ' << hist[i] << endl;
    }
    energy.close();
    position.close();
    return 0;
}

//Initialization, random velocity according to T. original positions.
void init(){
    double va = 0;
    for (int i = 0; i < P; i++){
        par[i].x = -4.0;
        par[i].f = 0.0;
        par[i].v = RanMars(P+i).gaussian()*sqrt(T);
        va += par[i].v;
    }
    for (auto &i : par) {
        i.v -= va/P;
    }
}

//Calculating the forces
void force(){
    for (int i = 0; i< P; i++){
        if (i == P-1){
            par[i].f = w*w*(par[0].x - 2*par[i].x + par[i-1].x) - par[i].x/P; //Harmonic
            //par[i].f = w*w*(par[0].x - 2*par[i].x + par[i-1].x) + (par[i].x-4*b*pow(par[i].x,3))/P; //Double well
            continue;
        }
        par[i].f = w*w*(par[i+1].x - 2*par[i].x + par[i-1].x) - par[i].x/P; //Harmonic
        //par[i].f = w*w*(par[i+1].x - 2*par[i].x + par[i-1].x) + (par[i].x-4*b*pow(par[i].x,3))/P; //double well
    }
}


double potential(){
    double d,pot = 0;

    //Internal spring potential
    for(int i = 0;i < P; i++){
        if (i == P-1) {
            d = par[i].x - par[0].x;
        }
        else {
            d = par[i].x - par[i + 1].x;
        }
        pot += 0.5*w*w*d*d;
    }
    //external potential
    for (auto &i : par){
        //pot += 0.5*i.x*i.x/P;
        pot += (-0.5*i.x*i.x + b*pow(i.x,4))/P;
    }
    return pot;
}

double kinetic(){
    double ke = 0;

    for (auto &i : par) {
        ke += 0.5*i.v*i.v;
    }
    return ke;
}

void renorm(double ke){
    for (auto &i : par){
        i.v *= sqrt(T*P/ke/2);
    }
}
