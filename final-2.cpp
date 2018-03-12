#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random_mars.cpp"

//Path Integral Metadynamics

using namespace std;

const int P = 16; //Number of beads
const double T = 0.1; //Temperature
const int Nt = 1000000; //Time Integration
const int SI = 500;//snapshot interval
const double dt = 0.02; //Time Intevals
const double gam = 5; //gamma
const double w0 = 2.0*T; //initial gaussian height
const double b = 1.0/64; //double well parameter

double sigma2 = 4*T*T; //square of gaussian width
double beta = 1.0/T; //inverse temperature
double w = P/beta/beta; //spring stiffness
double Tw = (gam - 1)*T; //a temperature for gaussian bias potentials

//particle class, position, velocity and force.
class particle {
public:
    double x;
    double v;
    double f;
};

particle par[P]; //P-bead ring polymer
vector <double> st; //Gaussian bias potential - collective variable
vector <double> wt;//Gaussian bias potential - height

void init(); //Initialization
void force(double); //calculating the forces, related to s
double kinetic(); //calculating the kinetic energy for renormaliztion
void renorm(double); //renormalization
double cv(); //calculating collective variable s
double heightg(double); //calculating gaussian height

int main() {
    double ke, s = 0, hg, u;
    double hist[100]; //position probabilty histogram
    double vb[100], si; //biased potential
    int count;

    for (auto &i : hist){
        i = 0;
    }

    init();//Initialization
    ofstream energy("energy.dat", ios::out);
    //Time Integration
    for (int t = 0; t < Nt; t++){
        //velocity verlet
        for (auto &i : par) {
            particle temp = i;
            i.x += i.v*dt + 0.5* i.f*dt*dt;
            s = cv();
            force(s);
            i.v += 0.5*(i.f+temp.f)*dt;
        }

        if (0 == t*100 % Nt) cout << 100 * t / Nt + 1 << "%" << endl; //monitoring the process

        //collection
        if (0 == t % SI){
            hg = heightg(s);
            st.push_back(s);
            wt.push_back(hg);
            ke  = kinetic();
            renorm(ke);
            energy << dt*t << " " << ke << ' ' << s/P/T << endl;
            /*if (t >= Nt/2){
                for (auto &i : par){
                    count = static_cast<int>((i.x + 8) / 16 * 100);
                    hist[count] += 1.0/P/Nt*SI*2*100/16;
                }
            }*/
        }
    }
    energy.close();
    ofstream position("position.dat", ios::out);
    ofstream FreeEnergy("FreeEnergy.dat", ios::out);

    for (int i = 0; i< 100; i++){
        position << i/100.0*16-8 << ' ' << hist[i] << endl;
    }
    for (int j = 0; j < 100; j++){
        si = j*2.0/100*P*T;
        vb[j] = 0;
        for (int i = 0; i < wt.size(); i++){
            vb[j] += wt[i]*exp(-1*pow(si-st[i],2)/2/sigma2);
        }
        FreeEnergy << si/P/T << ' '<< (1.0/gam-1)*vb[j] << endl;
    }
    position.close();
    FreeEnergy.close();
    return 0;
}

//Initialization, random velocity according to T. original positions.
void init(){
    double va = 0;
    for (int i = 0; i < P; i++){
        par[i].x = 0.0;
        par[i].f = 0.0;
        par[i].v = RanMars(1+i).gaussian()*sqrt(T);
        va += par[i].v;
    }
    for (auto &i : par) {
        i.v -= va/P;
    }
}

//Calculating the forces
void force(double s){
    double vs = 0;

    //an extra force coming from the biased potential, vs = dv/ds
    for (int i = 0; i < wt.size(); i++){
        vs += (s - st[i])/sigma2*wt[i]*exp(-1*pow(s-st[i],2)/2/sigma2);
    }
    for (int i = 0; i< P; i++){
        if (i == P-1){
            par[i].f = (1-vs)*w*(par[0].x - 2*par[i].x + par[i-1].x) - par[i].x/P; //harmonic
            //par[i].f = (1-vs)*w*(par[0].x - 2*par[i].x + par[i-1].x) + (par[i].x-4*b*pow(par[i].x,3))/P; //double well
            continue;
        }
        par[i].f = (1-vs)*w*(par[i+1].x - 2*par[i].x + par[i-1].x) - par[i].x/P; //harmonic
        //par[i].f = (1-vs)*w*(par[i+1].x - 2*par[i].x + par[i-1].x) + (par[i].x-4*b*pow(par[i].x,3))/P; //double well
    }
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

//calculating the collective variable s
double cv(){
    double s = 0, d;
    for(int i = 0;i < P; i++){
        if (i == P-1) {
            d = par[i].x - par[0].x;
        }
        else {
            d = par[i].x - par[i + 1].x;
        }
        s += 0.5*w*d*d;
    }
    return s;
}

//the time-dependent gaussian height
double heightg(double s){
    double wg, vg = 0;

    for (int i = 0; i < wt.size(); i++){
        vg += wt[i]*exp(-1*pow(s-st[i],2)/2/sigma2);
    }
    wg = w0*exp(-1*vg/Tw);
    return wg;
}
