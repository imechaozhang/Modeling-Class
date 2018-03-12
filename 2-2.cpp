#include <iostream>
#include <fstream>
#include <cmath>
#include "random_mars.cpp"

using namespace std;

const double dt = 0.004;
const int Nt = 5001; //loop
const int Ne = 100;
const int SI = 50; //snapshot interval
const int M = 200;

const int Na = 125; //Number of atoms
const double rho = 0.2; //density
const double T = 0.5;
const double rc2 = 9;

double V = Na/rho; //volume
double a = pow(V,1.0/3);

class particle
{
public:
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double fx;
    double fy;
    double fz;
    double xr;
    double yr;
    double zr;
};

particle par[Na];

void init();
void loop();
void force();
double potential();
double potential(int);
double kinetic();
double pressure();
void test();

int main() {
    init();
    test();
    loop();
    return 0;
}

void init(){
    int i;
    int naa = ceil(pow(Na,1.0/3)); //number of atoms per dimension
    double v1 = 0.0, v2 = 0.0, v3 = 0.0;
    ofstream initp("initp.dat", ios::out);
    for (i=0;i<Na;i++){
        par[i].x = i/naa/naa * a/naa;
        par[i].y = (i/naa)%naa * a/naa;
        par[i].z = i%naa * a/naa;
        initp << par[i].x << " "<< par[i].y << " "<< par[i].z << endl;
    }
    initp.close();
    ofstream initv("initv.dat", ios::out);
    for (i=0;i<Na;i++){
        par[i].vx = RanMars(Na+i).gaussian()*pow(3*T, 0.5);
        v1 += par[i].vx;
        par[i].vy = RanMars(2*Na+i).gaussian()*pow(3*T, 0.5);
        v2 += par[i].vy;
        par[i].vz = RanMars(3*Na+i).gaussian()*pow(3*T, 0.5);
        v3 += par[i].vz;
    }
    for (i=0;i<Na;i++){
        par[i].vx -= v1/Na;
        par[i].vy -= v2/Na;
        par[i].vz -= v3/Na;
        initv << par[i].vx << " "<< par[i].vy << " "<< par[i].vz << endl;
    }
    initv.close();

}

void loop(){
    int i, t;
    double u, pr, ke, tot;
    int j1, j2, m, count=0;
    double x12=0.0, y12=0.0, z12 = 0.0, dr=0.0, rdf[M];
    int j;
    double msd=0, vac = 0;

    for(m=0;m<M; m++){
        rdf[m] = 0.0;
    }

    force();
    ofstream energy("energy.dat", ios::out);
    ofstream P("pressure.dat", ios::out);
    ofstream RDF("RDF.dat", ios::out);
    ofstream MSD("MSD.dat", ios::out);
    for(t=0;t<Nt;t++){
        for (i = 0; i < Na; i++) {
            particle temp = par[i];
            par[i].xr += par[i].vx*dt + 0.5*par[i].fx*dt*dt;
            par[i].yr += par[i].vy*dt + 0.5*par[i].fy*dt*dt;
            par[i].zr += par[i].vz*dt + 0.5*par[i].fz*dt*dt;

            par[i].x += par[i].vx*dt + 0.5*par[i].fx*dt*dt;
            par[i].y += par[i].vy*dt + 0.5*par[i].fy*dt*dt;
            par[i].z += par[i].vz*dt + 0.5*par[i].fz*dt*dt;
            par[i].x -= floor(par[i].x/a)*a;
            par[i].y -= floor(par[i].y/a)*a;
            par[i].z -= floor(par[i].z/a)*a;
            force();
            par[i].vx += 0.5*(par[i].fx+temp.fx)*dt;
            par[i].vy += 0.5*(par[i].fy+temp.fy)*dt;
            par[i].vz += 0.5*(par[i].fz+temp.fz)*dt;
        }

        if(t%SI == 0){
            cout << 100*t/Nt << "%" << endl;
            count ++;
            u = potential();
            ke  = kinetic();
            tot = u+ke;
            energy << dt*t << " " << u << " " << ke << " " << tot << endl;
            pr = pressure();
            P << dt*t/SI << " " << pr << endl;
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
            if (t>Ne){
                msd = 0;
                for (j=0;j<Na;j++){
                    msd += pow(par[j].xr-par[Ne].xr,2)+pow(par[j].yr-par[Ne].yr,2)+pow(par[j].zr-par[Ne].zr,2);
                    vac += par[j].vx*par[Ne].vx+par[j].vy*par[Ne].vy+par[j].vz*par[Ne].vz;
                }
                msd /= Na;

                MSD<< dt*(t-Ne) << " "<<msd<<endl;
            }
        }
    }
    vac /= 3.0*Na/dt;
    cout << vac << endl;
    for(m=0;m<M; m++){
        rdf[m] /= count*Na*3.14/2*pow(m,2) *pow(a/M,3)*rho;
        RDF << a/2.0*m/M << " " << rdf[m] << endl;
    }
    energy.close();
    P.close();
    RDF.close();
    MSD.close();

    ofstream final("final.dat", ios::out);
    for (i=0;i<Na;i++){
        final << par[i].x << " "<< par[i].y << " "<< par[i].z << endl;
    }
    final.close();
}

void force(){
    int i, j;
    double dx=0, dy=0, dz=0, dr2 = 0, dr=0, f=0;

    for (i =0;i<Na;i++){
        par[i].fx = 0.0;
        par[i].fy = 0.0;
        par[i].fz = 0.0;
        for (j = 0;j<Na;j++){
            if(i==j) continue;
            dx = par[j].x - par[i].x;
            dy = par[j].y - par[i].y;
            dz = par[j].z - par[i].z;
            dx -= a*floor(dx/a + 0.5);
            dy -= a*floor(dy/a + 0.5);
            dz -= a*floor(dz/a + 0.5);
            dr2 = dx*dx+dy*dy+dz*dz;
            if (dr2<rc2) {
                dr = sqrt(dr2);
                f = 24.0 * (pow(dr, -7.0) - 2.0 * pow(dr, -13.0));
                par[i].fx += f * dx / dr;
                par[i].fy += f * dy / dr;
                par[i].fz += f * dz / dr;
            }

        }
    }
}

double potential(){
    int k;
    double u = 0;
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
/*
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
*/
double pressure(){
    int i,j;
    double dx, dy, dz, dr=0, f, pr = 0;

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

double kinetic(){
    int i;
    double ke=0;

    for (i = 0; i< Na; i++){
        ke += pow(par[i].vx,2)+pow(par[i].vy,2)+pow(par[i].vz,2);
    }
    ke /= 2.0;
    return ke;
}

void test(){
    int i, j;
    double k1=1, k2;
    particle temp = par[1];
    for (i =0;i<Na;i++){
        for (j = i+1;j<Na;j++){
            //cout << i << " " << j << endl;
        }
    }
}
