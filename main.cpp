#include <iostream>
#include <string>
#include "ctime"
#include "cstdlib"
#include <fstream>
#include <cmath>
#include "random_mars.h"
#include "random_mars.cpp"

using namespace std;

//Constant variables
const int nv = 2;
const int N = 100000; //loop times
const double d = 1.0 / 14;
double d0 = 0;  //d0 and alpha calcuted later
double alpha = 0;

//For RDFs
double K = 0; //Maximum distance
const int M = 6400;  //Number of areas divided

int potential(int, int, double[16][14], double[16][14]);
int main() {
    double a1, a2; // random variables
    int move, i, j, k, l, m; // for loop variables
    double x[16][14], y[16][14], xt, yt, xp, yp, r, xij, yij; //coordinates
    d0 = d*(1 - pow(2, nv - 8));
    K = 0.5/d0;
    alpha = d - d0;
    //RDF variables
    double rdf[M], r2[M];
    int count = 0;

    //RDF initialization
    for (i = 0; i < M; i++) {
        r2[i] = pow(d0, 2)*(1 + i *(pow(K, 2) - 1) / M);
        rdf[i] = 0.0;
    }

    //coordinate initialization, using rounding number to get the crisscrossing structure, and save to init.dat for ploting
    ofstream myfile("init.dat", ios::out);
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 14; j++) {
            x[i][j] = 1.0 / 14 * ((i % 2) / 2.0 + j);
            y[i][j] = 1.0 / 16 * i;
            myfile << x[i][j] << ' ' << y[i][j] << endl;
        }
    }
    myfile.close();

    //Main loop
    for (move = 0; move<N; move++) {
        //Two set of random numbers, one for moves, and one for picking one particle for movement
        a1 = RanMars(2 * move + 1).uniform();
        a2 = RanMars(2 * move + 2).uniform();
        k = rand() % 16;
        l = rand() % 14;

        //in case of rejecting the move, using temperary variables
        xt = x[k][l];
        yt = y[k][l];

        //move it
        xt += 2*alpha*(a1 - 0.5);
        yt += 2*alpha*(a2 - 0.5);

        //rounding, to make the particle in box
        xt -= floor(xt);
        yt -= floor(yt);

        //calculating the potential energy, accept if zero and reject if infinity
        for (i = 0; i < 16; i++) {
            for (j = 0; j < 14; j++) {
                xp = x[i][j];
                yp = y[i][j];

                if (i == k && j == l) {
                    continue;
                } //ignore itself

                if (xp - xt>0.5) {
                    xp--;
                }
                if (xp - xt<-0.5) {
                    xp++;
                }
                if (yp - yt>0.5) {
                    yp--;
                }
                if (yp - yt<-0.5) {
                    yp++;
                }

                r = sqrt(pow(xp - xt, 2) + pow(yp - yt, 2));

                //reject the move if infinite potential energy
                if (r<d0) {
                    xt = x[k][l];
                    yt = y[k][l];
                }
            }
        }
        x[k][l] = xt;
        y[k][l] = yt;

        //RDF loop
        //production after equilibrium, and collect data for each 16 cycles, as in the paper
        if ( move > 10000 && move % 2240 == 0) {
            count++; //counting the number of average times
                for (i = 0; i < 16; i++) {
                    for (j = 0; j < 14; j++) {
                        for (k = 0; k < 16; k++) {
                            for (l = 0; l < 14; l++) {
                                xp = x[i][j];
                                yp = y[i][j];
                                xt = x[k][l];
                                yt = y[k][l];
                                if (i == k && j == l) {
                                    continue;
                                }
                                if (xp - xt>0.5) {
                                    xp--;
                                }
                                if (xp - xt<-0.5) {
                                    xp++;
                                }
                                if (yp - yt>0.5) {
                                    yp--;
                                }
                                if (yp - yt<-0.5) {
                                    yp++;
                                }
                                m = (int)((pow(xp-xt, 2) + pow(yp-yt, 2)) /( pow(d0, 2)*(pow(K, 2) - 1) / M));
                                cout << m << endl;
                                if (m>=0 && m<M) {
                                    rdf[m] = rdf[m] + 1; //counting the particles in the bins
                                }
                            }
                        }
                    }
                }
        }
    }

    //Normalize the counting and save it to file
    ofstream myrdf("rdf.dat", ios::out);
    for (m = 0; m < M; m++) {
        rdf[m] /= count * 224 * 224 * 3.14 * pow(d0, 2) * (pow(K, 2) - 1) / M;
        myrdf << sqrt((1 + m*(pow(K, 2) - 1) / M))<<" " << rdf[m] << endl;
    }
    myrdf.close();

    //save the final configuration to the file
    ofstream myfinalfile("final.dat", ios::out);
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 14; j++) {

            myfinalfile << x[i][j] << " " << y[i][j] << endl;
        }
    }
    myfinalfile.close();
    return 0;
}
