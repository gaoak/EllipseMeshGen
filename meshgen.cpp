//generating a part rectangular part curve near field mesh and rectangle wake mesh
//2016-12-04
//Gao ak
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "utility.h"
#define ELLIPSE
//#define CYLINDER
using namespace std;

int main(int argc, char* argv[])
{
    ofstream outxml, outtec;
    string fname;
#ifdef ELLIPSE
    fname = "ellipse.xml";
#else
    fname = "cylinder.xml";
#endif
    outxml.open(fname);
    outtec.open("ellipse.dat");
	//parameters//////////////////////////////////
    //part 1
	int Np0 = 40;
    double firstH = 0.03;
    double R0 = 3.;
    double deltaTheta;
    int curvedpoints = 12;
    int layer0, ncedge;
    double **x0, **y0, ***cc;
    int *cedge;
    ///////////////////////////////////////////////////////////////////////////////////////
//generating points
#ifdef ELLIPSE
    int Nsg = 78, Ns0 = 8;//Np0 total circumferential mesh; Nsg total number used for smooth mesh; Ns0*4 is the real smooth mesh number
    Np0 = Nsg+30;//54 + 12
    double ella = 0.5, ellb = 0.0625;
    firstH = 0.002;
    ellbl_pointsgen(ella, ellb, Np0, Nsg, Ns0, R0, firstH, curvedpoints, deltaTheta, layer0, x0, y0, ncedge, cedge, cc);
    double alpha = 0.;//140./180.*M_PI;
    ///test
    alpha = 170./180.*M_PI;
    ///
    coord_trans(Np0, layer0, alpha, x0, y0);
    cedge_trans(ncedge, curvedpoints, alpha, cc);
#else
    Np0 = 50;
    firstH = 0.015;
    double rc = 0.5;
    cirbl_pointsgen(rc, Np0, Np0-12, 11, R0, firstH, curvedpoints, deltaTheta, layer0, x0, y0, ncedge, cedge, cc);
#endif

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double **point;
    int **cell, **edge;
    int npoint, ncell, nedge;
    int layers = layer0;
    npoint = layers*Np0 ;
    ncell = (layers-1)*Np0;
    nedge = (2*layers-1)*Np0;
    point = new double*[npoint];
    edge = new int*[nedge];
    cell = new int*[ncell];
    for(int i=0; i<npoint; i++) point[i] = new double[2];
    for(int i=0; i<nedge; i++) edge[i] = new int[2];
    for(int i=0; i<ncell; i++) cell[i] = new int[4];
    FEgen_cir(Np0, layers, x0, y0, point, 0, edge, 0, cell);
    ///////////////////////output session file////////////////////
    outFE(npoint, point, nedge, edge, curvedpoints, ncedge, cedge, cc, ncell, cell, outxml);
    outCOMPO(Np0, layers, cell, outxml);
    outxml.close();
    cout << "number of cells " << ncell << endl;
    cout << Np0 << ", " << layer0 << endl;
    return 0;
}
