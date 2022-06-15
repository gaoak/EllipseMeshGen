//generating a part rectangular part curve near field mesh and rectangle wake mesh
//2016-12-04
//Gao ak
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "utility.h"
//#define ELLIPSE
#define CYLINDER
using namespace std;

int main(int argc, char* argv[])
{
    ofstream outxml, outtec;
    outxml.open("ellipse.xml");
    outtec.open("ellipse.dat");
	//parameters//////////////////////////////////
//	//part 1
//	int Np0 = 70;
//    double firstH = 0.01;
//    double R0 = 5.;
//    double deltaTheta;
//    int curvedpoints = 12;
//    int layer0, ncedge;
//    double **x0, **y0, ***cc;
//    int *cedge;
//    //part 2
//    double R1 = 50.;
//    double downsx = 12.;//may be changed in runtime
//    int layer1, iupr, idor;
//    double **x1, **y1;
//    //part 3
//    double ysminus = -50., ysplus = 50., yeminus = -50., yeplus = 50., wakeslen = 80.;
//    double hflows = 1., hflowe = 2.;
    //part 1
	int Np0 = 40;
    double firstH = 0.03;
    double R0 = 4.;
    double deltaTheta;
    int curvedpoints = 12;
    int layer0, ncedge;
    double **x0, **y0, ***cc;
    int *cedge;
    //part 2
    double R1 = 16.;
    double downsx = 10.;//may be changed in runtime
    int layer1, iupr, idor;
    double **x1, **y1;
    //part 3
    double ysminus = -R1, ysplus = R1, yeminus = -R1*1.0, yeplus = R1*1.0, wakeslen = 8.;
    double hflows = 1., hflowe = 0.6;
    int is3, ie3, layer2, Np2;
    double wakeplus = (yeplus-yeminus)*0.4;
    int layer2plus;
    ///////////////////////////////////////////////////////////////////////////////////////
//generating points
#ifdef ELLIPSE
    //Np0 = 90;//78 + 12
    //int Nsg = 58, Ns0 = 9;//Np0 total circumferential mesh; Nsg total number used for smooth mesh; Ns0*4 is the real smooth mesh number
    Np0 = 90+2+4+2;//54 + 12
    int Nsg = 78, Ns0 = 8;//Np0 total circumferential mesh; Nsg total number used for smooth mesh; Ns0*4 is the real smooth mesh number
    //double ella = 0.5, ellb = 0.0625;
    double ella = 0.5, ellb = 0.03125;
    firstH = 0.005;
    ellbl_pointsgen(ella, ellb, Np0, Nsg, Ns0, R0, firstH, curvedpoints, deltaTheta, layer0, x0, y0, ncedge, cedge, cc);
    double alpha = 0.;//140./180.*M_PI;
    ///test
    //alpha = M_PI/2.;
    ///
    coord_trans(Np0, layer0, alpha, x0, y0);
    cedge_trans(ncedge, curvedpoints, alpha, cc);
#else
    Np0 = 104;
    firstH = 2.2E-4;
    double rc = 0.5;
    cirbl_pointsgen(rc, Np0, Np0-10, 11, R0, firstH, curvedpoints, deltaTheta, layer0, x0, y0, ncedge, cedge, cc);
#endif
    circle_rect(Np0, layer0, R0, deltaTheta, x0, y0, R1, downsx, layer1, iupr, idor, x1, y1);
    //merge first two layers
    double **x, **y;
    int layers = layer0+layer1;
    x = new double*[layers];
    y = new double*[layers];
    for(int i=0; i<layer0; i++)
    {
        x[i] = x0[i];
        y[i] = y0[i];
    }
    for(int i=0; i<layer1; i++)
    {
        x[layer0+i] = x1[i];
        y[layer0+i] = y1[i];
    }
    delete[] x0;
    delete[] y0;
    delete[] x1;
    delete[] y1;
    outtecplot_cir(Np0, layers, x, y, outtec);
    outtec.close();

    double **x2, **y2;
    wake_rect(downsx, layers, idor, iupr, x, y, is3, ie3, Np2, wakeslen, ysplus, ysminus, yeplus, yeminus, layer2, hflows, hflowe, x2, y2, wakeplus, layer2plus);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //layer2 = 1; layer2plus = 1;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double **point;
    int **cell, **edge;
    int npoint, ncell, nedge;
    npoint = layers*Np0 + (layer2+layer2plus)*Np2;
    ncell = (layers-1)*Np0 + (layer2+layer2plus)*(ie3-is3);
    nedge = (2*layers-1)*Np0 + (2*Np2-1)*(layer2+layer2plus);
    point = new double*[npoint];
    edge = new int*[nedge];
    cell = new int*[ncell];
    for(int i=0; i<npoint; i++) point[i] = new double[2];
    for(int i=0; i<nedge; i++) edge[i] = new int[2];
    for(int i=0; i<ncell; i++) cell[i] = new int[4];
    FEgen_cir(Np0, layers, x, y, point, 0, edge, 0, cell);
    FEgen_rec(Np0, layers, is3, ie3, Np2, (layer2+layer2plus), x2, y2, point+layers*Np0, edge+(2*layers-1)*Np0, cell+(layers-1)*Np0);
    ///////////////////////output session file////////////////////
    outFE(npoint, point, nedge, edge, curvedpoints, ncedge, cedge, cc, ncell, cell, outxml);
    outCOMPO(Np0, layers, is3, ie3, Np2, layer2, layer2plus, cell, outxml);
    outxml.close();
    cout << "number of cells " << ncell << endl;
    cout << Np0 << ", " << layer0 << endl;
    cout << Np0 << ", " << layer1 << endl;
    cout << Np2 << ", " << (layer2+layer2plus) << endl;
    return 0;
}
