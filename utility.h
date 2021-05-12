#include <cmath>
#include <fstream>

using namespace std;
int outtecplot_cir(int Np0, int layers, double **x, double **y, ofstream &outtec);
int outtecplot_rec(int Np0, int layers, double **x, double **y, ofstream &outtec);
int ellbl_pointsgen(double ella, double ellb, int Np0, int Nsg, int Ns0, double R0, double firstH, int curvedpoints,
                    double &deltaTheta, int &layers, double** &x, double** &y, int &ncedges, int* &cedges, double *** &cc);
int cirbl_pointsgen(double rc, int Np0, int Nsg, int Ns0, double R0, double firstH, int curvedpoints,
                    double &deltaTheta, int &layers, double** &x, double** &y, int &ncedges, int* &cedges, double *** &cc);
int circle_rect(int Np0, int layer0, double R0, double deltaTheta, double **x0, double **y0, double R1,
                double &downslen, int &layer1, int &iupr, int &idor, double **&x1, double **&y1);
int wake_rect(double downslen, int layers, int idor, int iupr, double **x, double **y,
                int &is3, int &ie3, int &Np2, double &wakeslen, double &ysplus, double &ysminus, double yeplus, double yeminus,
                int &layer2, double &hflows, double &hflowe, double **&x2, double **&y2, double &wakeplus, int &layer2plus);
int FEgen_rec(int Np0, int layers, int is3, int ie3, int Np2, int layer2, double **x2, double **y2, double **point, int **edge, int **cell);
int FEgen_cir(int Np0, int layers, double **x, double **y, double **point, int nspoint, int **edge, int nsedge, int **cell);
int outFE(int npoint, double **point, int nedge, int **edge, int curvedpoints, int ncedge, int *cedge, double ***cc, int ncell, int **cell, ofstream &outxml);
int outCOMPO(int Np0, int layers, int **cell, ofstream &outxml);
int coord_trans(int Np0, int layers, double alpha, double **x, double **y);
int cedge_trans(int ncedge, int curvedpoints, double alpha, double ***cc);
inline double fStretch(double mu, double mu0)
{
    double lambda = 2.;
    double d = exp(mu) - exp(-mu) - (exp(mu0) - exp(-mu0));
    return 1. - exp(-lambda*d);
}
inline double fs(double x)
{
    return pow((x/0.3),2);
}
inline double meshh(double Dr, double h0, double hi)
{
	return (h0+hi*fs(Dr))/(1.+fs(Dr));
}
