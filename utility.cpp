#include <iostream>
#include <fstream>
#include <iomanip>
#include "utility.h"

using namespace std;
static int refinethetas(int Np0, int Nplus, double alpha, double deltaa, double *thetas)
{
    double *thetaplus;
    int is, ie, Npa = Np0 + Nplus;
    for(is=0; is<Np0; is++)
    {
        if(thetas[is]>=alpha-deltaa) break;
    }
    is--;
    for(ie=0; ie<Np0; ie++)
    {
        if(thetas[ie]>=alpha+deltaa) break;
    }
    int Ntp = ie - is - 1 + Nplus;
    thetaplus = new double[Ntp];
    double subdelta = (thetas[ie] - thetas[is])/(Ntp + 1.);
    for(int i=1; i<=Ntp; i++) thetaplus[i-1] = thetas[is] + i*subdelta;
    for(int i=Np0-1; i>=ie; i--) thetas[Npa-1 - (Np0-1-i)] = thetas[i];
    for(int i=1; i<=Ntp; i++) thetas[is+i] = thetaplus[i-1];
    return 0;
}
int cirbl_pointsgen(double rc, int Np0, int Nsg, int Ns0, double R0, double firstH, int curvedpoints,
                    double &deltaTheta, int &layers, double** &x, double** &y, int &ncedge, int* &cedges, double *** &cc)
{//Np0, Nsg is even number, but Nsg is not quadruple, Ns0 means subdivide mesh numbers, which is odd number
    //circumferential mesh
    deltaTheta = 2.*M_PI/Nsg;
    double deltaTheta2 = 0.6*deltaTheta;
    double *thetas = new double[Np0+1];
    /*
    thetas[0] = M_PI*0.5; thetas[Np0] = M_PI*2.5;
    thetas[Nsg/2] = M_PI*1.5;
    for(int j=1; j<=(Nsg-2)/4; j++)
    {
        double dt = deltaTheta;
        thetas[j] = thetas[j-1] + dt;
        thetas[Nsg/2 - j] = thetas[Nsg/2 - j + 1] - dt;
        thetas[Nsg/2 + j] = thetas[Nsg/2 + j -1] + dt;
        thetas[Np0 - j] = thetas[Np0 - j + 1] - dt;
    }
    double subdeltaTheta = Ns0*deltaTheta/double(Np0-Nsg+Ns0);
    for(int i=1; i<=(Np0-Nsg + Ns0)/2; i++)
    {
        thetas[ (Nsg/2 - Ns0)/2 + Nsg/2 + i] = thetas[(Nsg/2 - Ns0)/2 + Nsg/2] + subdeltaTheta*i;
        thetas[-(Nsg/2 - Ns0)/2 + Np0 - i] = thetas[-(Nsg/2 - Ns0)/2 + Np0] - subdeltaTheta*i;
    }
    */
    thetas[0] = M_PI*0.5; thetas[Nsg] = M_PI*2.5; thetas[Np0] = M_PI*2.5;
    thetas[Nsg/2] = M_PI*1.5;
    for(int j=1; j<=(Nsg-2)/4; j++)
    {
        double dt = deltaTheta;
        thetas[j] = thetas[j-1] + dt;
        thetas[Nsg/2 - j] = thetas[Nsg/2 - j + 1] - dt;
        thetas[Nsg/2 + j] = thetas[Nsg/2 + j -1] + dt;
        thetas[Nsg - j] = thetas[Nsg - j + 1] - dt;
    }
    refinethetas(Nsg, Np0 - Nsg, 2.*M_PI, M_PI/5., thetas);

    double radius = 0., theta;
    //calculate mesh numbers
    layers = 0;
    double rs[2000];
    radius = rc; rs[layers] = radius;
    while(radius<R0)
    {
        layers++;
        double finalH = radius*deltaTheta2;
        radius += meshh(radius-rc, firstH, finalH); rs[layers] = radius;
    }
    layers+=2;//the last two layers is divied into three
    rs[layers-1] = R0;
    rs[layers-2] = (rs[layers-1] - rs[layers-4])*2./3. + rs[layers-4];
    rs[layers-3] = (rs[layers-1] - rs[layers-4])*1./3. + rs[layers-4];
    //generate points
    x = new double*[layers]; y = new double*[layers];
    for(int i=0; i<layers; i++)
    {
        x[i] = new double[Np0];
        y[i] = new double[Np0];
    }
    for(int i=0; i<layers; i++)
    {
	    for(int j=0; j<Np0; j++)
	    {
	        theta = thetas[j];
	        x[i][j] = rs[i]*cos(theta);
	        y[i][j] = rs[i]*sin(theta);
	    }
    }
    //curved edges
    ncedge = 0;
    cc = new double**[10000];
    cedges = new int[10000];
    for(int i=0; i<10000; i++)
    {
        cc[i] = new double*[curvedpoints];
        for(int j=0; j<curvedpoints; j++) cc[i][j] = new double[3];
    }
    for(int i=0; i<10; i++)
    {
        for(int j=0; j<Np0; j++)
        {
            cedges[ncedge] = 2*i*Np0+j;
            double s0 = thetas[j], s1 = thetas[j+1];
            double st = (s1-s0)/(curvedpoints-1.);
            for(int k=0; k<curvedpoints; k++)
            {
                theta = s0 + st*k;
                cc[ncedge][k][0] = rs[i]*cos(theta);
                cc[ncedge][k][1] = rs[i]*sin(theta);
                cc[ncedge][k][2] = 0.;
            }
            ncedge++;
        }
    }
    for(int j=0; j<Np0; j++)
    {
        cedges[ncedge] = 2*(layers-1)*Np0+j;
        double s0 = thetas[j], s1 = thetas[j+1];
        double st = (s1-s0)/(curvedpoints-1.);
        for(int k=0; k<curvedpoints; k++)
        {
            theta = s0 + st*k;
            cc[ncedge][k][0] = R0*cos(theta);
            cc[ncedge][k][1] = R0*sin(theta);
            cc[ncedge][k][2] = 0.;
        }
        ncedge++;
    }
    return 0;
}
int ellbl_pointsgen(double ella, double ellb, int Np0, int Nsg, int Ns0, double R0, double firstH, int curvedpoints,
                    double &deltaTheta, int &layers, double** &x, double** &y, int &ncedge, int* &cedges, double *** &cc)
{//Np0 total circumferential mesh; Nsg total number used for smooth mesh; Ns0*4 is the real smooth mesh number
    double ellratio = ellb/ella;
    double theta, mu, mu0, *thetas, firstmu;
    thetas = new double[Np0+1];
    int npl = 2+2+2, npr = 12+2;
    int Npt = Np0 - npr - npl;
    double ds = 2.*M_PI/Nsg, da = (M_PI_2 - ds*Ns0)/(Npt - 4.*Ns0)*4.;//smooth, angle
    deltaTheta = ds;
    mu0 = 0.5*log((1.+ellratio)/(1.-ellratio));
    double alpha = 2.*ella/(exp(mu0)+exp(-mu0));
    firstmu = firstH/ella;
    //double deltaMu = 0.4*deltaTheta;
    double deltaMu = 0.6*deltaTheta;
    //circumferential mesh
    thetas[0] = M_PI*0.5; thetas[Npt] = M_PI*2.5; thetas[Np0] = M_PI*2.5;
    thetas[Npt/2] = M_PI*1.5;
    for(int j=1; j<=(Npt-2)/4; j++)
    {
        double dt;
        if(j<=Ns0)
        {
            dt = ds;
        }
        else
        {
            dt = da;
        }
        thetas[j] = thetas[j-1] + dt;
        thetas[Npt/2 - j] = thetas[Npt/2 - j + 1] - dt;
        thetas[Npt/2 + j] = thetas[Npt/2 + j -1] + dt;
        thetas[Npt - j] = thetas[Npt - j + 1] - dt;
    }
    refinethetas(Npt, npl-2, M_PI, 0.25, thetas); Npt += npl-2;
    refinethetas(Npt, 2, M_PI, 0.1, thetas); Npt += 2;
    refinethetas(Npt, npr-2, 2.*M_PI, M_PI/3.5, thetas); Npt += npr-2;
    refinethetas(Npt, 2, 2.*M_PI, 0.1, thetas); Npt += 2;
    //refinethetas(Npt, Np0 - Npt, 1.5*M_PI, M_PI/3, thetas);
    double radius = 0.;
    //calculate mesh numbers
    layers = 0;
    double mus[2000];
    mu = mu0; mus[layers] = mu;
    while(radius<R0)
    {
        layers++;
        mu += firstmu + (deltaMu - firstmu)*fStretch(mu, mu0); mus[layers] = mu;
        radius = alpha*(exp(mu)+exp(-mu))*0.5;
    }
    layers+=2;//the last two layers is divied into three
    mus[layers-1] = log(R0/alpha + sqrt(R0*R0/(alpha*alpha)-1.));
    mus[layers-2] = (mus[layers-1] - mus[layers-4])*2./3. + mus[layers-4];
    mus[layers-3] = (mus[layers-1] - mus[layers-4])*1./3. + mus[layers-4];
    //generate points
    x = new double*[layers]; y = new double*[layers];
    for(int i=0; i<layers; i++)
    {
        x[i] = new double[Np0];
        y[i] = new double[Np0];
    }
    //points generation
    for(int i=0; i<layers; i++)
    {
    	mu = mus[i];
    	if(layers-1 == i)
    	{
		    for(int j=0; j<Np0; j++)
		    {
		        theta = thetas[j];
		        x[i][j] = R0*cos(theta);
		        y[i][j] = R0*sin(theta);
		    }
    	}
    	else
    	{
		    for(int j=0; j<Np0; j++)
		    {
		        theta = thetas[j];
		        x[i][j] = alpha*0.5*(exp(mu)+exp(-mu))*cos(theta);
		        y[i][j] = alpha*0.5*(exp(mu)-exp(-mu))*sin(theta);
		    }
        }
    }
    //curved edges
    ncedge = 0;
    cc = new double**[10000];
    cedges = new int[10000];
    for(int i=0; i<10000; i++)
    {
        cc[i] = new double*[curvedpoints];
        for(int j=0; j<curvedpoints; j++) cc[i][j] = new double[3];
    }
    for(int i=0; i<10; i++)
    {
    	mu = mus[i];
        for(int j=0; j<Np0; j++)
        {
            if(i>0 && fabs(thetas[j]-M_PI)>=M_PI/7. && fabs(thetas[j]-2.*M_PI)>=M_PI/9.) continue;
            cedges[ncedge] = 2*i*Np0+j;
            double s0 = thetas[j], s1 = thetas[j+1];
            double st = (s1-s0)/(curvedpoints-1.);
            for(int k=0; k<curvedpoints; k++)
            {
                theta = s0 + st*k;
                cc[ncedge][k][0] = alpha*0.5*(exp(mu)+exp(-mu))*cos(theta);
                cc[ncedge][k][1] = alpha*0.5*(exp(mu)-exp(-mu))*sin(theta);
                cc[ncedge][k][2] = 0.;
            }
            ncedge++;
        }
    }
    for(int j=0; j<Np0; j++)
    {
        cedges[ncedge] = 2*(layers-1)*Np0+j;
        double s0 = thetas[j], s1 = thetas[j+1];
        double st = (s1-s0)/(curvedpoints-1.);
        for(int k=0; k<curvedpoints; k++)
        {
            theta = s0 + st*k;
            cc[ncedge][k][0] = R0*cos(theta);
            cc[ncedge][k][1] = R0*sin(theta);
            cc[ncedge][k][2] = 0.;
        }
        ncedge++;
    }
    return 0;
}
int circle_rect(int Np0, int layer0, double R0, double deltaTheta, double **x0, double **y0, double R1,
                double &downslen, int &layer1, int &iupr, int &idor, double **&x1, double **&y1)
{//generating circle to rectangle points
    double tempm = 10*R1, temp;
    int Nplow, Nphig;
    for(Nplow=0; Nplow<Np0; Nplow++)
    {
        if(y0[layer0-1][Nplow+1]>=y0[layer0-1][Nplow]) break;
    }
    for(Nphig=Nplow; Nphig<Np0; Nphig++)
    {
        if(y0[layer0-1][Nphig+1]<=y0[layer0-1][Nphig]) break;
    }
    tempm = 10*R1;
    for(int i=Np0-1; i>=0; i--)
    {
        if(x0[layer0-1][i]<1.e-6) continue;
        if(y0[layer0-1][i]<1.e-6) continue;
        temp = (R1*x0[layer0-1][i]/y0[layer0-1][i] - downslen);
        if(fabs(temp)<tempm)
        {
            tempm = fabs(temp);
            iupr = i;
        }
        if(temp>0.) break;
    }
    tempm = 10*R1;
    for(int i=Np0-1; i>=0; i--)
    {
        if(x0[layer0-1][i]< 1.e-6) continue;
        if(y0[layer0-1][i]>-1.e-6) continue;
        temp = (-R1*x0[layer0-1][i]/y0[layer0-1][i] - downslen);
        if(fabs(temp)<tempm)
        {
            tempm = fabs(temp);
            idor = i;
        }
        if(temp<0.) break;
    }
    downslen = max(R1*x0[layer0-1][iupr]/y0[layer0-1][iupr], -R1*x0[layer0-1][idor]/y0[layer0-1][idor]);
    layer1 = 0;
    double rs[2000];
    double radius = R0;
    while(radius<R1)
    {
        double finalH = radius*deltaTheta;
        radius += finalH; rs[layer1] = radius;
        layer1++;
    }
    for(int i=0; i<layer1-1; i++)
    {
        rs[i] = (rs[i] - R0)/(rs[layer1-1] - R0);
    }
    rs[layer1-1] = 1.;
    //generate points
    x1 = new double*[layer1];
    y1 = new double*[layer1];
    for(int i=0; i<layer1; i++)
    {
        x1[i] = new double[Np0];
        y1[i] = new double[Np0];
    }
    double rLmax = sqrt(downslen*downslen+R1*R1);
    double rdelta;
    for(int j=0; j<Np0; j++)
    {
        double temprpre = R0;
	    for(int i=0; i<layer1; i++)
	    {
	        double Rt = 0.;
	        if(j<=Nplow || j>=Nphig)
	        {
                Rt = rs[i]*(R1-R0) + R0;
            }
            if((j>Nplow && j<idor) || (j>iupr && j<Nphig))
            {
                double tempr = fabs(R1*R0/y0[layer0-1][j]);
                Rt = rs[i]*(tempr-R0) + R0;
            }
            if(j>idor && j<iupr)
            {
                double tempr = downslen*R0/x0[layer0-1][j];
                rdelta = (tempr - R0)/(rLmax - R0);
                double r2 = min(pow(rs[i],rdelta), (double(i)+1.)/layer1);
                double compr = min(sqrt((tempr-R0)/(R1-R0)), 1.);
                Rt = min(r2*(tempr-R0) + R0, compr*rs[i]*(R1-R0)+R0);
                temprpre = Rt = temprpre + min(Rt - temprpre, (tempr - temprpre)/double(layer1 - i));
            }
            if(j==idor || j==iupr)
            {
                double tempr = max(0.9*sqrt(downslen*downslen+R1*R1), R1);
                Rt = rs[i]*(tempr-R0) + R0;
            }
            //if(j==idor || j==iupr)
            //{
            //    double tempr = downslen*R0/x0[layer0-1][j];
            //    Rt = rs[i]*(tempr-R0) + R0;
            //}
	        x1[i][j] = Rt*x0[layer0-1][j]/R0;
	        y1[i][j] = Rt*y0[layer0-1][j]/R0;

	    }
    }
    return 0;
}
int wake_rect(double downslen, int layers, int idor, int iupr, double **x, double **y,
                int &is3, int &ie3, int &Np2, double &wakeslen, double &ysplus, double &ysminus, double yeplus, double yeminus,
                int &layer2, double &hflows, double &hflowe, double **&x2, double **&y2, double &wakeplus, int &layer2plus)
{
    for(ie3=idor; ie3<iupr; ie3++)
    {
        if(y[layers-1][ie3] >= ysplus) break;
    }
    for(is3=iupr; is3>idor; is3--)
    {
        if(y[layers-1][is3] <= ysminus) break;
    }
    ysplus = y[layers-1][ie3];
    ysminus = y[layers-1][is3];
    if(hflows<1.e-7) hflows = 1.e7;
    hflows = min((y[layers-1][(ie3+is3)/2+3] - y[layers-1][(ie3+is3)/2-2])*0.2,hflows);
    layer2 = int(wakeslen*2./(hflowe+hflows)+0.5);
    double delta = (hflowe - hflows)/(layer2-1);
    wakeslen = layer2*(hflows+hflowe)*0.5;
    layer2plus = int(wakeplus/hflowe+0.5);
    wakeplus = layer2plus*hflowe;
    x2 = new double*[layer2+layer2plus];
    y2 = new double*[layer2+layer2plus];
    for(int i=0; i<layer2+layer2plus; i++)
    {
        x2[i] = new double[ie3-is3+1];
        y2[i] = new double[ie3-is3+1];
    }
    double xs, xe;
    xs = downslen; xe = downslen + wakeslen;
    int is3r, ie3r;
    for(is3r=is3; is3r<ie3; is3r++)
    {
        if(y[layers - 1][is3r+1] - y[layers - 1][is3r]<1.5*hflows) break;
    }
    for(ie3r=ie3; ie3r>is3; ie3r--)
    {
        if(y[layers - 1][ie3r] - y[layers - 1][ie3r-1]<1.5*hflows) break;
    }
    double ysmr = y[layers - 1][is3r], yspr = y[layers - 1][ie3r];
    //double yemr = -(ie3r - is3r)*hflowe*0.30, yepr = (ie3r - is3r)*hflowe*0.30;
    double yemr = ysmr*1.1, yepr = yspr*1.1;
    for(int j=is3; j<=ie3; j++)
    {
        double yept, yemt, yspt, ysmt;
        if(j<is3r)
        {
            ysmt = ysminus; yemt = yeminus;
            yspt = ysmr; yept = yemr;
        }
        else if(j>ie3r)
        {
            ysmt = yspr; yemt = yepr;
            yspt = ysplus; yept = yeplus;
        }
        else
        {
            ysmt = ysmr; yemt = yemr;
            yspt = yspr; yept = yepr;
        }
        double kp = (yept - yspt)/wakeslen; double bp = (yspt*xe - yept*xs)/wakeslen;
        double km = (yemt - ysmt)/wakeslen; double bm = (ysmt*xe - yemt*xs)/wakeslen;
        double kd = kp - km; double bd = bp - bm;
        double x0 = -bd/kd;
        for(int i=0; i<layer2; i++)
        {
            x2[i][j-is3] = xs + (1+i)*hflows + 0.5*i*(i+1)*delta;
            y2[i][j-is3] = (x2[i][j-is3]-x0)/(xs-x0)*(y[layers - 1][j] - ysmt) + (km*x2[i][j-is3]+bm);
        }
        for(int i=layer2; i<layer2+layer2plus; i++)
        {
            x2[i][j-is3] = x2[layer2-1][j-is3] + hflowe*(i-layer2+1);
            y2[i][j-is3] = y2[layer2-1][j-is3];
        }
    }
    Np2 = ie3 - is3 + 1;
    return 0;
}
int FEgen_cir(int Np0, int layers, double **x, double **y, double **point, int nspoint, int **edge, int nsedge, int **cell)
{
    //
    //point
    //
    int pointIndex = 0;
    for(int i=0; i<layers; i++)
    {
        for(int j=0; j<Np0; j++)
        {
            point[pointIndex][0] = x[i][j];
            point[pointIndex][1] = y[i][j];
            pointIndex++;
        }
    }
    //
    //edge, cell
    //
    int edgeIndex = 0;
    for(int i=0; i<2*layers-1; i++)
    {
    	if(0 == i%2)
			for(int j=0; j<Np0; j++)
			{
				edge[edgeIndex][0] = nspoint + Np0*(i/2) + j%Np0;
				edge[edgeIndex][1] = nspoint + Np0*(i/2) + (j+1)%Np0;
				if(2*layers-2 != i) cell[Np0*(i/2)+j][0] = nsedge + edgeIndex;
				if(i/2 != 0)        cell[Np0*(i/2-1)+j][2] = nsedge + edgeIndex;
				edgeIndex++;
			}
    	else
    		for(int j=0; j<Np0; j++)
			{
				edge[edgeIndex][0] = nspoint + Np0*(i/2) + j%Np0;
				edge[edgeIndex][1] = nspoint + Np0*(i/2+1) + j%Np0;
				cell[Np0*(i/2)+j][1] = nsedge + edgeIndex;
				cell[Np0*(i/2)+(j-1+Np0)%Np0][3] = nsedge + edgeIndex;
				edgeIndex++;
			}
    }
    return 0;
}
int FEgen_rec(int Np0, int layers, int is3, int ie3, int Np2, int layer2, double **x2, double **y2, double **point, int **edge, int **cell)
{
    int pointIndex = 0;
    for(int i=0; i<layer2; i++)
    {
        for(int j=0; j<Np2; j++)
        {
            point[pointIndex][0] = x2[i][j];
            point[pointIndex][1] = y2[i][j];
            pointIndex++;
        }
    }
    //
    //edge, cell
    //
    int edgeIndex = 0;
    int nspoint = layers*Np0;
    int nsedge = (2*layers-1)*Np0;
    for(int j=0; j<Np2-1; j++) (cell[j][0] = Np0*(2*layers-2) + is3 + j);
    for(int j=0; j<Np2; j++)
    {
        edge[edgeIndex][0] = Np0*(layers-1) + is3 + j;
        edge[edgeIndex][1] = nspoint + j;
        if(j != Np2-1) cell[j][1] = nsedge + edgeIndex;
        if(j != 0)     cell[j -1][3] = nsedge + edgeIndex;
        edgeIndex++;
    }
    for(int i=0; i<layer2; i++)
    {
        if(i>0)
        for(int j=0; j<Np2; j++)
        {
            edge[edgeIndex][0] = nspoint + (i-1)*Np2 + j;
            edge[edgeIndex][1] = nspoint + i*Np2 + j;
            if(j != Np2-1) cell[i*(Np2-1) + j][1] = nsedge + edgeIndex;
            if(j != 0)     cell[i*(Np2-1) + j -1][3] = nsedge + edgeIndex;
            edgeIndex++;
        }
        for(int j=0; j<Np2-1; j++)
        {
            edge[edgeIndex][0] = nspoint + i*Np2 + j;
            edge[edgeIndex][1] = nspoint + i*Np2 + j + 1;
            cell[i*(Np2-1) + j][2] = nsedge + edgeIndex;
            if(i != layer2 - 1) cell[(i+1)*(Np2-1) + j][0] = nsedge + edgeIndex;
            edgeIndex++;
        }
    }
    return 0;
}
int outtecplot_cir(int Np0, int layers, double **x, double **y, ofstream &outtec)
{
    //tecplot head
    outtec << "title = mesh\n";
    outtec << "variables = x, y\n";
    outtec << "zone i = " << Np0+1 << " j = " << layers << endl;
    for(int i=0; i<layers; i++)
    {
	    for(int j=0; j<=Np0; j++)
	    {
	        outtec << fixed << setprecision(7) << setw(20) << x[i][j%Np0] << " " << setw(20) << y[i][j%Np0] << endl;
	    }
    }
    return 0;
}
int outtecplot_rec(int Np0, int layers, double **x, double **y, ofstream &outtec)
{
    //tecplot head
    outtec << "title = mesh\n";
    outtec << "variables = x, y\n";
    outtec << "zone i = " << Np0+1 << " j = " << layers << endl;
    for(int i=0; i<layers; i++)
    {
	    for(int j=0; j<Np0; j++)
	    {
	        outtec << fixed << setprecision(7) << setw(20) << x[i][j] << " " << setw(20) << y[i][j] << endl;
	    }
    }
    return 0;
}
int outFE(int npoint, double **point, int nedge, int **edge, int curvedpoints, int ncedge, int *cedge, double ***cc, int ncell, int **cell, ofstream &outxml)
{
    //xml head
    outxml << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << "\n";
    outxml << "<NEKTAR>" << "\n";
    outxml << "    <GEOMETRY DIM=\"2\" SPACE=\"2\">" << "\n";
    outxml << "        <VERTEX>" << "\n";
    //
    //vertex
    //
    for(int i=0; i<npoint; i++)
    {
        outxml << "            <V ID=\"" << i << "\">" ;
        outxml << scientific << setprecision(17);
        outxml << setw(26) << point[i][0] << setw(26) << point[i][1] << setw(26) << 0.;
        outxml << "</V>" << "\n";
    }
    outxml << "        </VERTEX>" << "\n";
    outxml << "        <EDGE>" << "\n";
    //
    //edge
    //
    for(int i=0; i<nedge; i++)
    {
    	outxml << "            <E ID=\"" << i << "\">";
        outxml << setw(7) << edge[i][0] << setw(7) << edge[i][1] << "</E>" << "\n";
    }
    outxml << "        </EDGE>" << "\n";
    //
    //elements
    //
    outxml << "        <ELEMENT>" << "\n";
    for(int i=0; i<ncell; i++)
    {
        outxml << "            <Q ID=\"" << i << "\">"
            << setw(7) << cell[i][0] << setw(7) << cell[i][1]
            << setw(7) << cell[i][2] << setw(7) << cell[i][3] << " </Q>" << "\n";
    }
    outxml << "        </ELEMENT>" << "\n";
    //
    //curved edges
    //
    outxml << "        <CURVED>" << "\n";
    for(int i=0; i<ncedge; i++)
    {
        outxml << "            <E ID=\"" << i << "\" EDGEID=\"" << cedge[i]
        << "\" NUMPOINTS=\"" << curvedpoints << "\" TYPE=\"PolyEvenlySpaced\">";
        for(int k=0; k<curvedpoints; k++)
        {
            outxml << setw(26) << cc[i][k][0] << setw(26) << cc[i][k][1] << setw(26) << cc[i][k][2];
        }
        outxml << "</E>" << "\n";
    }
    outxml << "        </CURVED>" << "\n";
    return 0;
}
int outCOMPO(int Np0, int layers, int is3, int ie3, int Np2, int layer2, int layer2plus, int **cell, ofstream &outxml)
{
    //
    //composite
    //
    outxml << "        <COMPOSITE>" << "\n";
    //solid boundary
    outxml << "            <C ID=\"1\"> E[";
    for(int j=0; j<Np0-1; j++) outxml << cell[j][0] << ",";
    outxml << cell[Np0-1][0] << "] </C>" << "\n";
    //side boundary
    outxml << "            <C ID=\"2\"> E[";
    for(int j=0; j<Np0; j++)
    {
        if(j>=is3 && j<ie3) continue;
        outxml << cell[(layers-2)*Np0+j][2] << ",";
    }
    for(int i=0; i<layer2; i++)
    {
        outxml << cell[(layers-1)*Np0+i*(Np2-1)][1] << ",";
    }
    for(int i=0; i<layer2-1; i++)
    {
        outxml << cell[(layers-1)*Np0+(1+i)*(Np2-1)-1][3] << ",";
    }
    outxml << cell[(layers-1)*Np0+layer2*(Np2-1)-1][3] << "] </C>" << "\n";
    //rear side boundary
    outxml << "            <C ID=\"3\"> E[";
    for(int i=layer2; i<layer2+layer2plus; i++)
    {
        outxml << cell[(layers-1)*Np0+i*(Np2-1)][1] << ",";
    }
    for(int i=layer2; i<layer2+layer2plus-1; i++)
    {
        outxml << cell[(layers-1)*Np0+(1+i)*(Np2-1)-1][3] << ",";
    }
    outxml << cell[(layers-1)*Np0+(layer2+layer2plus)*(Np2-1)-1][3] << "] </C>" << "\n";
    //outflow boundary
    outxml << "            <C ID=\"4\"> E[";
    for(int i=0; i<Np2-2; i++)
    {
        outxml << cell[(layers-1)*Np0+(layer2+layer2plus-1)*(Np2-1) + i][2] << ",";
    }
    outxml << cell[(layers-1)*Np0+(layer2+layer2plus)*(Np2-1)-1][2] << "] </C>" << "\n";
    //the first layer cells near wall
    outxml << "            <C ID=\"5\"> Q[" << 0 << "-" << (layers-1)*Np0 + (layer2+layer2plus)*(Np2-1)-1 << "] </C>" << "\n";
    outxml << "        </COMPOSITE>" << "\n";
    outxml << "        <DOMAIN> C[5] </DOMAIN>" << "\n";
    outxml << "    </GEOMETRY>" << "\n";
    //expansion
    outxml << "    <EXPANSIONS>" << "\n";
    //outxml << "        <E COMPOSITE=\"C[5]\" NUMMODES=\"11\" TYPE=\"MODIFIED\" FIELDS=\"u,v,p\" />" << "\n";
    outxml << "        <E COMPOSITE=\"C[5]\" NUMMODES=\"10, 10\" BASISTYPE=\"Modified_A,Modified_A\" POINTSTYPE=\"GaussLobattoLegendre,GaussLobattoLegendre\" NUMPOINTS=\"12,12\" FIELDS=\"u,v,p\" />" << "\n";
    outxml << "    </EXPANSIONS>" << "\n";
    outxml << "</NEKTAR>" << "\n";
    return 0;
}
int coord_trans(int Np0, int layers, double alpha, double **x, double **y)
{//alpha less than M_PI
    double tx, ty;
    double ca = cos(alpha), sa = sin(alpha);
    for(int i=0; i<layers; i++)
    {
        for(int j=0; j<Np0; j++)
        {
            tx = x[i][j]; ty = y[i][j];
            x[i][j] = ca*tx - sa*ty;
            y[i][j] = sa*tx + ca*ty;
        }
    }
    return 0;
}
int cedge_trans(int ncedge, int curvedpoints, double alpha, double ***cc)
{//alpha less than M_PI
    double tx, ty;
    double ca = cos(alpha), sa = sin(alpha);
    for(int i=0; i<ncedge; i++)
    {
        for(int j=0; j<curvedpoints; j++)
        {
            tx = cc[i][j][0]; ty = cc[i][j][1];
            cc[i][j][0] = ca*tx - sa*ty;
            cc[i][j][1] = sa*tx + ca*ty;
        }
    }
    return 0;
}
