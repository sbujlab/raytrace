#include <TMinuit.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>

static double B[50][50][101][3];
static double expang;
static double step;
static double scale;
static double rmin[3];
static double rmax[3];
static int n[3];

static int    nlast;
static double lastray[2][200000];

static double mag_zoff;

bool testcol( double x, double y ){
    // Tyler's drawing Feb 9 ELOG 3361

    double crmin = 0.2125;
    double crmax = crmin + 0.09;
    double h = 0.0935;
    double chamf = 0.04;

    if( fabs(x) > h ) return false;
    if( y < -chamf ) return false;

    double azi = atan2(x,y+crmin);
//    printf("%f %f azi = %f\n", x,y, azi*180/3.14159);

    // Guess
    if( fabs(azi) > 0.4 ) return false;

    double r  =  sqrt( x*x + (y + 0.03 + crmin)*(y + 0.03 + crmin) );

    if( r > crmax ) return false;
    if( r < crmin ) return false;

//    printf("In acc\n");
    return true;
}

void getB( double rlab[3], double interpB[3]){
    double r[3];

    r[0] = rlab[0];
    r[1] = rlab[1];
    r[2] = rlab[2]-(1.7538 + mag_zoff);
//    r[2] = rlab[2]-1.8088;

    double x[3], res[3];
    int idx[3];

    interpB[0] = 0.0; interpB[1] = 0.0; interpB[2] = 0.0;

    int i;
    for(i = 0; i < 3; i++ ){
        x[i] = (r[i]-rmin[i])/(rmax[i]-rmin[i]);

        idx[i] = floor(x[i]*(n[i]-1));
        res[i] = x[i]*(n[i]-1) - idx[i];

        if( idx[i] < 0 || idx[i] >= n[i] ){
            return;
        }
    }

    double c00, c10, c01, c11;
    double c0, c1;

    for(i = 0; i < 3; i++ ){
        c00 = B[idx[0]][idx[1]][idx[2]][i]*(1.0-res[0]) + B[idx[0]+1][idx[1]][idx[2]][i]*res[0];
        c10 = B[idx[0]][idx[1]][idx[2]+1][i]*(1.0-res[0]) + B[idx[0]][idx[1]+1][idx[2]][i]*res[0];
        c01 = B[idx[0]][idx[1]][idx[2]+1][i]*(1.0-res[0]) + B[idx[0]+1][idx[1]][idx[2]+1][i]*res[0];
        c11 = B[idx[0]][idx[1]+1][idx[2]+1][i]*(1.0-res[0]) + B[idx[0]+1][idx[1]+1][idx[2]+1][i]*res[0];

        c0  = c00*(1.0-res[1]) + c10*res[1];
        c1  = c01*(1.0-res[1]) + c11*res[1];

        interpB[i] = scale*(  c0*(1.0-res[2]) + c1*res[2] );
    }

//    printf("B %f %f %f -> %f %f %f\n", r[0], r[1], r[2], interpB[0], interpB[1], interpB[2] );

    return;
}

void getF( double r[3], double v[3], double F[3] ){
    double B[3];
    getB(r,B);

    F[0] = v[1]*B[2] - v[2]*B[1];
    F[1] = v[2]*B[0] - v[0]*B[2];
    F[2] = v[0]*B[1] - v[1]*B[0];

//    printf("%g %g %g\n", F[0], F[1], F[2]);

    return;
}

void prop( double r0[3], double v0[3], double r1[3], double v1[3] ){
    double gamma = 1.063/0.000511;
    double m = 9.1e-31*gamma/1.6e-19;
    double K1[3], K2[3], K3[3], K4[3];

    getF( r0, v0, K1 );

    double tempx[3], tempv[3];
    int i;
    for( i = 0; i < 3; i++){
        tempx[i] = r0[i] + step*v0[i]/2.0 + step*step*K1[i]/8.0/m;
        tempv[i] = v0[i] + step*K1[i]/2.0/m;
    }

    getF( tempx, tempv, K2 );

    for( i = 0; i < 3; i++){
        tempx[i] = r0[i] + step*v0[i]/2.0 + step*step*K1[i]/8.0/m;
        tempv[i] = v0[i] + step*K2[i]/2.0/m;
    }
    
    getF( tempx, tempv, K3 );

    for( i = 0; i < 3; i++){
        tempx[i] = r0[i] + step*v0[i] + step*step*K3[i]/2.0/m;
        tempv[i] = v0[i] + step*K3[i]/m;
    }

    getF( tempx, tempv, K4 );

//    printf("%g %g %g %g\n", K1[0], K2[0], K3[0], K4[0] );

    for( i = 0; i < 3; i++){
        r1[i] = r0[i] + step*v0[i] + step*step*(K1[i] + K2[i] + K3[i])/6.0/m;
        v1[i] = v0[i] + step*(K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i])/6.0/m;
    }

//    printf("%f %f %f -> %f %f %f\n", r0[0], r0[1], r0[2], v0[0], v0[1], v0[2] );
//    printf("%f %g -> %g (%g)\n", r0[2], v0[2], r0[2] + step*v0[2], step );

    return;
}

// Test if it's on the wrong side of the quad boundary
bool testplane( double r1[3] ){
    double q1ent[3];
    double q1norm[3];

    q1ent[0] = 0.29847;
    q1ent[2] = 2.39597;
//    q1ent[2] = 2.45097;

    q1norm[0] = sin(12.5*3.14159/180.);
    q1norm[2] = cos(12.5*3.14159/180.);

    double dot = (r1[0]-q1ent[0])*q1norm[0] + (r1[2]-q1ent[2])*q1norm[2];

    if( dot > 0.0 ){
        return true;
    }

    return false;
}

double getqy( double r[3] ){
    double q1ent[3];
    q1ent[0] = 0.29847;
    q1ent[1] = 0.0;
    q1ent[2] = 2.39597;
//    q1ent[2] = 2.45097;

    double rtrans[3], rrot[3];
    int i;

    // Translate
    for( i = 0; i < 3; i++ ){
        rtrans[i] = r[i]-q1ent[i];
    }
    // Rotate
    rrot[0] = rtrans[0]*cos(12.5*3.14159/180) - rtrans[2]*sin(12.5*3.14159/180);
    rrot[1] = rtrans[1];
    rrot[2] = rtrans[2]*cos(12.5*3.14159/180) + rtrans[0]*sin(12.5*3.14159/180);

    return rrot[0];
}

double getqphi( double v[3] ){
    double vrot[3];
    // Rotate
    vrot[0] = v[0]*cos(12.5*3.14159/180) - v[2]*sin(12.5*3.14159/180);
    vrot[1] = v[1];
    vrot[2] = v[2]*cos(12.5*3.14159/180) + v[0]*sin(12.5*3.14159/180);

    return vrot[0]/vrot[2];
}


double getqx( double r[3] ){
    double q1ent[3];
    q1ent[0] = 0.29847;
    q1ent[1] = 0.0;
    q1ent[2] = 2.39597;
//    q1ent[2] = 2.45097;

    double rtrans[3], rrot[3];
    int i;

    // Translate
    for( i = 0; i < 3; i++ ){
        rtrans[i] = r[i]-q1ent[i];
    }
    // Rotate
    rrot[0] = rtrans[0]*cos(12.5*3.14159/180) - rtrans[2]*sin(12.5*3.14159/180);
    rrot[1] = rtrans[1];
    rrot[2] = rtrans[2]*cos(12.5*3.14159/180) + rtrans[0]*sin(12.5*3.14159/180);

    return rrot[1];
}

double getqtheta( double v[3] ){
    double vrot[3];
    // Rotate
    vrot[0] = v[0]*cos(12.5*3.14159/180) - v[2]*sin(12.5*3.14159/180);
    vrot[1] = v[1];
    vrot[2] = v[2]*cos(12.5*3.14159/180) + v[0]*sin(12.5*3.14159/180);

    return vrot[1]/vrot[2];
}


//  Go from z at angle a to quad, get y and phi

void toquad( double z, double ang, double outang, double rasx, double rasy, double &yq1, double &phq1, double &xq1, double &thq1 ){
    double r[3];
    r[0] = rasx;
    r[1] = rasy;
    r[2] = z;
    double v[3];
    v[0] = sin(ang*3.14159/180)*cos(outang*3.14159/180.0)*3e8;
    v[1] = sin(ang*3.14159/180)*sin(outang*3.14159/180.0)*3e8;
    v[2] = cos(ang*3.14149/180)*3e8;

    double r1[3], v1[3];

    int i;
    int nstep = 0;
    nlast = 0;
    while( !testplane(r) && nstep < 200000 ){
        lastray[0][nlast] = r[0];
        lastray[1][nlast] = r[2];
        prop( r,v,r1,v1 );
        for( i = 0; i < 3; i++ ){
            r[i] = r1[i];
            v[i] = v1[i];
        }
//        printf("prop step %d,  (%f %f %f) (%f %f %f)\n", nstep, r[0], r[1], r[2], v[0]/3e8, v[1]/3e8, v[2]/3e8 );
        nstep++;
        nlast++;
    }
    lastray[0][nlast] = r[0];
    lastray[1][nlast] = r[2];
    nlast++;

    yq1  = getqy( r );
    phq1 = getqphi( v );

    xq1  = getqx( r );
    thq1 = getqtheta( v );

    return;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double y, ph, x, th;
    scale = par[0];
    toquad(par[1], expang, 0.0, 0.0, 0.0, y, ph, x, th);

    double dfp = 0.003;
    double dy  = 0.002;
    double dph = 0.001;

//    f = y*y/(0.002*0.002) + ph*ph/(0.001*0.001);
    f = pow(y*2.54 - 5.56*ph,2.0)/(dfp*dfp) + y*y/(dy*dy) + ph*ph/(dph*dph);

    printf("Error components fp %f y %f ph %f\n", (y*2.54 - 5.56*ph)/dfp, y/dy, ph/dph);

    printf("z = %f, s = %f -> %f %f\n", par[1], scale, y, ph);
}


int opt_pos(){
    step  = 1e-12;
    scale = 1.0;

    // Read in field
    
    rmin[0] = -0.400;
    rmin[1] = -0.140;
    rmin[2] = -1.000;
    rmax[0] = 0.400;
    rmax[1] = 0.140;
    rmax[2] = 1.000;


    n[0] = 41;
    n[1] = 15;
    n[2] = 41;

    FILE *f = fopen("prex_septumfield.dat", "r");

    int nread = 6;

    int i, j, k, l;

    double dummy;

    for( j = 0; j < n[1]; j++ ){
        for( k = 0; k < n[0]; k++ ){
            for( i = 0; i < n[2]; i++ ){
                fscanf(f, "%lf%lf%lf%lf%lf%lf", &dummy, &dummy, &dummy, &B[i][j][k][0], &B[i][j][k][2], &B[i][j][k][1]);
    //            printf("%f %f %f\n", B[i][j][k][0],  B[i][j][k][1],  B[i][j][k][2] );
            }
        }
    }

    fclose(f);
    printf("B read \n");


    /*
    double y, ph;
    scale = 1.2;
    toquad(-0.45, 4, y, ph);

    printf("y = %f, ph = %f\n", y, ph);

    return;
    */
    double y,ph,zfit, x, th;

    TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 5 params

    gMinuit->SetFCN(fcn);

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    mag_zoff = 0.00;
    double beamE = 1.068;
    expang = 5.0;
    Double_t vstart[2] = {-1.18,  0.00};
    Double_t vstep[2] = {0.1, 0.00};



    gMinuit->mnparm(0, "scale", vstart[0], vstep[0], 0,0,ierflg);
    gMinuit->mnparm(1, "z", vstart[1], vstep[1], 0,0,ierflg);

    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);

    gMinuit->GetParameter(0, scale, dummy);
    gMinuit->GetParameter(1, zfit, dummy);

    toquad( zfit, expang, 0.0, 0.0, 0.0, y, ph, x, th);

    printf("y = %f ph = %f\n", y, ph);
    printf("y_{fp} = %f\n", 2.56*y - 5.54*ph);

  //  gStyle->SetOptStat(0);
   
    TH2F *frm = new TH2F("frm", "Ray Trace", 100, -0.5, 3.0, 100, -0.5, 3.0);
    TGraph *g = new TGraph(nlast, lastray[0], lastray[1]);
    g->SetLineColor(kRed);
    
    TCanvas *c = new TCanvas("c", "Ray Trace", 1024, 1024);

    TLine *l1 = new TLine(0.06, 1.7538 - 0.5 + mag_zoff, 0.42, 1.7538 - 0.5 + mag_zoff );
    TLine *l2 = new TLine(0.06, 1.7538 + 0.5 + mag_zoff, 0.42, 1.7538 + 0.5 + mag_zoff );
    TLine *l3 = new TLine(0.42, 1.7538 - 0.5 + mag_zoff, 0.42, 1.7538 + 0.5 + mag_zoff );
    TLine *l4 = new TLine(0.06, 1.7538 - 0.5 + mag_zoff, 0.06, 1.7538 + 0.5 + mag_zoff );

    double q1l = 0.12;
    TLine *q1 = new TLine( 0.29847 - 0.5*q1l*cos(12.5*3.14159/180) , 2.39597 + 0.5*q1l*sin(12.5*3.14159/180) , 0.29847 + 0.5*q1l*cos(12.5*3.14159/180) , 2.39597 - 0.5*q1l*sin(12.5*3.14159/180)  );

    frm->Draw();
    g->Draw("L");
    l1->Draw();
    l2->Draw();
    l3->Draw();
    l4->Draw();
    q1->Draw();

    TCanvas *cB = new TCanvas("cB", "Field");
    double r[3], vB[3];
    double By[2700];
    double By_center[2700];
    double z[2700];
    for( i = 0; i < 2700; i++ ){
        z[i] = i*1e-3;

        r[0] = z[i]*sin(5.0*3.14159/180);
        r[1] = 0.0;
        r[2] = z[i]*cos(5.0*3.14159/180);

        getB(r, vB);

        By[i] = vB[1];


        r[0] = 0.2;
        r[1] = 0.0;
        r[2] = z[i];
        getB(r, vB);
        By_center[i] = vB[1];

    }
    TGraph *gB = new TGraph(2700, z, By);
    TGraph *gBc = new TGraph(2700, z, By_center);
    gB->SetTitle("B_{y} on path");
    gB->GetXaxis()->SetTitle("z [m]");
    gB->GetXaxis()->CenterTitle();
    gB->GetYaxis()->SetTitle("B_{y} [T]");
    gB->GetYaxis()->CenterTitle();
    gB->Draw("AP");

    gBc->SetLineStyle(3);
    gBc->SetLineColor(kBlue);
//    gBc->Draw("L");


    TH1F *acc = new TH1F("acc", "acc", 100, 3.0, 7.0 );
    TH1F *hq2= new TH1F("q2", "q2", 100, 0.003, 0.017 );
    TH2F *hq1= new TH2F("q1", "q1", 100, -0.12, 0.12, 100, -0.12, 0.12);

    int nevt = 100000;

    printf("Simulating %d events\n", nevt);
    for( i = 0; i < nevt; i++ ){
        if( i % 1000 == 0 ){
            printf("\tEvent %8d\n", i);
        }
        double thlabrad = acos(gRandom->Uniform(cos(3.0*3.14159/180), cos(7.0*3.14159/180)));
        double thlab = thlabrad*180/3.14159;
        double phlab = gRandom->Uniform(-25.0,25.0);
        double rasx  = gRandom->Uniform(-0.002, 0.002);
        double rasy  = gRandom->Uniform(-0.002, 0.002);
//
//        double rasx  = 0.0;
//        double rasy  = 0.0;

        double q2 = 2.0*beamE*beamE*(1.0 - cos(thlabrad));

        double mott = pow(cos(thlabrad/2.0)/pow(sin(thlabrad/2.0),2.0),2.0)*(1.0 + q2*pow(tan(thlabrad/2.0)/(208*0.931),2.0)/2.0);

        toquad( zfit, thlab, phlab, rasx, rasy, y, ph, x, th);

//        printf("%f %f -> %f %f %f %f\n", thlab, phlab, x, y, th, ph);
        if(testcol(x,y)){
            acc->Fill(thlab, mott);
            hq2->Fill(q2, mott);
            hq1->Fill(y,x);
        }
    }

    TCanvas *cq2 = new TCanvas();
    hq2->Draw();
    new TCanvas();
    acc->Draw();
    new TCanvas();
    hq1->Draw();


    return 0;
}



