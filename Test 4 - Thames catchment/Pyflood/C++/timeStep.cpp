// TO DO LIST:
// SHOULD RUN ALL LOOPS WITH j FIRST TO IMPROVE CACHING?
// Split into sensible sub functions
// Check edges for drying
// Try initialisation/run split
// Fail gracefully when NaN/very high depth produced



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
//#include <omp.h>
using namespace std;


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

double totalPowerCompTime=0.0;


class array2d
{
	public:
		double getVal(int i, int j);
		void setVal(int i, int j, double val);
		void assign(int xszArg, int yszArg, double *dataArg);
           void zero(int xszArg, int yszArg);
          double sum();
          int cltz();
          double min();
        inline double &operator()(const int &i,const int &j);
        inline double element(int i, int j);
        friend void swapArrays(array2d &a, array2d &b),copyArray(array2d a, array2d b);

//private:
        double *data;
        int xsz, ysz;

};

inline void array2d::zero(int xszArg, int yszArg)
{
    xsz=xszArg;
    ysz=yszArg;

    data=(double *)malloc(xsz*ysz*sizeof(double));

    for(int i=0;i<xsz;i++)
        for(int j=0;j<ysz;j++)
            data[i*ysz+j]=0.;

    return;
}

inline void copyArray(array2d src, array2d dst)
{
    memcpy(dst.data,src.data,sizeof(double)*src.xsz*src.ysz);
    return;
}


inline void swapArrays(array2d &a, array2d &b)
{
    // Assume they're of equal size - no error checking here
    double *tmp=a.data;
    a.data=b.data;
    b.data=tmp;
}


inline void array2d::assign(int xszArg, int yszArg, double *dataArg)
{
	data=dataArg;
	xsz=xszArg;
	ysz=yszArg;
}

inline double &array2d::operator()(const int &i,const int &j)
{
    return data[i*ysz+j];
}

inline double array2d::sum()
{
    int i,j;
    double s=0.0;

    for(int i=0;i<xsz;i++)
        for(int j=0;j<ysz;j++)
            s+=data[i*ysz+j];

    return s;

}

inline double array2d::min()
{
    int i,j;
    double m=1e20;

    for(int i=0;i<xsz;i++)
        for(int j=0;j<ysz;j++)
            m=MIN(m,data[i*ysz+j]);

    return m;

}

inline int array2d::cltz()
{
    int i,j;
    double count=0;

    for(int i=0;i<xsz;i++)
        for(int j=0;j<ysz;j++)
            if(data[i*ysz+j]<0)
                count++;

    return count;
}

inline double timeElapsed(timeval &t1, timeval &t2)
{
    return (double)(t2.tv_sec-t1.tv_sec)+1e-6*(t2.tv_usec-t1.tv_usec);
}

inline double fastPowM43(const double &arg)
{
    float x0=arg;

    x0=x0*x0*x0*x0;

   union {int ix; float x;};

   x = x0;                      // x can be viewed as int.
   ix = 0x2a51067f + ix/3;      // Initial guess.

   return double(1./x);
}



inline double pow43(const double &arg)
{
    timeval time1,time2;
    double it1, it2;

    return 1./cbrt(arg*arg*arg*arg);

}


void xFlowCalc();
void yFlowCalc();
void hEdgeCalc();
double flowOutCalc();


array2d u, v, h, uNew, vNew, hNew, z, hEdgeX, hEdgeY, nGrid;
double depthThresh=1e-6;
const double g=9.81;
bool newAdvectionMethod=true;
double dt;
double dx,beta,crossTermUpwind;
int xsz, ysz, frictionLaw;
double Qout, integratedQout;







extern "C" double totaliser(double *arrArg,int xsz,int ysz)
{
	int i,j;
	double sum=0.;
	array2d arr;

	arr.assign(xsz,ysz,arrArg);

	for(i=0;i<xsz;i++) for (j=0;j<ysz;j++) sum+=arr(i,j);

	return sum;
}

//#define debug

#define profile

extern "C" void timeStep(double *uArg, double *vArg, double *hArg,
    double *uNewArg, double *vNewArg, double *hNewArg,
	double *zArg,
	double *ubcVals, int *ubcLocX, int *ubcLocY, int ubcN,
	double *vbcVals, int *vbcLocX, int *vbcLocY, int vbcN,
	double *hbcVals, int *hbcLocX, int *hbcLocY, int hbcN,
	double *sbcVals, int *sbcLocX, int *sbcLocY, int sbcN,
    int xszArg, int yszArg, double totalTime, double displayInterval,
    double (*calcTimeStep)(), void (*updateBCs)(double,double), void (*rainfall)(double,double),
    void (*display)(int,double,double,double),
    bool (*covergenceTest)(double),
    double dxArg, double *nGridArg, int frictionLawArg)
{
// Set globals by arguments
    dx=dxArg;
//    n=nArg;
//    gnn=g*n*n;
    frictionLaw=frictionLawArg;

// Other "local" variables
    int i, j, ni, nj,ts;

    double dhdx, dhdy, dudx, dudy, dvdy, dvdx, hEdge, u4, v4, q1, q2, q3, q4, dH, dryFactor;
    double wlij, wlNeighbourMax;

    double currentTime=0.,lastDisplayTime=0.,divdt;

    clockid_t tmpTime;
    double totalCallbackTime=0., totalFlowCompTime=0.;

    # ifdef profile
    double pt1, pt2, pt3, pt4, pt5;
    # endif


    xsz=xszArg;
    ysz=yszArg;

    double dx2=dx*dx;
    double dtdx2=dt/dx2;

// Set up array2d types
	u.assign(xsz+1,ysz,uArg);
	v.assign(xsz,ysz+1,vArg);
	h.assign(xsz,ysz,hArg);

	uNew.assign(xsz+1,ysz,uNewArg);
    vNew.assign(xsz,ysz+1,vNewArg);
	hNew.assign(xsz,ysz,hNewArg);

	z.assign(xsz,ysz,zArg);

    nGrid.assign(xsz,ysz,nGridArg);

	hEdgeX.zero(xsz+1,ysz);
	hEdgeY.zero(xsz,ysz+1);


    timeval time1,time2;


    display(0,0.0,0.0,0.0);

    integratedQout=0.0;

    for(ts=0;;ts++)
    {
        # ifdef debug
            printf("\nTick 1: %f %i %f\n",h.sum()*dx*dx, h.cltz(), h.min());
        # endif


        gettimeofday(&time1,NULL);
        dt=calcTimeStep();
        divdt=1./dt;

        updateBCs(currentTime,dt);
        rainfall(currentTime,dt);
        gettimeofday(&time2,NULL);
        totalCallbackTime+=timeElapsed(time1,time2);

    // Insert WL BCs
        for(i=0;i<hbcN;i++) h(hbcLocX[i],hbcLocY[i])=hbcVals[i]-z(hbcLocX[i],hbcLocY[i]);

    // Insert source BCs
        for(i=0;i<sbcN;i++) h(sbcLocX[i],sbcLocY[i])+=sbcVals[i]*dt/dx2;

    // Calculate flows
        gettimeofday(&time1,NULL);

        hEdgeCalc();
        xFlowCalc();
        yFlowCalc();

        gettimeofday(&time2,NULL);
        totalFlowCompTime+=timeElapsed(time1,time2);

    // Insert stuff for flow boundary conditions
        for(i=0;i<ubcN;i++) uNew(ubcLocX[i],ubcLocY[i])=ubcVals[i];
        for(i=0;i<vbcN;i++) vNew(vbcLocX[i],vbcLocY[i])=vbcVals[i];

# ifdef debug
    printf("Tick 2: %f %i %f\n",h.sum()*dx*dx, h.cltz(), h.min());
# endif

    // Check for drying points and scale flows so that depth->0
    // Don't do this for dry cells - these may have flows specified from BCs which we don't want to interfere with
        int nDC;
        for(nDC=0;nDC<10;nDC++)
        {

            int nChanged=0;
//            for(i=1;i<xsz-1;i++) for(j=1;j<ysz-1;j++) //if(h(i,j)>depthThresh) // TO DO: Need to check drying for edge points
            for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) //if(h(i,j)>depthThresh) // TO DO: Need to check drying for edge points
            {
                q1=hEdgeX(i,j)*uNew(i,j);
                q2=hEdgeX(i+1,j)*uNew(i+1,j);
                q3=hEdgeY(i,j)*vNew(i,j);
                q4=hEdgeY(i,j+1)*vNew(i,j+1);

                dH=(q1-q2+q3-q4)*dt/dx;

                if((h(i,j)+dH)<0 ) //-0.001)
                {
                    nChanged++;
                    dryFactor=-h(i,j)/dH;
                    uNew(i,j)*=dryFactor;
                    uNew(i+1,j)*=dryFactor;
                    vNew(i,j)*=dryFactor;
                    vNew(i,j+1)*=dryFactor;

//                    if(uNew(i,j)<0) uNew(i,j)*=dryFactor;
//                    if(uNew(i+1,j)>0) uNew(i+1,j)*=dryFactor;
//                    if(vNew(i,j)<0) vNew(i,j)*=dryFactor;
//                    if(vNew(i,j+1)>0) vNew(i,j+1)*=dryFactor;
                }
            }

            if (nChanged==0) break;

        }

//        printf("nDC=%i\n",nDC);

        if(nDC==10) printf("WARNING: Maximum drying iterations reached\n");

        Qout=flowOutCalc();
        integratedQout+=Qout*dt;


    // Update depths: interior points
    // For cells next to dry cells, need to allow for possibility we've got an internal BC adding water
        double volNeg=0.;

# ifdef debug
    printf("Tick 3: %f %i %f\n",h.sum()*dx*dx, h.cltz(), h.min());
# endif

        for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
        {
            q1=hEdgeX(i,j)*uNew(i,j);
            q2=hEdgeX(i+1,j)*uNew(i+1,j);
            q3=hEdgeY(i,j)*vNew(i,j);
            q4=hEdgeY(i,j+1)*vNew(i,j+1);

//            q3=q4=0.;

//            hNew(i,j)=MAX(h(i,j)+(q1-q2+q3-q4)*dt/dx,0);

            hNew(i,j)=h(i,j)+(q1-q2+q3-q4)*dt/dx;
        }

# ifdef debug
    printf("Tick 4: %f %i %f\n",hNew.sum()*dx*dx, hNew.cltz(), hNew.min());
# endif

    // Wetting points
        for(i=1;i<xsz-1;i++) for(j=1;j<ysz-1;j++) if(hNew(i,j)<depthThresh)
        {
            wlij=z(i,j)+hNew(i,j);
            wlNeighbourMax=-9999;


            if(hNew(i+1,j)>2*depthThresh && (hNew(i+1,j)+z(i+1,j))>wlNeighbourMax)
            {
                wlNeighbourMax=(hNew(i+1,j)+z(i+1,j));
                ni=i+1;
                nj=j;
            }
            else if(hNew(i-1,j)>2*depthThresh && (hNew(i-1,j)+z(i-1,j))>wlNeighbourMax)
            {
                wlNeighbourMax=(hNew(i-1,j)+z(i-1,j));
                ni=i-1;
                nj=j;
            }
            else if(hNew(i,j+1)>2*depthThresh && (hNew(i,j+1)+z(i,j+1))>wlNeighbourMax)
            {
                wlNeighbourMax=(hNew(i,j+1)+z(i,j+1));
                ni=i;
                nj=j+1;
            }
            else if(hNew(i,j-1)>2*depthThresh && (hNew(i,j-1)+z(i,j-1))>wlNeighbourMax)
            {
                wlNeighbourMax=(hNew(i,j-1)+z(i,j-1));
                ni=i;
                nj=j-1;
            }

            if(wlNeighbourMax>-9999 && wlNeighbourMax>wlij) // Wet cell from highest neighbours
            {
                hNew(i,j)=depthThresh;
                hNew(ni,nj)-=(depthThresh-hNew(i,j));
            }
        }

//        printf("Tick 3: %f\n",hNew.sum());
# ifdef debug
    printf("Tick 5: %f %i %f\n",hNew.sum()*dx*dx, hNew.cltz(), hNew.min());
# endif

    // Insert WL BCs
        for(i=0;i<hbcN;i++)
        {
            hNew(hbcLocX[i],hbcLocY[i])=hbcVals[i]-z(hbcLocX[i],hbcLocY[i]);
            if(hNew(hbcLocX[i],hbcLocY[i])<0) hNew(hbcLocX[i],hbcLocY[i])=0.;
        }


//        printf("Tick 4: %f\n",hNew.sum());
# ifdef debug
    printf("Tick 6: %f %i %f\n",hNew.sum()*dx*dx, hNew.cltz(), hNew.min());
# endif

    // Clean up array
        volNeg=0.;
        for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) hNew(i,j)=MAX(hNew(i,j),0);
//        {
//            if(hNew(i,j)<0.0)
//            {
//                printf("Negative depth found at %i,%i\n",i,j);
//                volNeg+=dx*dx*-hNew(i,j);
//                hNew(i,j)=0.;
//            }
//        }

//        printf("Tick 5: %f\n",hNew.sum());

# ifdef debug
    printf("Tick 7: %f %i %f\n",hNew.sum()*dx*dx, hNew.cltz(), hNew.min());
# endif
//        printf("volNeg=%f\n",volNeg);

        copyArray(hNew,h); // Probably can do this more efficiently in loops above
        copyArray(uNew,u);
        copyArray(vNew,v);

        currentTime+=dt;

        if((currentTime-lastDisplayTime)>=displayInterval)
        {
            tmpTime=clock();
            display(ts,currentTime,Qout,integratedQout);
            lastDisplayTime=displayInterval*(int)(currentTime/displayInterval);
            totalCallbackTime+=(double)(clock()-tmpTime)/CLOCKS_PER_SEC;
        }


        if(currentTime>totalTime)
        {
            display(ts,currentTime,Qout,integratedQout);
            break;
        }

//        printf("Tick 6: %f\n",h.sum()*dx*dx);

# ifdef debug
    printf("Tick 8: %f %i %f\n",h.sum()*dx*dx, h.cltz(), h.min());
# endif

        if(covergenceTest(currentTime)) break;

    }

    printf("Total callback time: %f\n",totalCallbackTime);
    printf("Total flow comp time: %f\n",totalFlowCompTime);
    printf("Total power comp time: %f\n",totalPowerCompTime);
    return;
}


// Can probably make this more efficient by using same variables as flow calculation? Worth it?
extern "C" void courantNumber(double *hArg,double *uArg,double *vArg,int xsz,int ysz,double dt,double dx, double *crArg)
{
    array2d h, u, v, cr;
    int i, j;
    double vel,g=9.81;

    h.assign(xsz,ysz,hArg);
    u.assign(xsz+1,ysz,uArg);
    v.assign(xsz,ysz+1,vArg);
    cr.assign(xsz,ysz,crArg);

    for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) if(h(i,j)>0.01)
        cr(i,j)=(sqrt(g*h(i,j)))*dt/dx;
    else cr(i,j)=0.;

    return;
}

void hEdgeCalc()
{
    int i,j;

    for(j=0;j<ysz;j++)
    {
        for(i=1;i<xsz;i++)
            hEdgeX(i,j)=MAX(0,MAX(h(i,j)+z(i,j),h(i-1,j)+z(i-1,j))-MAX(z(i,j),z(i-1,j)));
        hEdgeX(0,j)=h(0,j);
        hEdgeX(xsz,j)=h(xsz-1,j);
    }

    for(i=0;i<xsz;i++)
    {
        for(j=1;j<ysz;j++)
            hEdgeY(i,j)=MAX(0,MAX(h(i,j)+z(i,j),h(i,j-1)+z(i,j-1))-MAX(z(i,j),z(i,j-1)));
        hEdgeY(i,0)=h(i,0);
        hEdgeY(i,ysz)=h(i,ysz-1);
    }

    return;
}

void xFlowCalc()
{
    int i,j,ts;
    double dhdx,dudx,dudy,v4,magVel,gnn;

    double divdt;

    divdt=1./dt;

    for(i=1;i<xsz;i++) for(j=0;j<ysz;j++) if(hEdgeX(i,j)>depthThresh)
    {
        dhdx=(h(i,j)+z(i,j)-h(i-1,j)-z(i-1,j))/dx;
        v4=0.25*(v(i,j)+v(i,j+1)+v(i-1,j+1)+v(i-1,j));
        magVel=sqrt(v4*v4+u(i,j)*u(i,j));

        gnn=0.5*(nGrid(i,j)+nGrid(i-1,j));
        gnn=gnn*gnn*g;

        uNew(i,j)=(u(i,j)*divdt-g*dhdx)/(divdt+gnn*magVel*fastPowM43(hEdgeX(i,j))); // Assume Manning's
    }
    else uNew(i,j)=0.;

    // Edge cells - critical depth
    for(j=0;j<ysz;j++)
    {
        uNew(0,j)=-sqrt(g*h(0,j));
        uNew(xsz,j)=sqrt(g*h(xsz-1,j));
    }

    return;
}



void yFlowCalc()
{
    int i,j,ts;
    double dhdx,dvdy,dvdx,u4,magVel,gnn;

    double divdt;

    divdt=1./dt;


// Internal points - flow in y-direction
    for(i=0;i<xsz;i++) for(j=1;j<ysz;j++) if(hEdgeY(i,j)>depthThresh)
    {
        dhdx=(h(i,j)+z(i,j)-h(i,j-1)-z(i,j-1))/dx;
        u4=0.25*(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1));
        magVel=sqrt(u4*u4+v(i,j)*v(i,j));

        gnn=0.5*(nGrid(i,j)+nGrid(i,j-1));
        gnn=gnn*gnn*g;

        vNew(i,j)=(v(i,j)*divdt-g*dhdx)/(divdt+gnn*magVel*fastPowM43(hEdgeY(i,j))); // Assume Manning's
    }
    else vNew(i,j)=0.;

    // Edge cells - critical depth
    for(i=0;i<xsz;i++)
    {
        vNew(i,0)=-sqrt(g*h(i,0));
        vNew(i,ysz)=sqrt(g*h(i,ysz-1));
    }

    return;
}


double flowOutCalc()
{
    int i,j;
    double Qout=0.0, dQ;


    for(j=0;j<ysz;j++)
    {
        Qout+=(-uNew(0,j)*hEdgeX(0,j)+uNew(xsz,j)*hEdgeX(xsz,j))*dx;
//        if(uNew(0,j)!=0)
//            printf("uNew(0,j)=%f\n",uNew(0,j));
//        if(uNew(xsz,j)!=0)
//            printf("uNew(xsz,j)=%f\n",uNew(xsz,j));

    }

//    printf("Qout loop x=%f\n",Qout);


    for(i=0;i<xsz;i++)
    {
//        dQ=(-vNew(i,0)*hEdgeY(i,0)+vNew(i,ysz)*hEdgeY(i,ysz))*dx;
        Qout+=(-vNew(i,0)*hEdgeY(i,0)+vNew(i,ysz)*hEdgeY(i,ysz))*dx;
//        if(vNew(i,0)!=0)
//            printf("vNew(i,0)=%f\n",vNew(i,0));
//        if(vNew(i,ysz)!=0)
//            printf("vNew(i,ysz)=%f\n",vNew(i,ysz));
//        if(dQ!=0)
//            printf("dQ=%f %i\n",dQ,i);

    }

//    printf("Qout loop y=%f\n",Qout);


//    printf("Qout=%f\n",Qout);

    return Qout;
}


