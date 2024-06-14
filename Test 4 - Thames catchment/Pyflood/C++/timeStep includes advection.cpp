// TO DO LIST:
// SHOULD RUN ALL LOOPS WITH j FIRST TO IMPROVE CACHING?
// Split into sensible sub functions
// Check edges for drying
// Try initialisation/run split
// Fail gracefully when NaN/very high depth produced

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
//#include <omp.h>
using namespace std;


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

class array2d
{
	public:
		double getVal(int i, int j);
		void setVal(int i, int j, double val);
		void assign(int xszArg, int yszArg, double *dataArg);
        inline double &operator()(const int &i,const int &j);
        inline double element(int i, int j);
        friend void swapArrays(array2d &a, array2d &b),copyArray(array2d a, array2d b);

//private:
        double *data;
        int xsz, ysz;

};


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


inline double pow43(const double &arg)
{
    static const double ft=1.33333333333333333333333333333333333;
    return pow(arg,-ft);

//    double it1, it2;

//    it1=sqrt(arg);
//    it2=(arg/(it1*it1)+2.*it1)*0.333333333333;
//    return 3./((arg/(it2*it2)+2.*it2)*arg);

//    return 1./arg;
//    double it1,it2,p4;

//    if(arg<0.15) it1=arg*0.333333;
//    else if(arg>5) it1=arg*3.;
//    else it1=arg;

//    p4=arg*arg;
//    p4*=p4*p4;

//    it2=0.33333333*p4/(it1*it1)+0.66666667*it1;
//    return 1./(0.33333333*p4/(it2*it2)+0.66666667*it2);
}


void xFlowCalc();
void yFlowCalc();

array2d u, v, h, uNew, vNew, hNew, z;
double depthThresh=0.01;
const double g=9.81;
bool newAdvectionMethod=true;
double dt;
double dx,beta,n,gnn,crossTermUpwind;
int xsz, ysz, frictionLaw;





inline double timeElapsed(timeval &t1, timeval &t2)
{
    return (double)(t2.tv_sec-t1.tv_sec)+1e-6*(t2.tv_usec-t1.tv_usec);
}



extern "C" double totaliser(double *arrArg,int xsz,int ysz)
{
	int i,j;
	double sum=0.;
	array2d arr;

	arr.assign(xsz,ysz,arrArg);

	for(i=0;i<xsz;i++) for (j=0;j<ysz;j++) sum+=arr(i,j);

	return sum;
}


extern "C" void timeStep(double *uArg, double *vArg, double *hArg,
    double *uNewArg, double *vNewArg, double *hNewArg,
	double *zArg,
	double *ubcVals, int *ubcLocX, int *ubcLocY, int ubcN,
	double *vbcVals, int *vbcLocX, int *vbcLocY, int vbcN,
	double *hbcVals, int *hbcLocX, int *hbcLocY, int hbcN,
	double *sbcVals, int *sbcLocX, int *sbcLocY, int sbcN,
    int xszArg, int yszArg, double totalTime, double displayInterval,
    double (*calcTimeStep)(), void (*updateBCs)(double,double), void (*display)(int,double),
    double dxArg, double nArg, double betaArg, int frictionLawArg, double crossTermUpwindArg)
{
// Set globals by arguments
    dx=dxArg;
    beta=betaArg;
    n=nArg;
    gnn=g*n*n;
    frictionLaw=frictionLawArg;
    crossTermUpwind=crossTermUpwindArg;

// Other "local" variables
    int i, j, ni, nj,ts;

    double dhdx, dhdy, dudx, dudy, dvdy, dvdx, hEdge, u4, v4, q1, q2, q3, q4, dH, dryFactor;
    double wlij, wlNeighbourMax;

    double currentTime=0.,lastDisplayTime=0.,divdt;

    clockid_t tmpTime;
    double totalCallbackTime=0., totalFlowCompTime=0.;

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

    timeval time1,time2;

    for(ts=0;;ts++)
    {
        gettimeofday(&time1,NULL);
        dt=calcTimeStep();
        divdt=1./dt;
        updateBCs(currentTime,dt);
        gettimeofday(&time2,NULL);
        totalCallbackTime+=timeElapsed(time1,time2);

    // Insert WL BCs
        for(i=0;i<hbcN;i++) h(hbcLocX[i],hbcLocY[i])=hbcVals[i]-z(hbcLocX[i],hbcLocY[i]);

    // Insert source BCs
        for(i=0;i<sbcN;i++)
        {
//            printf("%i %i %i %f %f %f\n",i,sbcLocX[i],sbcLocY[i],sbcVals[i],dt,dx);
            h(sbcLocX[i],sbcLocY[i])+=sbcVals[i]*dt/dx2;

//            printf("%f\n",h(sbcLocX[i],sbcLocY[i]));
        }
    // Calculate flows

        gettimeofday(&time1,NULL);

        xFlowCalc();
        yFlowCalc();

        gettimeofday(&time2,NULL);
        totalFlowCompTime+=timeElapsed(time1,time2);


    // Insert stuff for flow boundary conditions
        for(i=0;i<ubcN;i++) uNew(ubcLocX[i],ubcLocY[i])=ubcVals[i];
        for(i=0;i<vbcN;i++) vNew(vbcLocX[i],vbcLocY[i])=vbcVals[i];

    // Check for drying points and scale flows so that depth->0
    // Don't do this for dry cells - these may have flows specified from BCs which we don't want to interfere with
        for(i=1;i<xsz-1;i++) for(j=1;j<ysz-1;j++) if(h(i,j)>depthThresh) // TO DO: Need to check drying for edge points
        {
            q1=0.5*(h(i,j)+h(i-1,j))*uNew(i,j);
            q2=0.5*(h(i,j)+h(i+1,j))*uNew(i+1,j);
            q3=0.5*(h(i,j)+h(i,j-1))*vNew(i,j);
            q4=0.5*(h(i,j)+h(i,j+1))*vNew(i,j+1);

            dH=(q1-q2+q3-q4)*dt/dx;

            if((h(i,j)+dH)<0)
            {
                dryFactor=-h(i,j)/dH;
                uNew(i,j)*=dryFactor;
                uNew(i+1,j)*=dryFactor;
                vNew(i,j)*=dryFactor;
                vNew(i,j+1)*=dryFactor;
            }
        }

    // Update depths: interior points
    // For cells next to dry cells, need to allow for possibility we've got an internal BC adding water
        for(i=1;i<xsz-1;i++) for(j=1;j<ysz-1;j++)
        {
            if(h(i,j)>depthThresh && h(i-1,j)>depthThresh) hEdge=MAX(0,MAX(h(i,j)+z(i,j),h(i-1,j)+z(i-1,j))-MAX(z(i,j),z(i-1,j)));
            else hEdge=h(i,j);
            q1=hEdge*uNew(i,j);

            if(h(i,j)>depthThresh && h(i+1,j)>depthThresh) hEdge=MAX(0,MAX(h(i,j)+z(i,j),h(i+1,j)+z(i+1,j))-MAX(z(i,j),z(i+1,j)));
            else hEdge=h(i,j);
            q2=hEdge*uNew(i+1,j);

            if(h(i,j)>depthThresh && h(i,j-1)>depthThresh) hEdge=MAX(0,MAX(h(i,j)+z(i,j),h(i,j-1)+z(i,j-1))-MAX(z(i,j),z(i,j-1)));
            else hEdge=h(i,j);
            q3=hEdge*vNew(i,j);

            if(h(i,j)>depthThresh && h(i,j+1)>depthThresh) hEdge=MAX(0,MAX(h(i,j)+z(i,j),h(i,j+1)+z(i,j+1))-MAX(z(i,j),z(i,j+1)));
            else hEdge=h(i,j);
            q4=hEdge*vNew(i,j+1);

            hNew(i,j)=MAX(h(i,j)+(q1-q2+q3-q4)*dt/dx,0);
        }

    // Update depths: edge points
        for(i=1;i<xsz-1;i++)
        {
    // Bottom edge
            q1=0.5*(h(i,0)+h(i-1,0))*uNew(i,0);
            q2=0.5*(h(i,0)+h(i+1,0))*uNew(i+1,0);
            q3=h(i,0)*vNew(i,0);
            q4=0.5*(h(i,0)+h(i,1))*vNew(i,1);
            hNew(i,0)=h(i,0)+(q1-q2+q3-q4)*dt/dx;

    // Top edge
            q1=0.5*(h(i,ysz)+h(i-1,ysz-1))*uNew(i,ysz-1);
            q2=0.5*(h(i,ysz)+h(i+1,ysz-1))*uNew(i+1,ysz-1);
            q3=0.5*(h(i,ysz)+h(i,ysz-2))*vNew(i,ysz-1);
            q4=h(i,ysz-1)*vNew(i,ysz);
            hNew(i,ysz-1)=h(i,ysz-1)+(q1-q2+q3-q4)*dt/dx;
        }

        for(j=1;j<ysz-1;j++)
        {
    // Left edge
            q1=h(0,j)*uNew(0,j);

            if(h(0,j)>depthThresh && h(1,j)>depthThresh) hEdge=MAX(0,MAX(h(0,j)+z(0,j),h(1,j)+z(1,j))-MAX(z(0,j),z(1,j)));
            else hEdge=h(0,j);
            q2=hEdge*uNew(1,j);

//            q2=0.5*(h(0,j)+h(1,j))*uNew(1,j);

            q3=0.5*(h(0,j)+h(0,j-1))*vNew(0,j);
            q4=0.5*(h(0,j)+h(0,j+1))*vNew(0,j+1);
            hNew(0,j)=h(0,j)+(q1-q2+q3-q4)*dt/dx;

    // Right edge
            q1=0.5*(h(xsz-1,j)+h(xsz-2,j))*uNew(xsz-1,j);
            q2=h(xsz-1,j)*uNew(xsz,j);
            q3=0.5*(h(xsz-1,j)+h(xsz-1,j-1))*vNew(xsz-1,j);
            q4=0.5*(h(xsz-1,j)+h(xsz-1,j+1))*vNew(xsz-1,j+1);
            hNew(xsz-1,j)=h(xsz-1,j)+(q1-q2+q3-q4)*dt/dx;
        }

    // SW Corner points
        q1=h(0,0)*uNew(0,0);
        q2=0.5*(h(0,0)+h(1,0))*uNew(1,0);
        q3=h(0,0)*vNew(0,0);
        q4=0.5*(h(0,0)+h(0,1))*vNew(0,1);
        hNew(0,0)=h(0,0)+(q1-q2+q3-q4)*dt/dx;

    // NW Corner points
        q1=h(0,ysz-1)*uNew(0,ysz-1);
        q2=0.5*(h(0,ysz-1)+h(1,0))*uNew(1,ysz-1);
        q3=0.5*(h(0,ysz-1)+h(0,ysz-2))*vNew(0,ysz-1);
        q4=h(0,ysz-1)*vNew(0,ysz);
        hNew(0,ysz-1)=h(0,ysz-1)+(q1-q2+q3-q4)*dt/dx;

    // SE Corner points
        q1=0.5*(h(xsz-1,0)+h(xsz-2,0))*uNew(xsz-1,0);
        q2=h(xsz-1,0)*uNew(xsz,0);
        q3=h(xsz-1,0)*vNew(xsz-1,0);
        q4=0.5*(h(xsz-1,0)+h(xsz-1,1))*vNew(0,ysz);
        hNew(xsz-1,0)=h(xsz-1,0)+(q1-q2+q3-q4)*dt/dx;

    // NE Corner points
        q1=0.5*(h(xsz-1,ysz-1)+h(xsz-2,ysz-1))*uNew(xsz-1,ysz);
        q2=h(xsz-1,ysz-1)*uNew(xsz,ysz-1);
        q3=0.5*(h(xsz-1,ysz-1)+h(xsz-1,ysz-2))*vNew(xsz-1,ysz-1);
        q4=h(xsz-1,ysz-1)*vNew(xsz-1,ysz);
        hNew(xsz-1,ysz-1)=h(xsz-1,ysz-1)+(q1-q2+q3-q4)*dt/dx;

    // Wetting points
        for(i=1;i<xsz-1;i++) for(j=1;j<ysz-1;j++) if(h(i,j)<depthThresh)
        {
            wlij=z(i,j)+h(i,j);
            wlNeighbourMax=-9999;


            if(h(i+1,j)>2*depthThresh && (h(i+1,j)+z(i+1,j))>wlNeighbourMax)
            {
                wlNeighbourMax=(h(i+1,j)+z(i+1,j));
                ni=i+1;
                nj=j;
            }
            else if(h(i-1,j)>2*depthThresh && (h(i-1,j)+z(i-1,j))>wlNeighbourMax)
            {
                wlNeighbourMax=(h(i-1,j)+z(i-1,j));
                ni=i-1;
                nj=j;
            }
            else if(h(i,j+1)>2*depthThresh && (h(i,j+1)+z(i,j+1))>wlNeighbourMax)
            {
                wlNeighbourMax=(h(i,j+1)+z(i,j+1));
                ni=i;
                nj=j+1;
            }
            else if(h(i,j-1)>2*depthThresh && (h(i,j-1)+z(i,j-1))>wlNeighbourMax)
            {
                wlNeighbourMax=(h(i,j-1)+z(i,j-1));
                ni=i;
                nj=j-1;
            }

            if(wlNeighbourMax>-9999 && wlNeighbourMax>wlij) // Wet cell from highest neighbours
            {
                hNew(i,j)=depthThresh;
                hNew(ni,nj)-=(depthThresh-h(i,j));
            }
        }

    // Insert WL BCs
        for(i=0;i<hbcN;i++)
        {
            hNew(hbcLocX[i],hbcLocY[i])=hbcVals[i]-z(hbcLocX[i],hbcLocY[i]);
            if(hNew(hbcLocX[i],hbcLocY[i])<0) hNew(hbcLocX[i],hbcLocY[i])=0.;
        }


    // Clean up array
        for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) hNew(i,j)=MAX(hNew(i,j),0);

        copyArray(hNew,h); // Probably can do this more efficiently in loops above
        copyArray(uNew,u);
        copyArray(vNew,v);

        if((currentTime-lastDisplayTime)>=displayInterval || currentTime==0.)
        {
            tmpTime=clock();
            display(ts,currentTime);
            lastDisplayTime=displayInterval*(int)(currentTime/displayInterval);
            totalCallbackTime+=(double)(clock()-tmpTime)/CLOCKS_PER_SEC;
        }

        currentTime+=dt;

        if(currentTime>totalTime)
        {
            display(ts,currentTime);
            break;
        }

//        if(ts>999) break;
    }

    printf("Total callback time: %f\n",totalCallbackTime);
    printf("Total flow comp time: %f\n",totalFlowCompTime);
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
    {
        if(beta==0) cr(i,j)=(sqrt(g*h(i,j)))*dt/dx;
        else
        {
            vel=sqrt(pow(0.5*(u(i,j)+u(i+1,j)),2)+pow(0.5*(v(i,j)+v(i,j+1)),2));
            cr(i,j)=(vel+sqrt(g*h(i,j)))*dt/dx;
        }
    }
    else cr(i,j)=0.;

    return;
}


void xFlowCalc()
{
    int i,j,ts;
    double dhdx,hEdge,dudx,dudy,v4,magVel;

    double divdt;

    divdt=1./dt;

    for(i=1;i<xsz;i++) for(j=0;j<ysz;j++) if(h(i,j)>=depthThresh || h(i-1,j)>=depthThresh)
    {
        dhdx=(h(i,j)+z(i,j)-h(i-1,j)-z(i-1,j))/dx;
//        hEdge=0.5*(h(i,j)+h(i-1,j));

        hEdge=MAX(h(i,j)+z(i,j),h(i-1,j)+z(i-1,j))-MAX(z(i,j),z(i-1,j));
        if(hEdge<=0) continue;

    // Calculate velocity gradient using upwinding
        if(u(i,j)>=0) dudx=(u(i,j)-u(i-1,j))/dx;
        else dudx=(u(i+1,j)-u(i,j))/dx;

    // Advection cross terms - a bit more complicated
        v4=0.25*(v(i,j)+v(i,j+1)+v(i-1,j+1)+v(i-1,j));

        if(j==0) dudy=(u(i,j+1)-u(i,j))/dx;
        else if(j==ysz-1) dudy=(u(i,j)-u(i,j-1))/dx;
        else if(crossTermUpwind==0) dudy=0.5*(u(i,j+1)-u(i,j-1))/dx;
        else
        {
            if(v4>0) dudy=(0.5/dx)*((1.+crossTermUpwind)*(u(i,j)-u(i,j-1))+(1.-crossTermUpwind)*(u(i,j+1)-u(i,j)));
            else dudy=(0.5/dx)*((1.-crossTermUpwind)*(u(i,j)-u(i,j-1))+(1.+crossTermUpwind)*(u(i,j+1)-u(i,j)));
        }

        magVel=sqrt(v4*v4+u(i,j)*u(i,j));

        if(frictionLaw==1) uNew(i,j)=(u(i,j)*divdt-beta*v4*dudy-g*dhdx)/(divdt+gnn*magVel*(1./hEdge)+beta*dudx); //Chezy
        else uNew(i,j)=(u(i,j)*divdt-beta*v4*dudy-g*dhdx)/(divdt+gnn*magVel*pow43(hEdge)+beta*dudx); // Assume Manning's
    }
    else uNew(i,j)=0.;

    for(j=0;j<ysz;j++)
    {
        uNew(0,j)=0.;
        uNew(xsz,j)=0.;
    }

    return;
}

void yFlowCalc()
{
    int i,j,ts;
    double dhdx,hEdge,dvdy,dvdx,u4,magVel;

    double divdt;

    divdt=1./dt;


// Internal points - flow in y-direction
    for(i=0;i<xsz;i++) for(j=1;j<ysz;j++) if(h(i,j)>=depthThresh || h(i,j-1)>=depthThresh)
    {
        dhdx=(h(i,j)+z(i,j)-h(i,j-1)-z(i,j-1))/dx;
//        hEdge=0.5*(h(i,j)+h(i,j-1));

        hEdge=MAX(h(i,j)+z(i,j),h(i,j-1)+z(i,j-1))-MAX(z(i,j),z(i,j-1));
        if(hEdge<=0) continue;


// Calculate velocity gradient using upwinding
        if(v(i,j)>=0) dvdy=(v(i,j)-v(i,j-1))/dx;
        else dvdy=(v(i,j+1)-v(i,j))/dx;

// Advection cross terms - a bit more complicated
        u4=0.25*(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1));

        if(i==0) dvdx=(v(i+1,j)-v(i,j))/dx;
        else if(i==xsz-1) dvdx=(v(i,j)-v(i-1,j))/dx;
        else if(crossTermUpwind==0) dvdx=0.5*(v(i+1,j)-v(i-1,j))/dx;
        else
        {
            if(u4>0) dvdx=(0.5/dx)*((1.+crossTermUpwind)*(v(i,j)-v(i-1,j))+(1.-crossTermUpwind)*(v(i+1,j)-v(i,j)));
            else dvdx=(0.5/dx)*((1.-crossTermUpwind)*(v(i,j)-v(i-1,j))+(1.+crossTermUpwind)*(v(i+1,j)-v(i,j)));
        }

        magVel=sqrt(u4*u4+v(i,j)*v(i,j));

        if(frictionLaw==1) vNew(i,j)=(v(i,j)*divdt-beta*u4*dvdx-g*dhdx)/(divdt+gnn*magVel*(1./hEdge)+beta*dvdy); //Chezy
        else vNew(i,j)=(v(i,j)*divdt-beta*u4*dvdx-g*dhdx)/(divdt+gnn*magVel*pow43(hEdge)+beta*dvdy); // Assume Manning's
    }
    else vNew(i,j)=0.;


    for(i=0;i<xsz;i++)
    {
        vNew(i,0)=0.;
        vNew(i,ysz)=0.;
    }

    return;
}





