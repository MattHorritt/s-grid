// Compile with:
//  g++ -O3 -c -fPIC sgridHydraulics.cpp -o sgridHydraulics.o
//  g++ -shared -Wl,-soname,sgridHydraulics.o -o sgridHydraulics.so sgridHydraulics.o

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <stdlib.h>
#include <ctime>
#include <algorithm>

using namespace std;

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define weirCoeff 1.705

class array2d
{
	public:
		float getVal(int i, int j);
		void setVal(int i, int j, float val);
		void assign(int xszArg, int yszArg, float *dataArg);
        inline float &operator()(const int &i,const int &j);
        inline float element(int i, int j);
        friend void swapArrays(array2d &a, array2d &b),copyArray(array2d a, array2d b);

//private:
        float *data;
        int xsz, ysz;

};

// Utility to check licence and print status - useful diagnostic, but don't
// want to run every time the licence is checked (e.g. every time step)
extern "C" void checkLicence()
{
    time_t timeNow=time(NULL);

    struct tm tm = {0};


    tm.tm_year=2024-1900;
    tm.tm_mon=12;
    tm.tm_mday=31;

    time_t timeExpire=mktime(&tm);

    if(difftime(timeExpire,timeNow)<0)
	{
		printf("\nLicence expired %02i/%02i/%i - stopping\n",tm.tm_mday,tm.tm_mon,tm.tm_year+1900);
           fflush(stdout);

		exit(1);
	}
    else
    {
        printf("\nLicence valid until %02i/%02i/%i\n",tm.tm_mday,tm.tm_mon,tm.tm_year+1900);
        return;
    }
}

inline void timeBomb()
{
    time_t timeNow=time(NULL);

    struct tm tm = {0};


    tm.tm_year=2024-1900;
    tm.tm_mon=12;
    tm.tm_mday=31;

    time_t timeExpire=mktime(&tm);

    if(difftime(timeExpire,timeNow)<0)
	{
		printf("\nLicence expired %i/%02i/%02i - stopping",tm.tm_year+1900,tm.tm_mon,tm.tm_mday);
           fflush(stdout);

		exit(1);
	}
}


inline void copyArray(array2d src, array2d dst)
{
    memcpy(dst.data,src.data,sizeof(float)*src.xsz*src.ysz);
    return;
}


inline void swapArrays(array2d &a, array2d &b)
{
    // Assume they're of equal size - no error checking here
    float *tmp=a.data;
    a.data=b.data;
    b.data=tmp;
}


inline void array2d::assign(int xszArg, int yszArg, float *dataArg)
{
	data=dataArg;
	xsz=xszArg;
	ysz=yszArg;
}

inline float &array2d::operator()(const int &i,const int &j)
{
    return data[i*ysz+j];
}

//array2d u, v, h, uNew, vNew, hNew, z;

extern "C" void maxVolGrid(float *v, float *mv,
	float *fx, float *mfx, float *fy, float *mfy,
	int xsz, int ysz)
{
	array2d vol,maxVol,fx2,fy2,mfx2,mfy2;
	int i,j;

	vol.assign(xsz,ysz,v);
	maxVol.assign(xsz,ysz,mv);
	fx2.assign(xsz+1,ysz,fx);
	mfx2.assign(xsz+1,ysz,mfx);
	fy2.assign(xsz,ysz+1,fy);
	mfy2.assign(xsz,ysz+1,mfy);

	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
		maxVol(i,j)=MAX(vol(i,j),maxVol(i,j));

	for(i=0;i<=xsz;i++) for(j=0;j<ysz;j++)
		if(fabs(fx2(i,j))>fabs(mfx2(i,j)))
			mfx2(i,j)=fx2(i,j);

	for(i=0;i<xsz;i++) for(j=0;j<=ysz;j++)
		if(fabs(fy2(i,j))>fabs(mfy2(i,j)))
			mfy2(i,j)=fy2(i,j);

	return;
}

//=================================================================================
extern "C" void calcStorageParameters(float x0, float y0, float x1, float y1,
	float *dtmArg, int xsz, int ysz, float xll, float yll, float dx,
	float *retArray,
    bool channel, float *cagArg,int cagXsz,int cagYsz,
       float cagXll, float cagYll, float cagDx,
    float chanMult, float chanExp, float chanAR, float maxChanD, float gridSize)
{
	int xi0, xi1, yi0, yi1;
	float zMin=1e20, zMax=-1e20, zBank;
	int i,j,k,wl;
	array2d dtm,catchmentAreaGrid;
    float *volIP, *wlList, *cellList;
    int numCells, numCellsChannel;
    float chanWidth, chanDepth, chanArea, maxArea;
    int cag_xi0,cag_xi1,cag_yi0,cag_yi1;

    float noDataValue=-9999.;
    int countNnd=0, countAll=0;

    timeBomb();

	dtm.assign(xsz,ysz,dtmArg);

	xi0=int((x0-xll)/dx);
	xi1=int((x1-xll)/dx);
	yi0=int((y0-yll)/dx);
	yi1=int((y1-yll)/dx);

	for(i=xi0;i<=xi1;i++) for(j=yi0;j<=yi1;j++)
    {
        countAll++;
        if(dtm(i,j)!=noDataValue)
        	{
        		zMin=MIN(zMin,dtm(i,j));
        		zMax=MAX(zMax,dtm(i,j));

                countNnd++;
        	}
    }

    if((countNnd/float(countAll))<0.5)
    {
       	*(retArray)=noDataValue;
       	*(retArray+1)=noDataValue;
       	*(retArray+2)=noDataValue;
       	*(retArray+3)=noDataValue;
       	*(retArray+4)=noDataValue;
       	*(retArray+5)=noDataValue;
       	*(retArray+6)=noDataValue;

//       	*(retArray+7)=noDataValue;

            return;
    }



    if(channel)
    {
// Channel parameters
	    catchmentAreaGrid.assign(cagXsz,cagYsz,cagArg);

       	cag_xi0=int((x0-cagXll)/cagDx);
       	cag_xi1=int((x1-cagXll)/cagDx);
       	cag_yi0=int((y0-cagYll)/cagDx);
       	cag_yi1=int((y1-cagYll)/cagDx);

        maxArea=-1e20;

        if((cag_xi0>=0)&&(cag_xi0<cagXsz)&&
           (cag_yi0>=0)&&(cag_yi0<cagYsz)&&
           (cag_xi1>=0)&&(cag_xi1<cagXsz)&&
           (cag_yi1>=0)&&(cag_yi1<cagYsz))
        {
	       	for(i=cag_xi0;i<cag_xi1;i++) for(j=cag_yi0;j<cag_yi1;j++)
	        maxArea=MAX(catchmentAreaGrid(i,j),maxArea);
        }
        else maxArea=0;

		if(maxArea<3.0) maxArea=0.0;


        chanWidth=chanMult*pow(maxArea,chanExp);
        chanDepth=MIN(chanWidth/chanAR,maxChanD);
        chanArea=chanWidth*gridSize;

//		printf("maxArea=%f\tchanWidth=%f\tchanDepth=%f\tchanArea=%f\n",maxArea,chanWidth,chanDepth,chanArea);


        numCellsChannel=int(chanArea/(dx*dx));

        zBank=zMin;
        zMin-=chanDepth;

        wlList=new float[5];

        if(zMax>(zBank+5))
        {
    		wlList[0]=zMin;
    		wlList[1]=zBank;
    		wlList[2]=zBank+1;
    		wlList[3]=zBank+5;
    		wlList[4]=zMax;
        }
        else
        {
    		wlList[0]=zMin;
    		wlList[1]=zBank;
    		wlList[2]=zBank+1;
    		wlList[3]=zBank+5;
    		wlList[4]=zBank+5;

    		zMax=zBank+5;
        }


        volIP=new float[5];
        volIP[0]=0.;

		numCells=(abs(xi1-xi0)+1)*(abs(yi1-yi0)+1);

        cellList=new float[numCells];

        k=0;
        for(i=xi0;i<=xi1;i++) for(j=yi0;j<=yi1;j++) cellList[k++]=dtm(i,j);

        std::sort(cellList,cellList+numCells);

		volIP[1]=chanArea*chanDepth/(dx*dx*numCells);

      	for(wl=2;wl<5;wl++)
       	{
       		volIP[wl]=chanArea*(wlList[wl]-zMin)/(dx*dx);
            for(i=0;i<numCells;i++)
				if(i==numCellsChannel)
				{
					volIP[wl]+=((numCellsChannel+1)-chanArea/(dx*dx))*MAX(0.,wlList[wl]-cellList[i]);
				//	printf("%f %f %f %f\n",(numCellsChannel+1)*dx*dx-chanArea,wl,cellList[i],volIP[wl]);
				}
				else if(i>numCellsChannel)
					volIP[wl]+=MAX(0.,wlList[wl]-cellList[i]);
       		volIP[wl]/=numCells;
       	}

       	*(retArray)=zMin;
       	*(retArray+1)=zBank;
       	*(retArray+2)=zMax;
       	*(retArray+3)=volIP[1];
       	*(retArray+4)=volIP[2];
       	*(retArray+5)=volIP[3];
       	*(retArray+6)=volIP[4];

	}
	else
    {
        	wlList=new float[4];

        	if(zMax>(zMin+5) && zMax>(zMin+1))
        	{
        		wlList[0]=zMin;
        		wlList[1]=zMin+1;
        		wlList[2]=zMin+5;
        		wlList[3]=zMax;

        	}
        	else
        	{
        		wlList[0]=zMin;
        		wlList[1]=zMin+1;
        		wlList[2]=zMin+5;
        		wlList[3]=zMin+5;

        		zMax=zMin+5;
        	}


        	volIP=new float[4];
            volIP[0]=0.;
        	for(wl=1;wl<4;wl++)
        	{
        		volIP[wl]=0.;
        		for(i=xi0;i<=xi1;i++) for(j=yi0;j<=yi1;j++)
        			if(dtm(i,j)<wlList[wl] && !isnan(dtm(i,j))) volIP[wl]+=(wlList[wl]-dtm(i,j));
        		volIP[wl]/=((abs(xi1-xi0)+1)*(abs(yi1-yi0)+1));
        	}

        	*(retArray)=zMin;
        	*(retArray+1)=zMax;
        	*(retArray+2)=volIP[1];
        	*(retArray+3)=volIP[2];
        	*(retArray+4)=volIP[3];
    }




	return;
}


void logLinearRegression(int n, float *x, float *y, float &s,float &i)
{
	float Sx=0, Sy=0, Sxy=0, Sxx=0;
	int ii;

	for(ii=0;ii<n;ii++)
	{
		x[ii]=log(x[ii]);
		y[ii]=log(y[ii]);

		Sx+=x[ii];
		Sxx+=x[ii]*x[ii];
		Sxy+=x[ii]*y[ii];
		Sy+=y[ii];
	}

	s=(Sxy-Sx*Sy/n)/(Sxx-Sx*Sx/n);
	i=(Sy-s*Sx)/n;

	return;
}

//==============================================================================
extern "C" void conveyanceParameters(int xi0, int yi0, int xi1, int yi1,
                    float *dtmArg,float dtmDx, int xsz, int ysz ,
                    float manningsN,
                    float *returnSext,
                    bool channel, float nChan, float *cagArg,int cagXsz,int cagYsz,
                    int cag_xi0, int cag_yi0, int cag_xi1, int cag_yi1,
                    float chanMult, float chanExp, float chanAR, float chanMaxD,
                    float *manningsGridArg)
{
	array2d dtm,catchmentAreaGrid,manningsGrid;
	int profileLength,cagProfileLength;
	float *profile,*nProfile,crossSectionLength,maxDepth=10.,maxZ=-1e20,dl,wl;
    float maxArea;
	int nConveyanceLevels=10;
	int i,j,pptr,n;
	float minZ, s1, s2, icpt1, icpt2, bankZ;

    float chanWidth, chanDepth, *nList;
    int numCellsChannel;

    timeBomb();

	dtm.assign(xsz,ysz,dtmArg);
    manningsGrid.assign(xsz,ysz,manningsGridArg);

	profileLength=MAX(abs(xi1-xi0),abs(yi1-yi0))+1;
	profile=new float[profileLength];
     nProfile=new float[profileLength];

    crossSectionLength=sqrt((xi1-xi0)*(xi1-xi0)+(yi1-yi0)*(yi1-yi0));

  	int dx = abs(xi1-xi0), sx = xi0<xi1 ? 1 : -1;
  	int dy = abs(yi1-yi0), sy = yi0<yi1 ? 1 : -1;
  	int err = (dx>dy ? dx : -dy)/2, e2;

 	pptr=0;


    // Create profile using only non-NaN values
    profileLength = 0;
	for(;;)
	{
	    if(!isnan(dtm(xi0,yi0)))
	    {

            profile[pptr]=dtm(xi0,yi0);
            nProfile[pptr]=manningsGrid(xi0,yi0);

               pptr++;
               profileLength++;
        }

		if (xi0==xi1 && yi0==yi1) break;
		e2 = err;
		if (e2 >-dx) { err -= dy; xi0 += sx; }
		if (e2 < dy) { err += dx; yi0 += sy; }
	}

//    printf("\n");
//    for(i=0;i<profileLength;i++) printf("%f %i ", profile[i], isnan(profile[i]));
//    printf("\n");

    dl=dtmDx*crossSectionLength/profileLength;

	minZ=1e20;
    maxZ=-1e20;
	for(i=0;i<profileLength;i++)
	{
		minZ=MIN(profile[i],minZ);
		maxZ=MAX(profile[i],maxZ);
	}


// Burn in channel cells
    if(channel)
    {
		std::sort(profile,profile+profileLength); // First sort the existing profile

        catchmentAreaGrid.assign(cagXsz,cagYsz,cagArg);

        cagProfileLength=MAX(abs(cag_xi1-cag_xi0),abs(cag_yi1-cag_yi0)); // Don't use last point

           // Check start and finish are within grid
        if((cag_xi0>=0)&&(cag_xi0<cagXsz)&&
            (cag_yi0>=0)&&(cag_yi0<cagYsz)&&
            (cag_xi1>=0)&&(cag_xi1<cagXsz)&&
            (cag_yi1>=0)&&(cag_yi1<cagYsz))
        {

        	dx = abs(cag_xi1-cag_xi0), sx = cag_xi0<cag_xi1 ? 1 : -1;
          	dy = abs(cag_yi1-cag_yi0), sy = cag_yi0<cag_yi1 ? 1 : -1;
          	err = (dx>dy ? dx : -dy)/2, e2;

            maxArea=-1;

            pptr=0;

            for(;;)
            {
            	maxArea=MAX(catchmentAreaGrid(cag_xi0,cag_yi0),maxArea);
                pptr++;
                if(pptr==cagProfileLength) break;
            	if (cag_xi0==cag_xi1 && cag_yi0==cag_yi1) break;
            	e2 = err;
            	if (e2 >-dx) { err -= dy; cag_xi0 += sx; }
            	if (e2 < dy) { err += dx; cag_yi0 += sy; }
            }
        }
        else maxArea=0.;

		if(maxArea<3.0) maxArea=0.0;

        chanWidth=chanMult*pow(maxArea,chanExp);
        chanDepth=MIN(chanWidth/chanAR,chanMaxD);

        numCellsChannel=int(chanWidth/dl);
        bankZ=minZ;
        minZ-=chanDepth;
    }
    else
    {
    	bankZ=minZ;
    }

	float *conveyance=new float[nConveyanceLevels];
	float *conveyanceUsingMinZ=new float[nConveyanceLevels];

    float mnt;

	// Slope and intercept for depths below bankfull
	if(channel and maxArea>0.)
	{
		s1=1.0;
		icpt1=-log((dl*profileLength/chanWidth)*(nChan/manningsN));
	}
	else
	{
		s1=-1; //Flag for no channel in this edge
	}


    float nAvg=0.0; // Average Manning's n as a reference and for calculating A later
    for(i=0;i<profileLength;i++) nAvg+=nProfile[i];
    nAvg/=profileLength;


	// Slope and intercept for depths above bankfull
	for(i=0;i<nConveyanceLevels;i++)
	{
		wl=bankZ+maxDepth*pow(10,(-2.+i*2./((nConveyanceLevels-1))));

        conveyanceUsingMinZ[i]=(dl*profileLength*pow(wl-minZ,1.666667)/nAvg);

		conveyance[i]=0.;
		if(channel)
		{
			conveyance[i]+=chanWidth*pow(wl-minZ,1.6666667)/nChan;

			for(j=numCellsChannel;j<profileLength;j++)
				if(j==numCellsChannel)
					conveyance[i]+=((numCellsChannel+1)*dl-chanWidth)*pow(MAX(0.,wl-profile[j]),1.6666667)/manningsN;
				else if(j>numCellsChannel)
					conveyance[i]+=dl*pow(MAX(0.,wl-profile[j]),1.6666667)/manningsN;
		}
		else for(j=0;j<profileLength;j++)
			conveyance[i]+=dl*pow(MAX(0.,wl-profile[j]),1.6666667)/nProfile[j];

//		if(chanDepth>=5) printf("i=%i \t wl=%f \t conveyance=%f\n",i,wl,conveyance[i]);
	}

	logLinearRegression(nConveyanceLevels,conveyanceUsingMinZ,conveyance,s2,icpt2);

//	if(chanDepth>=5) printf("Slope=%f Intercept=%f\n",slope,intercept);

	*(returnSext)=minZ;
	*(returnSext+1)=bankZ;
	*(returnSext+2)=s1;
	*(returnSext+3)=icpt1;
	*(returnSext+4)=s2;
	*(returnSext+5)=icpt2;
	*(returnSext+6)=nAvg;

	return;
}
//=================================================================================

extern "C" float sum(float *arrArg, int xsz, int ysz)
{
	array2d arr;
	int i,j;
	float sum=0.;

	printf("In sum()...\n");

	arr.assign(xsz,ysz,arrArg);

	printf("Array allocated...\n");

	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
		sum+=arr(i,j);

	printf("Done\n");

	return(sum);
}





extern "C" void timeStep(float *volArg, float *qxArg, float *qyArg, float dt,
	int xsz, int ysz,
	int *fpXi, int *fpYi, float *fpQ, int nfp)
{
	array2d vol, qx, qy;
	float q1, q2, q3, q4;
	int i,j;

    timeBomb();

// Set up array2d types
	vol.assign(xsz,ysz,volArg);
	qx.assign(xsz+1,ysz,qxArg);
	qy.assign(xsz,ysz+1,qyArg);

//    cout<<"nfp="<<nfp<<endl;

	if(nfp>0) for(i=0;i<nfp;i++)
    {
        vol(fpXi[i],fpYi[i])+=dt*fpQ[i];
 //       cout<<fpXi[i]<<" "<<fpYi[i]<<" "<<fpQ[i]*dt<<endl;
    }
	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		q1=qx(i,j);
		q2=-qx(i+1,j);
		q3=qy(i,j);
		q4=-qy(i,j+1);

		vol(i,j)+=dt*(q1+q2+q3+q4);
	}

	return;
}

extern "C" int dryCheck(float *volArg, float *qxArg, float *qyArg,
	float dt, int xsz, int ysz,
	int *fpXi, int *fpYi, float *fpQ, int nfp)
{
	int maxIt=20,it,i,j;
	int nChanged;
	float Q1t, Q2t, Q3t, Q4t,alpha;
	array2d vol, Qx, Qy;

// Set up array2d types
	vol.assign(xsz,ysz,volArg);
	Qx.assign(xsz+1,ysz,qxArg);
	Qy.assign(xsz,ysz+1,qyArg);

	for(it=0;it<maxIt;it++)
	{
		nChanged=0;

		if(nfp>0) for(i=0;i<nfp;i++)
		{
			vol(fpXi[i],fpYi[i])+=dt*fpQ[i];
		}

		for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
		{
		    Q1t=-Qx(i,j)*dt; // Flows out
            Q2t=Qx(i+1,j)*dt;
            Q3t=-Qy(i,j)*dt;
            Q4t=Qy(i,j+1)*dt;

		    if((Q1t+Q2t+Q3t+Q4t)>(vol(i,j)+1.))
		    {
                alpha=vol(i,j)/(Q1t+Q2t+Q3t+Q4t);

		        if(Q1t>0)
		        {
		            Qx(i,j)*=alpha;
		            nChanged+=1;
		        }
		        if(Q2t>0)
		        {
		            Qx(i+1,j)*=alpha;
		            nChanged+=1;
		        }
		        if(Q3t>0)
		        {
		            Qy(i,j)*=alpha;
		            nChanged+=1;
		        }
		        if(Q4t>0)
		        {
		            Qy(i,j+1)*=alpha;
		            nChanged+=1;
		        }
		    }
        }

        if(nfp>0) for(i=0;i<nfp;i++) vol(fpXi[i],fpYi[i])-=dt*fpQ[i];

 //       printf("In dryCheck, %i iterations, %i values changed.\n",it,nChanged);

        if(nChanged==0) break;
    }

    return(it+1);
}


extern "C" int dryCheckDiagnostic(float *volArg, float *qxArg, float *qyArg,
	float dt, int xsz, int ysz,
	int *fpXi, int *fpYi, float *fpQ, int nfp,float *nLimitArg)
{
	int maxIt=10,it,i,j;
	int nChanged;
	float Q1t, Q2t, Q3t, Q4t,alpha;
	array2d vol, Qx, Qy,nLimit;

// Set up array2d types
	vol.assign(xsz,ysz,volArg);
	Qx.assign(xsz+1,ysz,qxArg);
	Qy.assign(xsz,ysz+1,qyArg);
     nLimit.assign(xsz,ysz,nLimitArg);

	for(it=0;it<maxIt;it++)
	{
		nChanged=0;

		if(nfp>0) for(i=0;i<nfp;i++)
		{
			vol(fpXi[i],fpYi[i])+=dt*fpQ[i];
		}

		for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
		{
		    Q1t=-Qx(i,j)*dt; // Flows out
            Q2t=Qx(i+1,j)*dt;
            Q3t=-Qy(i,j)*dt;
            Q4t=Qy(i,j+1)*dt;

		    if((Q1t+Q2t+Q3t+Q4t)>(vol(i,j)+1.))
		    {
                alpha=vol(i,j)/(Q1t+Q2t+Q3t+Q4t);

                nLimit(i,j)+=1;

		        if(Q1t>0)
		        {
		            Qx(i,j)*=alpha;
		            nChanged+=1;
		        }
		        if(Q2t>0)
		        {
		            Qx(i+1,j)*=alpha;
		            nChanged+=1;
		        }
		        if(Q3t>0)
		        {
		            Qy(i,j)*=alpha;
		            nChanged+=1;
		        }
		        if(Q4t>0)
		        {
		            Qy(i,j+1)*=alpha;
		            nChanged+=1;
		        }
		    }
        }

        if(nfp>0) for(i=0;i<nfp;i++) vol(fpXi[i],fpYi[i])-=dt*fpQ[i];

 //       printf("In dryCheck, %i iterations, %i values changed.\n",it,nChanged);

        if(nChanged==0) break;
    }

    return(it+1);
}

extern "C" float calcFlowInertial(float Q, float wl1, float wl2, float zEdge, float zBank,
	float s1, float i1, float s2, float i2,
	float dxGrid, float n, float dt)
{
	float wl,kMin,k,S,flow;
	float g=9.81,h,A,flowNI;
	float slope,intercept;
    float steepFlowLim=1e20;


    if(MAX(wl1,wl2) <=zEdge || zEdge==-9999) return 0.;

    wl=0.5*(wl1+wl2);
    if(wl<zEdge)
    {
        wl=MAX(wl1,wl2);
//        steepFlowLim=(wl-zEdge)*dxGrid*dxGrid/dt;
    }

//	wl=MAX(wl1,wl2);

	if(wl<=zBank && s1>=0.) // s1 indicates whether channel is present
	{
		slope=s1;
		intercept=i1;
	}
	else
	{
		slope=s2;
		intercept=i2;
	}

    kMin=(dxGrid*(pow((wl-zEdge),1.666667))/n);
    k=exp(slope*log(kMin)+intercept);

	h=wl-zEdge;
	A=MAX(1.,k*n/pow(h,0.6667));

    S=(wl2-wl1)/dxGrid;

    flow=(Q-g*S*A*dt)/(1.+g*A*dt*fabs(Q)/(k*k));

    if(S>0) flowNI=-k*sqrt(S);
    else flowNI=k*sqrt(-S);

	if(isnan(flow)) flow=flowNI;

//    flow=MAX(flow,steepFlowLim);

    if(isnan(flow))
    {
        cout<<"Q"<<Q<<endl;
        cout<<"wl1"<<wl1<<endl;
        cout<<"wl2"<<wl2<<endl;
        cout<<"zEdge"<<zEdge<<endl;
        cout<<"zBank"<<zBank<<endl;
        cout<<"s1"<<s1<<endl;
        cout<<"i1"<<i1<<endl;
        cout<<"s2"<<s2<<endl;
        cout<<"i2"<<i2<<endl;
        cout<<"dxGrid"<<dxGrid<<endl;
        cout<<"n"<<n<<endl;
        cout<<"dt"<<dt<<endl;
    }

    return(flow);
}

extern "C" float calcFlow(float wl1, float wl2, float zEdge, float slope, float intercept, float dxGrid, float n)
{
	float wl,kMin,k,S,flow;

	wl=MAX(wl1,wl2);

    if(wl<=zEdge || zEdge==-9999) return 0.;

    kMin=(dxGrid*(pow((wl-zEdge),1.666667))/n);
    k=exp(slope*log(kMin)+intercept);

    S=(wl2-wl1)/dxGrid;

    if(S>0) flow=-k*sqrt(S);
    else flow=k*sqrt(-S);

    return(flow);
}

extern "C" float calcFlowEdges(float *wlArg,
	float *cvpX,
     float *cvpY,
      float *sp,
	float cellSize, int xsz, int ysz,
	float *qxArg, float *qyArg, float *dryMaskArg)
{
    array2d wl,qx,qy,dryMask;

	wl.assign(xsz,ysz,wlArg);
	qx.assign(xsz+1,ysz,qxArg);
	qy.assign(xsz,ysz+1,qyArg);
     dryMask.assign(xsz,ysz,dryMaskArg);

        float flow,wl1,wl2,z1,z2,ze,flowOut=0;
    int i,j;

    flowOut=0.;

	for(i=1;i<xsz;i++) for(j=0;j<ysz;j++)
    {
           wl1=wl(i-1,j);
		wl2=wl(i,j);

            z1=sp[((i-1)*ysz+j)*5+0];
            z2=sp[(i*ysz+j)*5+0];

//		ze=cvpX[(i*ysz+j)*6+0]+1.;



        // Flow to the right
            if(z1>-9999 && z2==-9999 && dryMask(i-1,j)!=0 && wl1>(z1+1.0))
            {
                ze=z1+1.0;

                flow=cellSize*weirCoeff*pow(wl1-ze,1.5);


                qx(i,j)=flow;
                flowOut+=flow;

//                printf("%i,%i\t%f\t%f\t%f\t%f\n",i,j,wl1,ze,flow,flowOut);
            }

        // Flow to the left
            if(z2>-9999 && z1==-9999 && dryMask(i,j)!=0 && wl2>(z2+1.0))
            {
                ze=z2+1.0;

                flow=-cellSize*weirCoeff*pow(wl2-ze,1.5);
                qx(i,j)=flow;
                flowOut-=flow;

//                    printf("%i,%i\t%f\t%f\t%f\t%f\n",i,j,wl2,ze,flow,flowOut);


            }
    }


	for(i=0;i<xsz;i++) for(j=1;j<ysz;j++)
    {
           wl1=wl(i,j-1);
		wl2=wl(i,j);

            z1=sp[(i*ysz+j-1)*5+0];
            z2=sp[(i*ysz+j)*5+0];

//    		ze=cvpY[(i*(ysz+1)+j)*6+0]+1.;


        // Flow up
            if(z1>-9999 && z2==-9999 && dryMask(i,j-1)!=0 && wl1>(z1+1.0))
            {
                ze=z1+1.0;

                flow=cellSize*weirCoeff*pow(wl1-ze,1.5);
                qy(i,j)=flow;
                flowOut+=flow;


//                    printf("%i,%i\t%f\t%f\t%f\t%f\n",i,j,wl1,ze,flow,flowOut);


            }

        // Flow down
            if(z2>-9999 && z1==-9999 && dryMask(i,j)!=0 && wl2>(z2+1.0))
            {
                ze=z2+1.0;

                flow=-cellSize*weirCoeff*pow(wl2-ze,1.5);
                qy(i,j)=flow;
                flowOut-=flow;

//                    printf("%i,%i\t%f\t%f\t%f\t%f\n",i,j,wl2,ze,flow,flowOut);

            }


    }


    return flowOut;

}


extern "C" void calcFlowGrid(float *wlArg,
	float *cvpX,
	float *cvpY,
	float dx, float dt, int xsz, int ysz,
	float *qxArg, float *qyArg, float *dryMaskArg)
{
	array2d wl,qx,qy,dryMask;
	int i,j;
	float wl1, wl2, ze, sl1, intercept1, sl2, intercept2, zb;
     float nAvg;

    timeBomb();

	wl.assign(xsz,ysz,wlArg);
	qx.assign(xsz+1,ysz,qxArg);
	qy.assign(xsz,ysz+1,qyArg);
     dryMask.assign(xsz,ysz,dryMaskArg);

	for(i=1;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		wl1=wl(i-1,j);
		wl2=wl(i,j);
		ze=cvpX[(i*ysz+j)*7+0];

		zb=cvpX[(i*ysz+j)*7+1];

		sl1=cvpX[(i*ysz+j)*7+2];
		intercept1=cvpX[(i*ysz+j)*7+3];
		sl2=cvpX[(i*ysz+j)*7+4];
		intercept2=cvpX[(i*ysz+j)*7+5];

          nAvg=cvpX[(i*ysz+j)*7+6];

		if(ze==-9999.) qx(i,j)=0.;
		else qx(i,j)=calcFlowInertial(qx(i,j),wl1,wl2,ze,zb,sl1,intercept1,sl2,intercept2,dx,nAvg,dt);

            if (qx(i,j)>0 and dryMask(i-1,j)==0) // Flow from dry cell
                qx(i,j)=0;

            if (qx(i,j)<0 and dryMask(i,j)==0) // Flow from dry cell
                qx(i,j)=0;

//		qx(i,j)=calcFlow(wl1,wl2,ze,sl,intercept,dx,n);
	}

	for(i=0;i<xsz;i++) for(j=1;j<ysz;j++)
	{
		wl1=wl(i,j-1);
		wl2=wl(i,j);
		ze=cvpY[(i*(ysz+1)+j)*7+0];

		zb=cvpY[(i*(ysz+1)+j)*7+1];

		sl1=cvpY[(i*(ysz+1)+j)*7+2];
		intercept1=cvpY[(i*(ysz+1)+j)*7+3];
		sl2=cvpY[(i*(ysz+1)+j)*7+4];
		intercept2=cvpY[(i*(ysz+1)+j)*7+5];

          nAvg=cvpY[(i*(ysz+1)+j)*7+6];

		if(ze==-9999.) qx(i,j)=0.;
		else qy(i,j)=calcFlowInertial(qy(i,j),wl1,wl2,ze,zb,sl1,intercept1,sl2,intercept2,dx,nAvg,dt);
//		qy(i,j)=calcFlow(wl1,wl2,ze,sl,intercept,dx,n);

            if (qy(i,j)>0 and dryMask(i,j-1)==0) // Flow from dry cell
                qy(i,j)=0;

            if (qy(i,j)<0 and dryMask(i,j)==0) // Flow from dry cell
                qy(i,j)=0;

	}

	return;
}

float interp4(float x,float x1,float x2,float x3,float x4,
	float y1,float y2,float y3,float y4)
{
	if(x<=x1) return(y1);
	if(x>=x4) return(y4);

	if(x>x1 and x<=x2) return(y1+(y2-y1)*(x-x1)/(x2-x1));
	if(x>x2 and x<=x3) return(y2+(y3-y2)*(x-x2)/(x3-x2));
	if(x>x3 and x<=x4) return(y3+(y4-y3)*(x-x3)/(x4-x3));

	return -9999;
}

float interp5(float x,float x1,float x2,float x3,float x4, float x5,
	float y1,float y2,float y3,float y4, float y5)
{
	if(x<=x1) return(y1);
	if(x>=x5) return(y5);

	if(x>x1 and x<=x2) return(y1+(y2-y1)*(x-x1)/(x2-x1));
	if(x>x2 and x<=x3) return(y2+(y3-y2)*(x-x2)/(x3-x2));
	if(x>x3 and x<=x4) return(y3+(y4-y3)*(x-x3)/(x4-x3));
	if(x>x4 and x<=x5) return(y4+(y5-y4)*(x-x4)/(x5-x4));

	return -9999;
}

float wlFromVol4(float v,float zMin, float zMax, float vip1, float vip2, float vip3, float dx)
{
    if(v/(dx*dx)>=vip3) return((v/(dx*dx)-vip3)+zMax);
    else return(interp4(v/(dx*dx),0,vip1,vip2,vip3,zMin,zMin+1,zMin+5,zMax));
}

float wlFromVol5(float v,float zMin, float zBank, float zMax,
    float vip1, float vip2, float vip3, float vip4, float dx)
{
    if(v/(dx*dx)>=vip4) return((v/(dx*dx)-vip4)+zMax);
    else return(interp5(v/(dx*dx),0,vip1,vip2,vip3,vip4,zMin,zBank,zBank+1,zBank+5,zMax));
}


extern "C" void wlFromVolGrid(float *v, float *wl, float *sp, int xsz, int ysz,
    bool channel, float dx)
{
	int i,j;

//    printf("channel = %s\n",channel  ? "true" : "false"  );

    if(channel)
    {

        	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
        		wl[i*ysz+j]=wlFromVol5(v[i*ysz+j],
        			sp[(i*ysz+j)*7+0],
        			sp[(i*ysz+j)*7+1],
        			sp[(i*ysz+j)*7+2],
        			sp[(i*ysz+j)*7+3],
        			sp[(i*ysz+j)*7+4],
        			sp[(i*ysz+j)*7+5],
        			sp[(i*ysz+j)*7+6],
        			dx);
    }
    else
    {
        	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
        		wl[i*ysz+j]=wlFromVol4(v[i*ysz+j],
        			sp[(i*ysz+j)*5+0],
        			sp[(i*ysz+j)*5+1],
        			sp[(i*ysz+j)*5+2],
        			sp[(i*ysz+j)*5+3],
        			sp[(i*ysz+j)*5+4],
        			dx);
    }

    return;
}



inline int nint(float f)
{
	return(int(f+0.5));
}

float bilinear(float x, float y, float z1, float z2, float z3, float z4)
{
	float l1, l2;

	l1=(y+1-x)/2;
	l2=(-y-x+1)/2;

	return (l1*z2+(1-l1)*z1)*(1-l2)+(l1*z3+(1-l1)*z4)*l2;
}

extern "C" void fillWlGrid(float *wlArg,float *wlOutArg, int xsz, int ysz, float dx,
    int nExpansionIts,int nSmoothingIts,float expansionSlope)
{
	int it,i,j;
	int nnzNeighbours;
	float sumNeighbours;

	array2d wlIn, wlOut, tmp;

	wlIn.assign(xsz,ysz,wlArg);
	wlOut.assign(xsz,ysz,wlOutArg);

	copyArray(wlIn,wlOut);

    // Create a bit set and set to 1 where wlIn!=-9999
    // This will be used to remember where our original data is, and where has been expanded
    // we'll then only smooth the expanded bits, maintaining the original cells' levels
    char *maskArr=new char[xsz*ysz];
    for(i=0;i<(xsz*ysz);i++) maskArr[i]=0;

    for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) if(wlIn(i,j)!=-9999) maskArr[i*ysz+j]=1;

	printf("Expanding");

	for(it=0;it<nExpansionIts;it++)
	{
		for(i=1;i<(xsz-1);i++) for(j=0;j<(ysz-1);j++) if(wlIn(i,j)==-9999.)
		{
			nnzNeighbours=0;
			sumNeighbours=0.;

			if(wlIn(i-1,j-1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i-1,j-1);}
			if(wlIn(i  ,j-1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i  ,j-1);}
			if(wlIn(i+1,j-1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i+1,j-1);}
			if(wlIn(i-1,j  )!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i-1,j  );}
			if(wlIn(i+1,j  )!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i+1,j  );}
			if(wlIn(i-1,j+1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i-1,j+1);}
			if(wlIn(i  ,j+1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i  ,j+1);}
			if(wlIn(i+1,j+1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i+1,j+1);}

//			if(nnzNeighbours!=0) wlOut(i,j)=sumNeighbours/nnzNeighbours;

			if(nnzNeighbours!=0) wlOut(i,j)=(sumNeighbours-expansionSlope*1.5*dx)/nnzNeighbours;
		}

		copyArray(wlOut,wlIn);

		printf(".");
	}
	printf(" done.\nSmoothing");

	for(it=0;it<nSmoothingIts;it++)
	{
		for(i=1;i<(xsz-1);i++) for(j=0;j<(ysz-1);j++) if(wlIn(i,j)!=-9999. && maskArr[i*ysz+j]==0)
		{
			nnzNeighbours=0;
			sumNeighbours=0.;

			if(wlIn(i-1,j-1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i-1,j-1);}
			if(wlIn(i  ,j-1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i  ,j-1);}
			if(wlIn(i+1,j-1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i+1,j-1);}
			if(wlIn(i-1,j  )!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i-1,j  );}
			if(wlIn(i+1,j  )!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i+1,j  );}
			if(wlIn(i-1,j+1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i-1,j+1);}
			if(wlIn(i  ,j+1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i  ,j+1);}
			if(wlIn(i+1,j+1)!=-9999) {nnzNeighbours++;sumNeighbours+=wlIn(i+1,j+1);}

			if(nnzNeighbours!=0) wlOut(i,j)=sumNeighbours/nnzNeighbours;
		}

		copyArray(wlOut,wlIn);

		printf(".");
	}
	printf(" done.\n");

    delete[] maskArr;

	return;
}

extern "C" void burnFlowPaths(float *depthArg,float *fpArg,int xsz, int ysz)
{
    int i,j;

	array2d depth, fp;

	depth.assign(xsz,ysz,depthArg);
	fp.assign(xsz,ysz,fpArg);

    for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) if(depth(i,j)==0 && fp(i,j)>0) depth(i,j)=fp(i,j);

    return;
}

extern "C" void makeWlGrid(float *dtmArg,float *depthArg,float *wlArg,int xsz,int ysz,float depthThreshold)
{
    int i,j;

	array2d dtm, depth, wl;

	dtm.assign(xsz,ysz,dtmArg);
	depth.assign(xsz,ysz,depthArg);
	wl.assign(xsz,ysz,wlArg);

    for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
    if(depth(i,j)>=depthThreshold)
        wl(i,j)=dtm(i,j)+depth(i,j);
    else wl(i,j)=-9999;

    return;
}

extern "C" void clipZero(float *inArg,int xsz,int ysz)
{
    int i,j;

	array2d in;

	in.assign(xsz,ysz,inArg);

    for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) if(in(i,j)<0) in(i,j)=0.;

    return;
}


extern "C" void resample3(float *wlArg, float *volGridArg,
	float *qxArg, float *qyArg,
     float *zMinArg, float *zMaxArg,
    float flowThresh,
	float xll, float yll, float dx,
	int xsz, int ysz,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll, float *wlOutArg, float *depthOutArg)
{
	int i,j,xi0,yi0,xi1,yi1,ii,jj;
	array2d wl, dtm, wlOut, depthOut, volGrid, qx, qy, zMin, zMax;
	float xCell1,yCell1,xCell2,yCell2,xc,yc,z1,z2,z3,z4,xi,yi,wli;
    float q1, q2, q3, q4;

// Set up array2d types
	wl.assign(xsz,ysz,wlArg);
	volGrid.assign(xsz,ysz,volGridArg);
	dtm.assign(dtmXsz,dtmYsz,dtmArg);
	wlOut.assign(dtmXsz,dtmYsz,wlOutArg);
	depthOut.assign(dtmXsz,dtmYsz,depthOutArg);
	qx.assign(xsz+1,ysz,qxArg);
	qy.assign(xsz,ysz+1,qyArg);
     zMin.assign(xsz,ysz,zMinArg);
     zMax.assign(xsz,ysz,zMaxArg);

//    printf("zMin:\n");
//    for(i=0;i<xsz;i++) for(j=0;j<ysz;j++) printf(" %0.3f",zMin(i,j));


	for(i=0;i<dtmXsz;i++) for(j=0;j<dtmYsz;j++)
	{
		depthOut(i,j)=0.;
		wlOut(i,j)=-9999.;
	}

	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		if(volGrid(i,j)<1e-3) continue;

            q1=qx(i+1,j);
            q2=qy(i,j+1);
            q3=qx(i,j);
            q4=qy(i,j);

//            if(i==19 && j==41)
//            {
//                printf("wl=%f zMin=%f zMax=%f \n",wl(i,j), zMin(i,j), zMax(i,j));
//            }

// Criteria for filling a cell are:
//      It's got flow in or out above threshold
//      Or ... WL is above max ground level in cell
//      Or ... Depth is above a threshold

//            printf("zMin=%f wl=%f \n",zMin(i,j),wl(i,j));
            if( //(wl(i,j)-zMin(i,j))<0.5) continue;

                !((fabs(q1)>=flowThresh || fabs(q2)>=flowThresh || fabs(q3)>=flowThresh || fabs(q4)>=flowThresh)||
                ((wl(i,j)-zMin(i,j))>0.5)||
                (wl(i,j)>zMax(i,j))))
                continue;

		xCell1=xll+i*dx;
		yCell1=yll+j*dx;
		xCell2=xCell1+dx;
		yCell2=yCell1+dx;

        xc=xCell1+dx/2;
        yc=yCell1+dx/2;

		if(i<(xsz-1) && q1<-flowThresh) z1=wl(i+1,j);
		else z1=wl(i,j);

		if(j<(ysz-1) && q2<-flowThresh) z2=wl(i,j+1);
		else z2=wl(i,j);

		if(i>0 && q3>flowThresh) z3=wl(i-1,j);
		else z3=wl(i,j);

		if(j>0 && q4>flowThresh) z4=wl(i,j-1);
		else z4=wl(i,j);

/*		printf("%i %i %f %f %f\n",i,j,qy(i,j),wl(i,j),wl(i,j-1));
		z1=z2=z3=wl(i,j);*/

        xi0=nint((xCell1-dtmXll)/dtmCellSize);
        xi1=nint((xCell2-dtmXll)/dtmCellSize);
        yi0=nint((yCell1-dtmYll)/dtmCellSize);
        yi1=nint((yCell2-dtmYll)/dtmCellSize);

//        if(xi0<0 || xi1>=dtmXsz || yi0<0 || yi1>=dtmYsz) continue; // Some of cell outside DTM

		for(ii=xi0;ii<xi1;ii++) for(jj=yi0;jj<yi1;jj++)
		{
                if(ii<0 || ii>=dtmXsz || jj<0 || jj>=dtmYsz) continue; // Cell outside DTM


			xi=(dtmXll+ii*dtmCellSize-xc)/(0.5*dx);
			yi=(dtmYll+jj*dtmCellSize-yc)/(0.5*dx);

			wli=bilinear(xi,yi,z1,z2,z3,z4);

			if(wli>dtm(ii,jj))
			{
				depthOut(ii,jj)=MAX(depthOut(ii,jj),wli-dtm(ii,jj));
				wlOut(ii,jj)=MAX(wlOut(ii,jj),wli);
			}
//			wlOut(ii,jj)=MAX(wlOut(ii,jj),wli);
		}
	}

//	printf("%f   %f   %f\n",dtm(1234,5678),wlOut(1234,5678),depthOut(1234,5678));

    return;
}




extern "C" void resample2(float *wlArg, float *volGridArg,
	float xll, float yll, float dx,
	int xsz, int ysz,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll, float *wlOutArg, float *depthOutArg)
{
	int i,j,xi0,yi0,xi1,yi1,ii,jj;
	array2d wl, dtm, wlOut, depthOut, volGrid;
	float xCell1,yCell1,xCell2,yCell2;

// Set up array2d types
	wl.assign(xsz,ysz,wlArg);
	volGrid.assign(xsz,ysz,volGridArg);
	dtm.assign(dtmXsz,dtmYsz,dtmArg);
	wlOut.assign(dtmXsz,dtmYsz,wlOutArg);
	depthOut.assign(dtmXsz,dtmYsz,depthOutArg);

	for(i=0;i<dtmXsz;i++) for(j=0;j<dtmYsz;j++)
	{
		depthOut(i,j)=0.;
		wlOut(i,j)=-9999.;
	}

	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		if(volGrid(i,j)<1e-3) continue;

		xCell1=xll+i*dx;
		yCell1=yll+j*dx;
		xCell2=xCell1+dx;
		yCell2=yCell1+dx;

        xi0=nint((xCell1-dtmXll)/dtmCellSize);
        xi1=nint((xCell2-dtmXll)/dtmCellSize);
        yi0=nint((yCell1-dtmYll)/dtmCellSize);
        yi1=nint((yCell2-dtmYll)/dtmCellSize);

/*        xi0=int((xCell1-dtmXll)/dtmCellSize);
        xi1=int((xCell2-dtmXll)/dtmCellSize);
        yi0=int((yCell1-dtmYll)/dtmCellSize);
        yi1=int((yCell2-dtmYll)/dtmCellSize);*/


		for(ii=xi0;ii<xi1;ii++) for(jj=yi0;jj<yi1;jj++)
		{
			if(wl(i,j)>dtm(ii,jj))
			{
				depthOut(ii,jj)=MAX(depthOut(ii,jj),wl(i,j)-dtm(ii,jj));
				wlOut(ii,jj)=MAX(wlOut(ii,jj),wl(i,j));
			}
		}
	}
    return;
}





extern "C" void findLowPoints(
	float xll, float yll, float dx,
	int xsz, int ysz,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll,
	float *lpx, float *lpy, float *lpc)
{
	int i,j,xi0,yi0,xi1,yi1,ii,jj;
	float xCell1,yCell1,xCell2,yCell2,minz;
	array2d dtm;

	dtm.assign(dtmXsz,dtmYsz,dtmArg);

	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		xCell1=xll+i*dx;
		yCell1=yll+j*dx;
		xCell2=xCell1+dx;
		yCell2=yCell1+dx;

        xi0=int((xCell1-dtmXll)/dtmCellSize);
        xi1=int((xCell2-dtmXll)/dtmCellSize);
        yi0=int((yCell1-dtmYll)/dtmCellSize);
        yi1=int((yCell2-dtmYll)/dtmCellSize);

		minz=1e20;

        if(xi0<0 || xi1>=dtmXsz || yi0<0 || yi1>=dtmYsz)
        {
        	lpc[(i*ysz+j)*3]=(xCell1+xCell2)/2;
			lpc[(i*ysz+j)*3+1]=(yCell1+yCell2)/2;

			lpx[(i*ysz+j)*3]=xCell1;
			lpx[(i*ysz+j)*3+1]=(yCell1+yCell2)/2;
			lpx[(i*ysz+j)*3+2]=-9999;

			lpy[(i*(ysz+1)+j)*3]=yCell1;
			lpy[(i*(ysz+1)+j)*3+1]=(xCell1+xCell2)/2;
			lpy[(i*(ysz+1)+j)*3+2]=-9999;

			continue;
		}

		for(ii=xi0;ii<=xi1;ii++) for(jj=yi0;jj<=yi1;jj++)
		{
			if(dtm(ii,jj)<minz)
			{
				lpc[(i*ysz+j)*3]=dtmXll+ii*dtmCellSize+dtmCellSize/2;
				lpc[(i*ysz+j)*3+1]=dtmYll+jj*dtmCellSize+dtmCellSize/2;
				minz=dtm(ii,jj);
				lpc[(i*ysz+j)*3+2]=minz;
			}
		}

		minz=1e20;
		for(jj=yi0;jj<=yi1;jj++)
		{
			if(dtm(xi0,jj)<minz)
			{
				minz=dtm(xi0,jj);
				lpx[(i*ysz+j)*3]=dtmXll+xi0*dtmCellSize+dtmCellSize/2;
				lpx[(i*ysz+j)*3+1]=dtmYll+jj*dtmCellSize+dtmCellSize/2;
				lpx[(i*ysz+j)*3+2]=minz;
			}
		}

		minz=1e20;
		for(ii=xi0;ii<=xi1;ii++)
		{
			if(dtm(ii,yi0)<minz)
			{
				minz=dtm(ii,yi0);
				lpy[(i*(ysz+1)+j)*3]=dtmXll+ii*dtmCellSize+dtmCellSize/2;
				lpy[(i*(ysz+1)+j)*3+1]=dtmYll+yi0*dtmCellSize+dtmCellSize/2;
				lpy[(i*(ysz+1)+j)*3+2]=minz;
			}
		}
// If top or right hand row/col, also do edges
		if(i==(xsz-1))
		{
			minz=1e20;
			for(jj=yi0;jj<=yi1;jj++)
			{
				if(dtm(xi1,jj)<minz)
				{
					minz=dtm(xi1,jj);
					lpx[((i+1)*ysz+j)*3]=dtmXll+xi1*dtmCellSize+dtmCellSize/2;
					lpx[((i+1)*ysz+j)*3+1]=dtmYll+jj*dtmCellSize+dtmCellSize/2;
					lpx[((i+1)*ysz+j)*3+2]=minz;
				}
			}
		}
		if(j==(ysz-1))
		{
			minz=1e20;
			for(ii=xi0;ii<=xi1;ii++)
			{
				if(dtm(ii,yi1)<minz)
				{
					minz=dtm(ii,yi1);
					lpy[(i*(ysz+1)+j+1)*3]=dtmXll+ii*dtmCellSize+dtmCellSize/2;
					lpy[(i*(ysz+1)+j+1)*3+1]=dtmYll+yi1*dtmCellSize+dtmCellSize/2;
					lpy[(i*(ysz+1)+j+1)*3+2]=minz;
				}
			}
		}

	}
    return;
}

inline float distance(float x1, float y1, float x2, float y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

void addHalfwayPointX(float &xh,float &yh,
	float xc,float yc,float x,float y,
	float yCell1, float yCell2,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll)
{
	int xih,yi1,yi2,k;
	float minZ;
	array2d dtm;

	dtm.assign(dtmXsz,dtmYsz,dtmArg);

	xh=0.5*(xc+x);
	yh=0.5*(yc+y);

	xih=int((xh-dtmXll)/dtmCellSize);

//	yi1=int((yCell1-dtmYll)/dtmCellSize);
//    yi2=int((yCell2-dtmYll)/dtmCellSize);

	yi1=int((MIN(y,yc)-dtmYll)/dtmCellSize);
    yi2=int((MAX(y,yc)-dtmYll)/dtmCellSize);


    minZ=1e20;
    for(k=yi1;k<=yi2;k++) if(dtm(xih,k)<minZ)
    {
	   	xh=dtmXll+xih*dtmCellSize;
       	yh=dtmYll+k*dtmCellSize;
   		minZ=dtm(xih,k);
    }

    return;
}

void addHalfwayPointY(float &xh,float &yh,
	float xc,float yc,float x,float y,
	float xCell1, float xCell2,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll)
{
	int yih,xi1,xi2,k;
	float minZ;
	array2d dtm;

	dtm.assign(dtmXsz,dtmYsz,dtmArg);

	xh=0.5*(xc+x);
	yh=0.5*(yc+y);

//	xi1=int((xCell1-dtmXll)/dtmCellSize);
//	xi2=int((xCell2-dtmYll)/dtmCellSize);

	xi1=int((MIN(x,xc)-dtmXll)/dtmCellSize);
	xi2=int((MAX(x,xc)-dtmXll)/dtmCellSize);

    yih=int((yh-dtmYll)/dtmCellSize);

    minZ=1e20;
    for(k=xi1;k<=xi2;k++) if(dtm(k,yih)<minZ)
    {
	   	xh=dtmXll+k*dtmCellSize;
       	yh=dtmYll+yih*dtmCellSize;
   		minZ=dtm(k,yih);
    }

    return;
}

extern "C" void writeLowPointsFile(const char *fileName,
	float *lpx, float *lpy, float *lpc, float *qx, float *qy,
	int xsz, int ysz, float dx, float xll, float yll,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll,int option)
{
	FILE *fptr;
	int i,j;
	float xc, yc, x1, x2, x3, x4, y1, y2, y3, y4, q1, q2, q3, q4;
	float xCell1, xCell2, yCell1, yCell2, xh, yh;

	fptr=fopen(fileName,"w");

	fprintf(fptr,"WKT;Q\n");

	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		xCell1=xll+i*dx;
		yCell1=yll+j*dx;
		xCell2=xCell1+dx;
		yCell2=yCell1+dx;

		xc=lpc[(i*ysz+j)*3];
		yc=lpc[(i*ysz+j)*3+1];

        x1=lpx[(i*ysz+j)*3];
        y1=lpx[(i*ysz+j)*3+1];
        x2=lpy[(i*(ysz+1)+j)*3];
        y2=lpy[(i*(ysz+1)+j)*3+1];
        x3=lpx[((i+1)*ysz+j)*3];
        y3=lpx[((i+1)*ysz+j)*3+1];
        x4=lpy[(i*(ysz+1)+j+1)*3];
        y4=lpy[(i*(ysz+1)+j+1)*3+1];

		q1=qx[i*(ysz)+j];
		q2=qy[i*(ysz+1)+j];
		q3=qx[(i+1)*(ysz)+j];
		q4=qy[i*(ysz+1)+j+1];

		if(option==1)
		{
			// Add halfway kinks
			if(distance(xc,yc,x1,y1)>dx/3)
			{
				addHalfwayPointX(xh,yh,xc,yc,x1,y1,
					yCell1,yCell2,dtmArg,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll);
				fprintf(fptr,"LINESTRING(%f %f, %f %f, %f %f);%f\n",
					x1,y1,xh,yh,xc,yc,fabs(q1));
			}
			else fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x1,y1,xc,yc,fabs(q1));

			if(distance(xc,yc,x3,y3)>dx/3)
			{
				addHalfwayPointX(xh,yh,xc,yc,x3,y3,
					yCell1,yCell2,dtmArg,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll);

				fprintf(fptr,"LINESTRING(%f %f, %f %f, %f %f);%f\n",
					x3,y3,xh,yh,xc,yc,fabs(q1));
			}
			else fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x3,y3,xc,yc,fabs(q3));

			if(distance(xc,yc,x2,y2)>dx/3)
			{
				addHalfwayPointY(xh,yh,xc,yc,x2,y2,
					xCell1,xCell2,dtmArg,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll);

				fprintf(fptr,"LINESTRING(%f %f, %f %f, %f %f);%f\n",
					x2,y2,xh,yh,xc,yc,fabs(q1));
			}
			else fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x2,y2,xc,yc,fabs(q2));

			if(distance(xc,yc,x4,y4)>dx/3)
			{
				addHalfwayPointY(xh,yh,xc,yc,x4,y4,
					xCell1,xCell2,dtmArg,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll);
				fprintf(fptr,"LINESTRING(%f %f, %f %f, %f %f);%f\n",
					x4,y4,xh,yh,xc,yc,fabs(q1));
			}
			else fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x4,y4,xc,yc,fabs(q4));
		}
		else
		{
			fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x1,y1,xc,yc,fabs(q1));
			fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x3,y3,xc,yc,fabs(q3));
			fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x2,y2,xc,yc,fabs(q2));
			fprintf(fptr,"LINESTRING(%f %f, %f %f);%f\n",x4,y4,xc,yc,fabs(q4));
		}
	}

	fclose(fptr);
	return;
}



extern "C" void flowPaths(float xll, float yll, float dx,
	int xsz, int ysz,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll,
	const char *fileName,
	float *qx, float *qy,
	int option)
{
	int i,j,xi0,yi0,xi1,yi1,ii,jj;
	float xCell1,yCell1,xCell2,yCell2,minz;
	float *lpx, *lpy, *lpc;

	array2d dtm;

	dtm.assign(dtmXsz,dtmYsz,dtmArg);

	lpx=new float[(xsz+1)*ysz*3];
	lpy=new float[xsz*(ysz+1)*3];
	lpc=new float[xsz*ysz*3];

	findLowPoints(xll,yll,dx,xsz,ysz,
		dtmArg,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll,
		lpx,lpy,lpc);

	writeLowPointsFile(fileName,lpx,lpy,lpc,qx,qy,xsz,ysz,dx,xll,yll,
		dtmArg,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll,option);


	delete[] lpx;
	delete[] lpy;
	delete[] lpc;

	return;
}





extern "C" void findLowPointsInt(
	float xll, float yll, float dx,
	int xsz, int ysz,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll,
	int *lpxi, int *lpyi, int *lpci)
{
	int i,j,xi0,yi0,xi1,yi1,ii,jj;
	float xCell0,yCell0,xCell1,yCell1,minz;
	array2d dtm;

	dtm.assign(dtmXsz,dtmYsz,dtmArg);

	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		xCell0=xll+i*dx;
		yCell0=yll+j*dx;
		xCell1=xCell0+dx;
		yCell1=yCell0+dx;

        xi0=int((xCell0-dtmXll)/dtmCellSize);
        xi1=int((xCell1-dtmXll)/dtmCellSize);
        yi0=int((yCell0-dtmYll)/dtmCellSize);
        yi1=int((yCell1-dtmYll)/dtmCellSize);

		minz=1e20;


		// Cell completely outside DTM - skip on to next cell
        if(xi0<0 || xi1>=dtmXsz || yi0<0 || yi1>=dtmYsz)
//        if(xi0>(dtmXsz-1) || xi1<0 || yi0>(dtmYsz-1) || yi1<0)
        {
        	lpci[(i*ysz+j)*2]=-9999; //(xi0+xi1)/2;
			lpci[(i*ysz+j)*2+1]=-9999; //(yi0+yi1)/2;


			lpxi[(i*ysz+j)*2]=-9999; //xi0;
			lpxi[(i*ysz+j)*2+1]=-9999; //(yi0+yi1)/2;

			lpyi[(i*(ysz+1)+j)*2]=-9999; //yi0;
			lpyi[(i*(ysz+1)+j)*2+1]=-9999; //(xi0+xi1)/2;

			continue;
		}

		// Low point within cell
		for(ii=xi0;ii<=xi1;ii++) for(jj=yi0;jj<=yi1;jj++)
		{
//		    if(ii<0 || ii>=dtmXsz || jj<0 || jj>=dtmYsz) continue;

			if(dtm(ii,jj)<minz)
			{
				lpci[(i*ysz+j)*2]=ii;
				lpci[(i*ysz+j)*2+1]=jj;
				minz=dtm(ii,jj);
			}
		}


		// ALong cell edges
		minz=1e20;
		for(jj=yi0;jj<=yi1;jj++)
		{
//		    if(xi0<0 || xi0>=dtmXsz || jj<0 || jj>=dtmYsz) continue;

			if(dtm(xi0,jj)<minz)
			{
				minz=dtm(xi0,jj);
				lpxi[(i*ysz+j)*2]=xi0;
				lpxi[(i*ysz+j)*2+1]=jj;
			}
		}

		minz=1e20;
		for(ii=xi0;ii<=xi1;ii++)
		{
//		    if(ii<0 || ii>=dtmXsz || yi0<0 || yi0>=dtmYsz) continue;

			if(dtm(ii,yi0)<minz)
			{
				minz=dtm(ii,yi0);
				lpyi[(i*(ysz+1)+j)*2]=ii;
				lpyi[(i*(ysz+1)+j)*2+1]=yi0;
			}
		}

// If top or right hand row/col, also do edges
		if(i==(xsz-1))
		{
			minz=1e20;
			for(jj=yi0;jj<=yi1;jj++)
			{
//				if(xi1<0 || xi1>=dtmXsz || jj<0 || jj>=dtmYsz) continue;

				if(dtm(xi1,jj)<minz)
				{
					minz=dtm(xi1,jj);
					lpxi[((i+1)*ysz+j)*2]=xi1;
					lpxi[((i+1)*ysz+j)*2+1]=jj;
				}
			}

		}
		if(j==(ysz-1))
		{
			minz=1e20;
			for(ii=xi0;ii<=xi1;ii++)
			{
//				if(ii<0 || ii>=dtmXsz || yi1<0 || yi1>=dtmYsz) continue;

				if(dtm(ii,yi1)<minz)
				{
					minz=dtm(ii,yi1);
					lpyi[(i*(ysz+1)+j+1)*2]=ii;
					lpyi[(i*(ysz+1)+j+1)*2+1]=yi1;
				}
			}

		}

	}
    return;
}

int lazyFlowPath(float *dtm,int xsz,int ysz,float *depthGrid,
	int i0,int i1,int j0,int j1,int iStart,int jStart,int iEnd,int jEnd,
	float defaultDepth)
{
	int iCurrent, jCurrent, i, j, iNew, jNew;
	int jump=0;

	int xWinSz=(i1-i0+1);
	int yWinSz=(j1-j0+1);

	// Can probably do this more efficiently, but this probably just grabs and
	// returns the same memory from the heap, which should be quick(ish).
	char *tmpArr=new char[xWinSz*yWinSz];

	for(i=0;i<xWinSz*yWinSz;i++) tmpArr[i]=(char)0;

	tmpArr[(iStart-i0)*yWinSz+jStart-j0]=(char)1;
//	tmpArr[(iEnd-i0)*yWinSz+jEnd-j0]=(char)1;


	iCurrent=iStart;	// Do we need iStart again?
	jCurrent=jStart;

	float minZ;

	while(iCurrent!=iEnd && jCurrent!=jEnd)
	{
		minZ=1e20;

		iNew=iCurrent;
		jNew=jCurrent;

		if(iCurrent>i0 && dtm[(iCurrent-1)*ysz+jCurrent]<minZ && tmpArr[(iCurrent-1-i0)*yWinSz+jCurrent-j0]!=1)
		{
			iNew=iCurrent-1;
			minZ=dtm[(iCurrent-1)*ysz+jCurrent];
		}
		if(iCurrent<i1 && dtm[(iCurrent+1)*ysz+jCurrent]<minZ && tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent-j0]!=1)
		{
			iNew=iCurrent+1;
			minZ=dtm[(iCurrent+1)*ysz+jCurrent];
		}
		if(jCurrent>j0 && dtm[iCurrent*ysz+jCurrent-1]<minZ && tmpArr[(iCurrent-i0)*yWinSz+jCurrent-1-j0]!=1)
		{
			jNew=jCurrent-1;
			minZ=dtm[iCurrent*ysz+jCurrent-1];
		}
		if(jCurrent<j1 && dtm[iCurrent*ysz+jCurrent+1]<minZ && tmpArr[(iCurrent-i0)*yWinSz+jCurrent+1-j0]!=1)
		{
			jNew=jCurrent+1;
			minZ=dtm[iCurrent*ysz+jCurrent+1];
		}

		if(minZ>1e10)  // Painted into a corner - jump out
		{
//			printf("In a corner - jumping!\n");
			jump=1;
			break;
		}

		tmpArr[(iNew-i0)*yWinSz+jNew-j0]=(char)1;

		iCurrent=iNew;
		jCurrent=jNew;
	}

	for(i=i0;i<=i1;i++) for(j=j0;j<=j1;j++)
		if(tmpArr[(i-i0)*yWinSz+j-j0]==1 && depthGrid[i*ysz+j]==0.)
			depthGrid[i*ysz+j]=defaultDepth;

//	if(jump==1) depthGrid[iCurrent*ysz+jCurrent]=10.*defaultDepth;


	delete[] tmpArr;

	return jump;
}

int lazyFlowPath8(float *dtm,int xsz,int ysz,float *depthGrid,
	int i0,int i1,int j0,int j1,int iStart,int jStart,int iEnd,int jEnd,
	float defaultDepth)
{
	int iCurrent, jCurrent, i, j, iNew, jNew;
	int jump=0;

	int xWinSz=(i1-i0+1);
	int yWinSz=(j1-j0+1);


	// Can probably do this more efficiently, but this probably just grabs and
	// returns the same memory from the heap, which should be quick(ish).
	char *tmpArr=new char[xWinSz*yWinSz];

	for(i=0;i<xWinSz*yWinSz;i++) tmpArr[i]=(char)0;

	tmpArr[(iStart-i0)*yWinSz+jStart-j0]=(char)1;


	iCurrent=iStart;	// Do we need iStart again?
	jCurrent=jStart;

	float minZ;

	while(iCurrent!=iEnd || jCurrent!=jEnd)
	{
		minZ=1e20;

		iNew=iCurrent;
		jNew=jCurrent;

		// Diagonal neighbours
		// NE
		if(iCurrent<i1 && jCurrent<j1 &&
			dtm[(iCurrent+1)*ysz+jCurrent+1]<minZ && tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent+1-j0]!=1)
		{
			iNew=iCurrent+1;
			jNew=jCurrent+1;
			minZ=dtm[(iCurrent+1)*ysz+jCurrent+1];
		}

		// SE
		if(iCurrent<i1 && jCurrent>j0 &&
			dtm[(iCurrent+1)*ysz+jCurrent-1]<minZ && tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent-1-j0]!=1)
		{
			iNew=iCurrent+1;
			jNew=jCurrent-1;
			minZ=dtm[(iCurrent+1)*ysz+jCurrent-1];
		}

		// NW
		if(iCurrent>i0 && jCurrent<j1 &&
			dtm[(iCurrent-1)*ysz+jCurrent+1]<minZ && tmpArr[(iCurrent-1-i0)*yWinSz+jCurrent+1-j0]!=1)
		{
			iNew=iCurrent-1;
			jNew=jCurrent+1;
			minZ=dtm[(iCurrent-1)*ysz+jCurrent+1];
		}

		// SW
		if(iCurrent>i0 && jCurrent>j0 &&
			dtm[(iCurrent-1)*ysz+jCurrent-1]<minZ && tmpArr[(iCurrent-1-i0)*yWinSz+jCurrent-1-j0]!=1)
		{
			iNew=iCurrent-1;
			jNew=jCurrent-1;
			minZ=dtm[(iCurrent-1)*ysz+jCurrent-1];
		}

		// Four neighbours
		// W
		if(iCurrent>i0 && dtm[(iCurrent-1)*ysz+jCurrent]<minZ && tmpArr[(iCurrent-1-i0)*yWinSz+jCurrent-j0]!=1)
		{
			iNew=iCurrent-1;
			jNew=jCurrent;
			minZ=dtm[(iCurrent-1)*ysz+jCurrent];
		}
		// E
		if(iCurrent<i1 && dtm[(iCurrent+1)*ysz+jCurrent]<minZ && tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent-j0]!=1)
		{
			iNew=iCurrent+1;
			jNew=jCurrent;
			minZ=dtm[(iCurrent+1)*ysz+jCurrent];
		}
		// S
		if(jCurrent>j0 && dtm[iCurrent*ysz+jCurrent-1]<minZ && tmpArr[(iCurrent-i0)*yWinSz+jCurrent-1-j0]!=1)
		{
			iNew=iCurrent;
			jNew=jCurrent-1;
			minZ=dtm[iCurrent*ysz+jCurrent-1];
		}
		// N
		if(jCurrent<j1 && dtm[iCurrent*ysz+jCurrent+1]<minZ && tmpArr[(iCurrent-i0)*yWinSz+jCurrent+1-j0]!=1)
		{
			iNew=iCurrent;
			jNew=jCurrent+1;
			minZ=dtm[iCurrent*ysz+jCurrent+1];
		}



		if(minZ>1e10)  // Painted into a corner - jump out
		{
//			printf("In a corner - jumping!\n");
			jump=1;
			break;
		}

		tmpArr[(iNew-i0)*yWinSz+jNew-j0]=(char)1;

		iCurrent=iNew;
		jCurrent=jNew;
	}

	for(i=i0;i<=i1;i++) for(j=j0;j<=j1;j++)
		if(tmpArr[(i-i0)*yWinSz+j-j0]==1 && depthGrid[i*ysz+j]==0.)
			depthGrid[i*ysz+j]=defaultDepth;

//	if(jump==1) depthGrid[iCurrent*ysz+jCurrent]=10.*defaultDepth;

	delete[] tmpArr;

	return jump;
}

float cellScore(float zCurrent,float zNew,int iCurrent, int jCurrent, int iNew, int jNew, int iTarget, int jTarget, float dWeight)
{
	float dCurrent=sqrt((iCurrent-iTarget)*(iCurrent-iTarget)+(jCurrent-jTarget)*(jCurrent-jTarget));
	float dNew=sqrt((iNew-iTarget)*(iNew-iTarget)+(jNew-jTarget)*(jNew-jTarget));

	return (1.-dWeight)*(zNew-zCurrent)+dWeight*(dNew-dCurrent);
}

int homingPigeon(float *dtm,int xsz,int ysz,float *depthGrid,
	int i0,int i1,int j0,int j1,int iStart,int jStart,int iEnd,int jEnd,
	float defaultDepth, float dWeight)
{
	int iCurrent, jCurrent, i, j, iNew, jNew;
	int jump=0;

	int xWinSz=(i1-i0+1);
	int yWinSz=(j1-j0+1);


	// Can probably do this more efficiently, but this probably just grabs and
	// returns the same memory from the heap, which should be quick(ish).
	char *tmpArr=new char[xWinSz*yWinSz];

	for(i=0;i<xWinSz*yWinSz;i++) tmpArr[i]=(char)0;

	tmpArr[(iStart-i0)*yWinSz+jStart-j0]=(char)1;


	iCurrent=iStart;	// Do we need iStart again?
	jCurrent=jStart;

	float minScore,zCurrent;

	while(iCurrent!=iEnd || jCurrent!=jEnd)
	{
		minScore=1e20;

		iNew=iCurrent;
		jNew=jCurrent;

		zCurrent=dtm[iCurrent*ysz+jCurrent];

		// Diagonal neighbours
		// NE
		if(iCurrent<i1 && jCurrent<j1 &&
			cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent+1],iCurrent,jCurrent,iCurrent+1,jCurrent+1,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent+1-j0]!=1)
		{
			iNew=iCurrent+1;
			jNew=jCurrent+1;
			minScore=cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent+1],iCurrent,jCurrent,iCurrent+1,jCurrent+1,iEnd,jEnd,dWeight);
		}

		// SE
		if(iCurrent<i1 && jCurrent>j0 &&
			cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent-1],iCurrent,jCurrent,iCurrent+1,jCurrent-1,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent-1-j0]!=1)
		{
			iNew=iCurrent+1;
			jNew=jCurrent-1;
			minScore=cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent-1],iCurrent,jCurrent,iCurrent+1,jCurrent-1,iEnd,jEnd,dWeight);
		}

		// NW
		if(iCurrent>i0 && jCurrent<j1 &&
			cellScore(zCurrent,dtm[(iCurrent-1)*ysz+jCurrent+1],iCurrent,jCurrent,iCurrent-1,jCurrent+1,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent-1-i0)*yWinSz+jCurrent+1-j0]!=1)
		{
			iNew=iCurrent-1;
			jNew=jCurrent+1;
			minScore=cellScore(zCurrent,dtm[(iCurrent-1)*ysz+jCurrent+1],iCurrent,jCurrent,iCurrent-1,jCurrent+1,iEnd,jEnd,dWeight);
		}

		// SW
		if(iCurrent>i0 && jCurrent>j0 &&
			cellScore(zCurrent,dtm[(iCurrent-1)*ysz+jCurrent-1],iCurrent,jCurrent,iCurrent-1,jCurrent-1,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent-1-i0)*yWinSz+jCurrent-1-j0]!=1)
		{
			iNew=iCurrent-1;
			jNew=jCurrent-1;
			minScore=cellScore(zCurrent,dtm[(iCurrent-1)*ysz+jCurrent-1],iCurrent,jCurrent,iCurrent-1,jCurrent-1,iEnd,jEnd,dWeight);
		}

		// Four neighbours
		// W
		if(iCurrent>i0 &&
			cellScore(zCurrent,dtm[(iCurrent-1)*ysz+jCurrent],iCurrent,jCurrent,iCurrent-1,jCurrent,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent-1-i0)*yWinSz+jCurrent-j0]!=1)
		{
			iNew=iCurrent-1;
			jNew=jCurrent;
			minScore=cellScore(zCurrent,dtm[(iCurrent-1)*ysz+jCurrent],iCurrent,jCurrent,iCurrent-1,jCurrent,iEnd,jEnd,dWeight);
		}
		// E
		if(iCurrent<i1 &&
			cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent],iCurrent,jCurrent,iCurrent+1,jCurrent,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent-j0]!=1)
		{
			iNew=iCurrent+1;
			jNew=jCurrent;
			minScore=cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent],iCurrent,jCurrent,iCurrent+1,jCurrent,iEnd,jEnd,dWeight);
		}
		// S
		if(jCurrent>j0 &&
			cellScore(zCurrent,dtm[(iCurrent)*ysz+jCurrent-1],iCurrent,jCurrent,iCurrent,jCurrent-1,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent-i0)*yWinSz+jCurrent-1-j0]!=1)
		{
			iNew=iCurrent;
			jNew=jCurrent-1;
			minScore=cellScore(zCurrent,dtm[(iCurrent)*ysz+jCurrent-1],iCurrent,jCurrent,iCurrent,jCurrent-1,iEnd,jEnd,dWeight);
		}
		// N
		if(jCurrent<j1 &&
			cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent],iCurrent,jCurrent,iCurrent+1,jCurrent,iEnd,jEnd,dWeight)<minScore &&
			tmpArr[(iCurrent+1-i0)*yWinSz+jCurrent-j0]!=1)
		{
			iNew=iCurrent;
			jNew=jCurrent+1;
			minScore=cellScore(zCurrent,dtm[(iCurrent+1)*ysz+jCurrent],iCurrent,jCurrent,iCurrent+1,jCurrent,iEnd,jEnd,dWeight);
		}


		if(minScore>1e10)  // Painted into a corner - jump out
		{
//			printf("In a corner - jumping!\n");
			jump=1;
			break;
		}

		tmpArr[(iNew-i0)*yWinSz+jNew-j0]=(char)1;

		iCurrent=iNew;
		jCurrent=jNew;
	}

	for(i=i0;i<=i1;i++) for(j=j0;j<=j1;j++)
		if(tmpArr[(i-i0)*yWinSz+j-j0]==1 && depthGrid[i*ysz+j]==0.)
			depthGrid[i*ysz+j]=defaultDepth;

//	if(jump==1) depthGrid[iCurrent*ysz+jCurrent]=10.*defaultDepth;

	delete[] tmpArr;

	return jump;
}


extern "C" void fillLazyFlowPaths(
	int *lpxi, int *lpyi, int *lpci, float *qx, float *qy,
	int xsz, int ysz,
	float *dtmArg, int dtmXsz, int dtmYsz,
	float *depthGrid, float threshQ, float defaultDepth, float dWeight)
{
	int i,j;
	int xci, yci, i1, i2, i3, i4, j1, j2, j3, j4;
	float q1, q2, q3, q4;

	int jumpCount=0;
	int totalCount=0;


    // Use a function pointer to make it easier to change flow path algorithm
     int (*flowPathFnPointer)(float *,int,int,float *,int,int,int,int,int,int,int,int,float,float);
    flowPathFnPointer=&homingPigeon;


	for(i=0;i<xsz;i++) for(j=0;j<ysz;j++)
	{
		xci=lpci[(i*ysz+j)*2];
		yci=lpci[(i*ysz+j)*2+1];

        i1=lpxi[(i*ysz+j)*2];
        j1=lpxi[(i*ysz+j)*2+1];
        i2=lpyi[(i*(ysz+1)+j)*2];
        j2=lpyi[(i*(ysz+1)+j)*2+1];
        i3=lpxi[((i+1)*ysz+j)*2];
        j3=lpxi[((i+1)*ysz+j)*2+1];
        i4=lpyi[(i*(ysz+1)+j+1)*2];
        j4=lpyi[(i*(ysz+1)+j+1)*2+1];

		q1=qx[i*(ysz)+j];
		q2=qy[i*(ysz+1)+j];
		q3=qx[(i+1)*(ysz)+j];
		q4=qy[i*(ysz+1)+j+1];

		if(i1==-9999 || i3==-9999 || j2==-9999 || j4==-9999) continue;

		if(q1>threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,i1,j1,xci,yci,defaultDepth,dWeight);
		}
		if(q1<-threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,xci,yci,i1,j1,defaultDepth,dWeight);
		}

		if(q2>threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,i2,j2,xci,yci,defaultDepth,dWeight);
		}
		if(q2<-threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,xci,yci,i2,j2,defaultDepth,dWeight);
		}

		if(q3<-threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,i3,j3,xci,yci,defaultDepth,dWeight);
		}
		if(q3>threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,xci,yci,i3,j3,defaultDepth,dWeight);
		}

		if(q4<-threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,i4,j4,xci,yci,defaultDepth,dWeight);
		}
		if(q4>threshQ)
		{
			totalCount++;
			jumpCount+=flowPathFnPointer(dtmArg,dtmXsz,dtmYsz,depthGrid,i1,i3,j2,j4,xci,yci,i4,j4,defaultDepth,dWeight);
		}
	}

//	printf("Error count: failed to find flow path in %2f%% of cells\n",(100.*jumpCount)/(totalCount));

	return;
}



extern "C" void lazyFlowPaths(float xll, float yll, float dx,
	int xsz, int ysz,
	float *dtmArg, float dtmCellSize, int dtmXsz, int dtmYsz,
	float dtmXll, float dtmYll,
	float *depthGrid,
	float *qx, float *qy, float flowThreshold, float defaultDepth,
    float dWeight)
{
	int i,j,xi0,yi0,xi1,yi1,ii,jj;
	float xCell1,yCell1,xCell2,yCell2,minz;
	int *lpxi, *lpyi, *lpci;

	array2d dtm;

	dtm.assign(dtmXsz,dtmYsz,dtmArg);

	lpxi=new int[(xsz+1)*ysz*2];
	lpyi=new int[xsz*(ysz+1)*2];
	lpci=new int[xsz*ysz*2];

	findLowPointsInt(xll,yll,dx,xsz,ysz,
		dtmArg,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll,
		lpxi,lpyi,lpci);


	fillLazyFlowPaths(lpxi,lpyi,lpci,qx,qy,xsz,ysz,
		dtmArg,dtmXsz,dtmYsz,depthGrid,flowThreshold,defaultDepth,dWeight);


	delete[] lpxi;
	delete[] lpyi;
	delete[] lpci;

	return;
}


//-------------------------------------------------------------------------
// Take in rainfall, additional rainfall and curve number grids and calculate
// additional runoff


float scsRunoff(float rainfall,float cn)
{
    float rainfallInches=rainfall/25.4;
    float S=1000./cn-10.;

    if(rainfallInches<=0.2*S) return 0;
    else return 25.4*pow(rainfallInches-0.2*S,2)/(rainfallInches+0.8*S);
}


float scsAdditionalRunoff(float r, float dr, float cn)
{
    float S=25.4*(1000./cn-10.);

    if(r<=0.2*S) return 0;
    else
    {
        float dQdI=2.*(r-0.2*S)/(r+0.8*S)-pow((r-0.2*S)/(r+0.8*S),2);
        return dQdI*dr;
    }
}


extern "C" void scsAdditionalRunoff(float *totalRainfallArg, float *additionalRainfallArg,
    float *curveNumberArg, float *additionalRunoffArg, int xsz, int ysz)
{
    array2d totalRainfall,additionalRainfall,curveNumber,additionalRunoff;

    totalRainfall.assign(xsz,ysz,totalRainfallArg);
    additionalRainfall.assign(xsz,ysz,additionalRainfallArg);
    curveNumber.assign(xsz,ysz,curveNumberArg);
    additionalRunoff.assign(xsz,ysz,additionalRunoffArg);

    for(int i=0;i<xsz;i++) for (int j=0;j<ysz;j++)
        additionalRunoff(i,j)=scsAdditionalRunoff(totalRainfall(i,j),additionalRainfall(i,j),curveNumber(i,j));
}



