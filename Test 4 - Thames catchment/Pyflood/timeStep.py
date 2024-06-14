import numpy

def timeStep(u,v,h,unew,vnew,hnew,z,\
    ubcVals,ubcLocX,ubcLocY,ubcN,\
    vbcVals,vbcLocX,vbcLocY,vbcN,\
    hbcVals,hbcLocX,hbcLocY,hbcN,\
    nstepx,nstepy,\
    totalTime,displayInterval,\
    calcTimeStep,updateBCs,display,dx,n,beta):

    newAdvectionMethod=True
    depthThresh=0.01
    g=9.81

    ts=0
    currentTime=0.
    lastDisplayTime=0.

    while True:
        dt=calcTimeStep()
        updateBCs(currentTime)

#        print h[-1,nstepy/2]

###################################################################################


    # Flow in x-direction - interior points
        for i in range(1,nstepx):
            for j in range(0,nstepy):
                if(h[i,j]>=depthThresh or h[i-1,j]>=depthThresh):
                    dhdx=(h[i,j]+z[i,j]-h[i-1,j]-z[i-1,j])/dx

#                    hEdge=0.5*(h[i,j]+h[i-1,j])
                    hEdge=max(h[i,j]+z[i,j],h[i-1,j]+z[i-1,j])-max(z[i,j],z[i-1,j])

    #                print "hEdge_x=", hEdge

                    if hEdge<=0: continue

        # Calculate velocity gradient using upwinding
                    if u[i,j]>=0:
                        dudx=(u[i,j]-u[i-1,j])/dx
                    else:
                        dudx=(u[i+1,j]-u[i,j])/dx

        # Cross terms - a bit more complicated
                    if i<nstepx-1:
                        v4=0.25*(v[i,j]+v[i,j+1]+v[i-1,j+1]+v[i-1,j])
                    else:
                        v4=0.

                    if j>0 and j<nstepy-1:
                        if v4>0:
                            dudy=(u[i,j]-u[i,j-1])/dx
                        else:
                            dudy=(u[i,j+1]-u[i,j])/dx
                    else:
                        dudy=0

                    absVel=numpy.sqrt(v4*v4+u[i,j]*u[i,j])

                    if newAdvectionMethod:
                        unew[i,j]=(u[i,j]/dt-beta*v4*dudy-g*dhdx)/(1./dt+g*n*n*absVel/(hEdge**(4./3))+beta*dudx)
                    else:
                        unew[i,j]=(u[i,j]/dt-beta*(u[i,j]*dudx+v4*dudy)-g*dhdx)/(1./dt+g*n*n*absVel/(hEdge**(4./3)))

                else:
                    unew[i,j]=0.


    # Flow on boundaries - set to zero - can be updated by BCs
        unew[0,:]=0.
        unew[-1,:]=0.
        vnew[:,0]=0.
        vnew[:,-1]=0.

    # Insert stuff for boundary conditions
        for i in range(ubcN):
            unew[ubcLocX[i],ubcLocY[i]]=ubcVals[i]

        for i in range(vbcN):
            vnew[vbcLocX[i],vbcLocY[i]]=vbcVals[i]



    # Flow in y-direction - interior points
        for i in range(0,nstepx):
            for j in range(1,nstepy):
                if(h[i,j]>=depthThresh or h[i,j-1]>=depthThresh):

                    dhdx=(h[i,j]+z[i,j]-h[i,j-1]-z[i,j-1])/dx
#                    hEdge=0.5*(h[i,j]+h[i,j-1])

                    hEdge=max(h[i,j]+z[i,j],h[i,j-1]+z[i,j-1])-max(z[i,j],z[i,j-1])

#                    if i==3:
#                        print "hEdge_y=", hEdge, h[i,j]+z[i,j],h[i,j-1]+z[i,j-1],z[i,j],z[i,j-1]

                    if hEdge<=0: continue

        # Calculate velocity gradient using upwinding
                    if v[i,j]>=0:
                        dvdy=(v[i,j]-v[i,j-1])/dx
                    else:
                        dvdy=(v[i,j+1]-v[i,j])/dx

        # Cross terms - a bit more complicated
                    if j>0:
                        u4=0.25*(u[i,j]+u[i+1,j]+u[i,j-1]+u[i+1,j-1])
                    else:
                        u4=0.

                    if i>0 and i<nstepx-1:
                        if u4>0:
                            dvdx=(v[i,j]-v[i-1,j])/dx
                        else:
                            dvdx=(v[i+1,j]-v[i,j])/dx
                    else:
                        dvdx=0

                    absVel=numpy.sqrt(u4*u4+v[i,j]*v[i,j])

                    if newAdvectionMethod:
                        vnew[i,j]=(v[i,j]/dt-beta*u4*dvdx-g*dhdx)/(1./dt+g*n*n*absVel/(hEdge**(4./3))+beta*dvdy)
                    else:
                        vnew[i,j]=(v[i,j]/dt-beta*(v[i,j]*dvdy+u4*dvdx)-g*dhdx)/(1./dt+g*n*n*absVel/(hEdge**(4./3)))
                else:
                    vnew[i,j]=0.

    # Flow in y-direction - points on upper+lower boundaries
        vnew[:,0]=0.
        vnew[:,-1]=0.

        # Check for drying points and modify flows accordingly
        for i in range(1,nstepx-1):
            for j in range(1,nstepy-1):
                if h[i,j]>depthThresh and h[i-1,j]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i-1,j]+z[i-1,j])-max(z[i,j],z[i-1,j])
                else:
                    hEdge=0.
                q1=hEdge*unew[i,j]

                if h[i,j]>depthThresh and h[i+1,j]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i+1,j]+z[i+1,j])-max(z[i,j],z[i+1,j])
                else:
                    hEdge=0.
                q2=hEdge*unew[i+1,j]

                if h[i,j]>depthThresh and h[i,j-1]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i,j-1]+z[i,j-1])-max(z[i,j],z[i,j-1])
                else:
                    hEdge=0.
                q3=hEdge*vnew[i,j]

                if h[i,j]>depthThresh and h[i,j+1]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i,j+1]+z[i,j+1])-max(z[i,j],z[i,j+1])
                else:
                    hEdge=0.
                q4=hEdge*vnew[i,j+1]

                dH=(q1-q2+q3-q4)*dt/dx
                hnewTest=h[i,j]+(q1-q2+q3-q4)*dt/dx

                if (h[i,j]+dH)<0:
                    dryFactor=-h[i,j]/dH
#                    print "dryFactor=", dryFactor
                    unew[i,j]*=dryFactor
                    unew[i+1,j]*=dryFactor
                    vnew[i,j]*=dryFactor
                    vnew[i,j+1]*=dryFactor

        # Update depths: interior points
        for i in range(1,nstepx-1):
            for j in range(1,nstepy-1):

                if h[i,j]>depthThresh or h[i-1,j]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i-1,j]+z[i-1,j])-max(z[i,j],z[i-1,j])
                else:
                    hEdge=h[i,j]
                q1=hEdge*unew[i,j]

                if h[i,j]>depthThresh or h[i+1,j]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i+1,j]+z[i+1,j])-max(z[i,j],z[i+1,j])
                else:
                    hEdge=h[i,j]
                q2=hEdge*unew[i+1,j]

                if h[i,j]>depthThresh or h[i,j-1]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i,j-1]+z[i,j-1])-max(z[i,j],z[i,j-1])
                else:
                    hEdge=h[i,j]
                q3=hEdge*vnew[i,j]

                if h[i,j]>depthThresh or h[i,j+1]>depthThresh:
                    hEdge=max(h[i,j]+z[i,j],h[i,j+1]+z[i,j+1])-max(z[i,j],z[i,j+1])
                else:
                    hEdge=h[i,j]
                q4=hEdge*vnew[i,j+1]

#                dh=(q1-q2+q3-q4)*dt/dx
#
#                if abs(q1)>0 or abs(q2)>0 or abs(q3)>0 or abs(q4)>0:
#                    print q1, q2, q3, q4
#
#                if abs(dh)>1e-3:
#                    print i,j,"dh=",dh

                hnew[i,j]+=(q1-q2+q3-q4)*dt/dx

        # Update depths: edge points
        for i in range(1,nstepx-1):
            # Bottom edge
            q1=0.5*(h[i,0]+h[i-1,0])*unew[i,0]
            q2=0.5*(h[i,0]+h[i+1,0])*unew[i+1,0]
            q3=h[i,0]*vnew[i,0]
            q4=0.5*(h[i,0]+h[i,1])*vnew[i,1]
            hnew[i,0]+=(q1-q2+q3-q4)*dt/dx

            # Top edge
            q1=0.5*(h[i,-1]+h[i-1,-1])*unew[i,-1]
            q2=0.5*(h[i,-1]+h[i+1,-1])*unew[i+1,-1]
            q3=0.5*(h[i,-1]+h[i,-2])*vnew[i,-2]
            q4=h[i,-1]*vnew[i,-1]
            hnew[i,-1]+=(q1-q2+q3-q4)*dt/dx

        for j in range(1,nstepy-1):
            # Left edge
            q1=h[0,j]*unew[0,j]

            hEdge=max(h[0,j]+z[0,j],h[1,j]+z[1,j])-max(z[0,j],z[1,j])
            hEdge=max(hEdge,0)
            q2=hEdge*unew[1,j]

            q3=0.5*(h[0,j]+h[0,j-1])*vnew[0,j]
            q4=0.5*(h[0,j]+h[0,j+1])*vnew[0,j+1]
            hnew[0,j]+=(q1-q2+q3-q4)*dt/dx

            # Right edge
            q1=0.5*(h[-1,j]+h[-2,j])*unew[-2,j]
            q2=h[-1,j]*unew[-1,j]
            q3=0.5*(h[-1,j]+h[-1,j-1])*vnew[-1,j]
            q4=0.5*(h[-1,j]+h[-1,j+1])*vnew[-1,j+1]
            hnew[-1,j]+=(q1-q2+q3-q4)*dt/dx

        # Update depths: corner points
        # SW
        q1=h[0,0]*unew[0,0]
        q2=0.5*(h[0,0]+h[1,0])*unew[1,0]
        q3=h[0,0]*vnew[0,0]
        q4=0.5*(h[0,0]+h[0,1])*vnew[0,1]
        hnew[0,0]+=(q1-q2+q3-q4)*dt/dx

        # NW
        q1=h[0,-1]*unew[0,-1]
        q2=0.5*(h[0,-1]+h[1,0])*unew[1,-1]
        q3=0.5*(h[0,-1]+h[0,-2])*vnew[0,-2]
        q4=h[0,-1]*vnew[0,-1]
        hnew[0,-1]+=(q1-q2+q3-q4)*dt/dx

        # SE
        q1=0.5*(h[-1,0]+h[-2,0])*unew[-2,0]
        q2=h[-1,0]*unew[-1,0]
        q3=h[-1,0]*vnew[-1,0]
        q4=0.5*(h[-1,0]+h[-1,1])*vnew[0,-1]
        hnew[-1,0]+=(q1-q2+q3-q4)*dt/dx

        # NE
        q1=0.5*(h[-1,-1]+h[-2,-1])*unew[-2,-1]
        q2=h[-1,-1]*unew[-1,-1]
        q3=0.5*(h[-1,-1]+h[-1,-2])*vnew[-1,-2]
        q4=h[-1,-1]*vnew[-1,-1]
        hnew[-1,-1]+=(q1-q2+q3-q4)*dt/dx

# Boundary conditions
        for i in range(hbcN):
            hnew[hbcLocX[i],hbcLocY[i]]=max(0,hbcVals[i]-z[hbcLocX[i],hbcLocY[i]])

        # Wetting points
#        for i in range(1,nstepx-1):
#            for j in range(1,nstepy-1):
#                if h[i,j]<depthThresh:
#                    wlij=z[i,j]+h[i,j]
#
#                    wlNeighbourMax=-9999
#
#                    if h[i+1,j]>2*depthThresh and (h[i+1,j]+z[i+1,j])>wlNeighbourMax:
#                        wlNeighbourMax=(h[i+1,j]+z[i+1,j])
#                        ni=i+1
#                        nj=j
#                    elif h[i-1,j]>2*depthThresh and (h[i-1,j]+z[i-1,j])>wlNeighbourMax:
#                        wlNeighbourMax=(h[i-1,j]+z[i-1,j])
#                        ni=i-1
#                        nj=j
#                    elif h[i,j+1]>2*depthThresh and (h[i,j+1]+z[i,j+1])>wlNeighbourMax:
#                        wlNeighbourMax=(h[i,j+1]+z[i,j+1])
#                        ni=i
#                        nj=j+1
#                    elif h[i,j-1]>2*depthThresh and (h[i,j-1]+z[i,j-1])>wlNeighbourMax:
#                        wlNeighbourMax=(h[i,j-1]+z[i,j-1])
#                        ni=i
#                        nj=j-1
#
#                    if wlNeighbourMax>-9999 and wlNeighbourMax>wlij: # Wet cell from highest neighbours
#                        print "Wetting cell", i, j, "from neighbour", ni, nj
#                        print h[ni,nj]
#                        hnew[i,j]=depthThresh
#                        hnew[ni,nj]-=(depthThresh-h[i,j])


###################################################################################
#        print unew.min(), unew.max()
#        print vnew.min(), vnew.max()
#        print "h[1,75]=", hnew[1,75]

        h[:]=hnew[:]
        u[:]=unew[:]
        v[:]=vnew[:]

        currentTime+=dt
        if (currentTime-lastDisplayTime)>=displayInterval:
            display(ts,currentTime)
            lastDisplayTime=displayInterval*(int)(currentTime/displayInterval)

        if currentTime>=totalTime:
            display(ts,currentTime)
            break

        ts+=1

    return