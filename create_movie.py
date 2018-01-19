#!/usr/bin/env python2.7

##Note - attempting to create a move with more than 97 frames will throw an
##UnboundLocalError. This should be a IndexError but I supressed it as in 
##other situations the code would throw and index error but otherwise work fine.
##This occurs on lines 261-264. The rest of the code in that area is sub-par
##and may need re-writing as I commented a large section of it out due to issues
##with the paramater hvel. Running this program with the paramater hvel produces
##a movie in which all frames are identical. 

import os,sys
import netCDF4
import numpy
import matplotlib.animation as animation
#import matplotlib.tri as mtri #Commented out as unused
#from matplotlib.collections import PolyCollection,LineCollection #Commented out as unused
from matplotlib.dates import YearLocator, DayLocator, DateFormatter, AutoDateLocator
import matplotlib.pyplot as plt
import json
from scipy.interpolate import griddata
#import pylab #Commented out as unused

global r, h
r=h=1

def lat_lon_proportion(plt,ax):
	#Calculate the distances along the axis
	xlim=plt.xlim()
#	ylim=plt.ylim() #Commented out as variable unused
#	x_dist = numpy.diff(xlim); #Commented out as variable unused
#	y_dist = numpy.diff(ylim); #Commented out as variable unused

	#Adjust the aspect ratio
	c_adj = numpy.cos(numpy.mean(numpy.deg2rad(xlim)));
#	dar = [1,c_adj,1]; #Commented out as variable unused
#	pbar = [x_dist[0]*c_adj/y_dist[0],1,1 ]; #Commented out as variable unused
	ax.set_aspect(abs(c_adj))
	#ax.set_aspect(abs(x_dist[0]*c_adj/y_dist[0]))

def near2d_selfe(x, y, x0, y0, nnodes=1, nreturn='multiple'):

    """
    Find the indexes of the grid point that is
    nearest a chosen (x0, y0).
    Usage: line, col = near2d(x, y, x0, y0)
    """   
    dx = numpy.abs(x - x0); dx = dx / dx.max()
    dy = numpy.abs(y - y0); dy = dy / dy.max()
    dn = dx + dy    

    if nnodes > 1:
        line = []
        for node in range(nnodes):
            fn = numpy.where(dn == dn.min())[0][0]
            line.append(fn)
            dn[fn] = 9999

    else:
        fn = numpy.where(dn == dn.min())[0][0]
        line = [int(fn)]  

    if nreturn == 'multiple':

        if nnodes > 1: return line
        else: return line[0]

    elif nreturn == 'single':
        return line[nnodes-1]


def plotpatch(bnd):

	data=numpy.loadtxt(bnd)
	x=data[:,0];
	y=data[:,1];

	splits=numpy.bitwise_or(numpy.isnan(x),x>99999999999,y==1).nonzero()[0]
	splits=numpy.insert(splits,0,-1)
	splits=numpy.insert(splits,len(splits),len(splits))
	pols=[]
	for ist in range(0,len(splits)-2):
		ind=numpy.arange(splits[ist]+1,splits[ist+1])
		plt.fill( x[ind],y[ind],'silver')

def get_min_max(dirin,params,Istart,Iend,level):
    global Z,b
    Zmin=numpy.inf
    Zmax=numpy.inf*-1
    for k in range(Istart,Iend+1):
        for vars in params:
            fullfile=os.path.join(dirin,'schout_'+str(k)+'.nc')
            nc=netCDF4.Dataset(fullfile)
            nt=len(nc.variables['time'])
          
            for t in range(0,nt):
                if 'two'  in nc.variables[vars].dimensions:
                    tmp=nc.variables[vars][t]
                    u=tmp[...,0]
                    v=tmp[...,1]
                    Z=numpy.sqrt(u**2+v**2)
                else:
                    Z=nc.variables[vars][t,:]
	    
                if len(Z.shape)>1:
                    Z=Z[:,level]
                elif len(Z.shape)==3:
                    Z=Z[:,t,level]
            if Zmin==numpy.inf:
#                z=Z #Commented out as variable unused
                pass
        Zmin=min(Zmin,min(Z))
        Zmax=max(Zmax,max(Z))
    return Z,numpy.round(Zmax,2),numpy.round(Zmin,2)

def extract_ts(Istart,Iend,node,dirin):
    E=[]
    T=[]

    for k in range(Istart,Iend+1):
        fullfile=os.path.join(dirin,'schout_'+str(k)+'.nc')
        nc=netCDF4.Dataset(fullfile)
        dtime = netCDF4.num2date(nc.variables['time'][:],nc.variables['time'].units)
        T=numpy.hstack((T,dtime))
        E=numpy.hstack((E,nc.variables['elev'][:,node]))
    return E,T



def process(moviefile,dirin,Istart,Iend,dpi,params,quiver,quiver_res,quiver_scale,fps,zmin,zmax,bnd,level,lim):
    global r,X, h    
    figdir, file = os.path.split(moviefile)

    fig = plt.figure(figsize=(30,18))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    tt1 = ax.text(.5, 1.05, '', transform = ax.transAxes, va='center',fontsize = 30)
    ax.tick_params(labelsize=15)


    
    # get the static variable
    staticfile=os.path.join(dirin,'schout_'+str(Istart)+'.nc')
    
    ncs=netCDF4.Dataset(staticfile)
    X=ncs.variables['SCHISM_hgrid_node_x'][:]
    Y=ncs.variables['SCHISM_hgrid_node_y'][:]
    Z=ncs.variables['depth'][:]
    ele=ncs.variables['SCHISM_hgrid_face_nodes'][...,:3]-1
    nt=len(ncs.variables['time'])
#    t=ncs.variables['time'] #Commented out as variable unused
    
    if quiver is not None:
        if lim is not None:
            Xreg, Yreg = numpy.meshgrid(numpy.linspace(lim[0],lim[1], quiver_res), numpy.linspace(lim[2], lim[3], quiver_res))
        else:
            Xreg, Yreg = numpy.meshgrid(numpy.linspace(min(X[:]),max(X[:]), quiver_res), numpy.linspace(min(Y[:]), max(Y[:]), quiver_res))
        XY=numpy.array((X,Y)).T

    ncs.close()
    if lim is not None:
        node   = near2d_selfe(X[:],Y[:],(lim[1]-lim[0])/2, (lim[3]-lim[2])/2,  nnodes=1, nreturn='single')
    else:
        node   = near2d_selfe(X[:],Y[:],numpy.mean(X[:]),numpy.mean(Y[:]),  nnodes=1, nreturn='single')
    
    elev,time=extract_ts(Istart,Iend,node,dirin)

    ZZ,ZZmax,ZZmin=get_min_max(dirin,params,Istart,Iend,level)
    if zmin is None:
        Zmin=ZZmin
    else:
        Zmin=zmin
    if zmax is None:
        Zmax=ZZmax
    else:
        Zmax=zmax
        
    nc = None
    Q = None

    
    def init_img():
            global r
            if r >1: ##When this function is run more than once, it causes errors.
                return() #This code prevents that.
            r+=1
            global nc,Q,F,tide
            
            ZZ[ZZ>Zmax]=Zmax
            ZZ[ZZ<Zmin]=Zmin
            levels = numpy.linspace(Zmin, Zmax, 60)
            F=plt.tricontourf(X,Y,ele,ZZ,vmin=Zmin,vmax=Zmax,cmap=plt.cm.Spectral_r,levels=levels)
            plt.clim(Zmin,Zmax)
            plt.xlabel('Easting (meters)',fontsize = 30)
            plt.ylabel('Northing (meters)',fontsize = 30)
            
            cbar=plt.colorbar(F,ax=ax)#numpy.linspace(Zmin,Zmax,10))
            if params == ['elev']: #Find out what units are used.
                cbar.set_label(r"Wave height [m]", size=30)#Elevation probably meters
            elif params == ['temp']:
                cbar.set_label(r"Water Temperature [°C]", size=30)#Temperature probably °C
            elif params == ['hvel'] or params == ['dahv']:
                cbar.set_label(r"Current speed [m.s^-1]", size=30)#3D & 2D water current probably m.s^-1
            else:
                cbar.set_label(r"Current speed [m.s^-1]", size=30)
            cbar.ax.tick_params(labelsize=15) 
            plt.draw()
            nc = None
            Q= None
            if bnd is not None: 
                plotpatch(bnd)

            if lim is not None:
                ax.set_xlim([lim[0], lim[1]])
                ax.set_ylim([lim[2], lim[3]])
            else:
                ax.set_xlim([X.min(), X.max()])
                ax.set_ylim([Y.min(), Y.max()])

            ax.set_facecolor('black')
            # ADD ELEVATION
            rect = [0.1,0.1,0.3,0.2] # l,b,w,h
            ax2 = fig.add_axes(rect)

            ax2.plot(time,elev,color='b', lw=2)
#            zeros = elev*0 #Commented out as variable unused
            tide=ax2.plot([time[0],time[0]],[-1,1], color='k')
            ax2.set_ylim([-1,1])
            ax2.set_ylabel('elevation [m]',fontsize = 30)
            ax2.xaxis.set_major_locator(   DayLocator() )
            ax2.xaxis.set_major_formatter( DateFormatter( '%d ' ) )
            ax2.tick_params(labelsize=15)

    def update_img(k):
            global nc,Q,F,tide,h

            N=(k/nt)
#            K=k-(N*nt) #Commented out as variable unused

            for vars in params:
                fullfile=os.path.join(dirin,'schout_'+str(Istart)+'.nc')
                
                if nc is not None:
                    if nc.filepath() != fullfile :
                        nc.close()
                        nc=netCDF4.Dataset(fullfile)
                        print ('Reading %s' % fullfile)
                else:
                    nc=netCDF4.Dataset(fullfile)
                    print ('Reading %s' % fullfile)
                        
                ##The following section of code prabably needs to be re-written,
                ##out a section that was causing problems with the parameter hvel
                ##and suppressed and error messaage that did not seem to affect the output.
                
#                if 'two'  in nc.variables[vars].dimensions:
#                    tmp=nc.variables[vars][k-1]
#                    u=tmp[...,0]
#                    v=tmp[...,1]
#
#                    if len(v.shape)>1:
#                        u=u[level,:]
#                        v=v[level,:]
#
#                    Z=numpy.sqrt(u**2+v**2)
#                else:
                try:
                    Z=nc.variables[vars][k-1,:] ##This code will occasionally
                except IndexError: #throw index error that appear to have no
                    pass           #effect on the file outputted by the program
                if len(Z.shape)==2:#so I suppressed the error message.
                    Z=Z[:,level]
                elif len(Z.shape)==3:
                    Z=Z[:,level,1]                
            
                for coll in F.collections:
                    coll.remove() 

                Z[Z>Zmax]=Zmax
                Z[Z<Zmin]=Zmin
                levels = numpy.linspace(Zmin, Zmax, 60)
                F=ax.tricontourf(X,Y,ele,Z,cmap=plt.cm.Spectral_r,vmin=Zmin,vmax=Zmax,levels=levels)


                
                if quiver is not None:
                    
                    for filename_quiver,vars_quiver in quiver[0].iteritems():
                        fullfile=os.path.join(dirin,str(Istart+N)+'_'+filename_quiver+'.nc')

                        if nc is not None:
                            if nc.filepath() != fullfile :
                                nc.close()
                                nc=netCDF4.Dataset(fullfile)
                                print ('Reading %s' % fullfile)
                        else:
                            nc=netCDF4.Dataset(fullfile)
                            print ('Reading %s' % fullfile)
                    u,v=vars_quiver.split(' ')
                                    
                    u=nc.variables[u][k,:]
                    v=nc.variables[v][k,:]

                    if len(u.shape)>1:
                        u=u[level,:]
                        v=v[level,:]

                    U = griddata(XY, u,(Xreg,Yreg),method='linear')
                    V = griddata(XY, v,(Xreg,Yreg),method='linear')

#                    mag=numpy.sqrt(U**2+V**2) #Commented out as variable unused
                    
                    U[U==0]=numpy.NaN
                    V[V==0]=numpy.NaN
                    
           
#                    ix=3 #Commented out as variable unused
                   
                    if Q is not None:
                
                       Q.set_UVC(U,V)
                       Q.set_zorder(10)


                    else:
                        Q = ax.quiver(Xreg, Yreg, U, V,scale=quiver_scale,zorder=1,color='k')
#                        qk = ax.quiverkey(Q,0.1,0.1,1,'1m/s') #Commented out as variable unused


                dtime = netCDF4.num2date(nc.variables['time'][k-1],ncs.variables['time'].units)
              
                tide[0].set_data([dtime,dtime],[-1,1])
                tt1.set_text(dtime.strftime("%Y-%m-%d %H:%M"))
                plt.draw()

                plt.savefig( os.path.join( figdir,'current_snapshot_%03i.png'%k), dpi=dpi )
                
               
            return F,Q,tide

     
    
    
        
#    tot= (Iend+1-Istart)*nt #I was unable to determine the use of this code,
#    init_img()              #and it was causing Index Errors, so I commented it out.
#    for kk in range(0,tot):
#        update_img(kk)
    ani = animation.FuncAnimation(fig,update_img,frames=95,interval=fps,init_func=init_img,blit=False)

    writer = animation.FFMpegWriter()
    ani.save(moviefile,writer=writer)#,dpi=dpi)
    figdir=figdir+'/'
    os.system('convert -loop 0 -delay %i %scurrent_snapshot_*.png %s'%(fps,figdir,moviefile) )
    os.system('rm %scurrent_snapshot_*.png'%(figdir) )
    print ('file: %s created'%moviefile)
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='create_movie.py', usage='%(prog)s fileout dirout INDstart INDend params [options]')
    ## main arguments
    parser.add_argument('fileout', type=str,help='name of the output GIF file')
    parser.add_argument('dirout', type=str,help='name of the output where are the SELFE files')
    parser.add_argument('INDstart', type=int,help='First file to take')
    parser.add_argument('INDend', type=int,help='Last file to take')
    parser.add_argument('params', type=str,nargs='+',help='name of the parameter to plot')
    ## options
    parser.add_argument('-lim', type=float,help='Xlim,Ylim',nargs='+')
    parser.add_argument('-zmin', type=float,help='minimum value')
    parser.add_argument('-zmax', type=float,help='maximum value')
    parser.add_argument('-dpi', type=int,help='numper of pixel (default 50)',default=50)
    parser.add_argument('-quiver', type=json.loads,nargs='+',help='name of the quiver variable to plot')
    parser.add_argument('-quiver_res', type=int,help='Quiver resolution (default 10)',default=10)
    parser.add_argument('-quiver_scale', type=int,help='Quiver scale (default 1)',default=1)
    parser.add_argument('-fps', type=int,help='delay',default=60)
    parser.add_argument('-bnd', type=str,help='Blanking file',default=None)
    parser.add_argument('-level',type=int,help='level to plot for 3D data',default=1)
    args = parser.parse_args()



    ### PRINT IT ALL
    print ('output name : %s' % (args.fileout))
    print ('Directory : %s' % (args.dirout))
    print ('From file #%i and #%i' % (args.INDstart,args.INDend))
    print ('Do parameters : %s' % (args.params))



    process(args.fileout,args.dirout,args.INDstart,args.INDend,args.dpi,args.params,args.quiver,args.quiver_res,args.quiver_scale,args.fps,args.zmin,args.zmax,args.bnd,args.level-1,args.lim)
