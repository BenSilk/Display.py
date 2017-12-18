#!/usr/bin/env python2.7

import os,sys
import netCDF4
import numpy
import matplotlib.animation as animation
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.dates import YearLocator, DayLocator, DateFormatter, AutoDateLocator
import matplotlib.pyplot as plt
import json
from scipy.interpolate import griddata
import pylab

def lat_lon_proportion(plt,ax):
	#Calculate the distances along the axis
	xlim=plt.xlim()
	ylim=plt.ylim()
	x_dist = numpy.diff(xlim);
	y_dist = numpy.diff(ylim);

	#Adjust the aspect ratio
	c_adj = numpy.cos(numpy.mean(numpy.deg2rad(xlim)));
	dar = [1,c_adj,1];
	pbar = [x_dist[0]*c_adj/y_dist[0],1,1 ];
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

def get_min_max(nc,vars,level,t):
    
    if len(vars.split(',')) >1:
        	  u,v=vars.split(',')
                  u=nc.variables[u][t,:]
                  v=nc.variables[v][t,:]
                  Z=numpy.sqrt(u**2+v**2)
    else:
          	  Z=nc.variables[vars][t,:]
	    
    if len(Z.shape)>1:
                  Z=Z[level,:]


    Zmin=min(Z)
    Zmax=max(Z)
    return Z,numpy.round(Zmax,2),numpy.round(Zmin,2)

def process(filein,ts,params,quiver,quiver_res,quiver_scale,zmin,zmax,bnd,level,lim):

    fig = plt.figure(figsize=(15,9))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    tt1 = ax.text(.5, 1.05, '', transform = ax.transAxes, va='center',fontsize = 30)
    ax.tick_params(labelsize=25)


    
    # get the static variable
   
    
    ncs=netCDF4.Dataset(filein)
    X=ncs.variables['x'][:]
    Y=ncs.variables['y'][:]
    Z=ncs.variables['depth'][:]
    ele=ncs.variables['ele'][:]-1
    nt=len(ncs.variables['time'])
    
    if quiver is not None:
        if lim is not None:
            Xreg, Yreg = numpy.meshgrid(numpy.linspace(lim[0],lim[1], quiver_res), numpy.linspace(lim[2], lim[3], quiver_res))
        else:
            Xreg, Yreg = numpy.meshgrid(numpy.linspace(min(X[:]),max(X[:]), quiver_res), numpy.linspace(min(Y[:]), max(Y[:]), quiver_res))
        XY=numpy.array((X,Y)).T

    
    if lim is not None:
        node   = near2d_selfe(X[:],Y[:],(lim[1]-lim[0])/2, (lim[3]-lim[2])/2,  nnodes=1, nreturn='single')
    else:
        node   = near2d_selfe(X[:],Y[:],numpy.mean(X[:]),numpy.mean(Y[:]),  nnodes=1, nreturn='single')
    
    
    time = netCDF4.num2date(ncs.variables['time'][ts],ncs.variables['time'].units)

    
    ZZ,ZZmax,ZZmin=get_min_max(ncs,params,level,ts)
    if zmin is None:
        Zmin=ZZmin
    else:
        Zmin=zmin
    if zmax is None:
        Zmax=ZZmax
    else:
        Zmax=zmax
        
          
    ZZ[ZZ>Zmax]=Zmax
    ZZ[ZZ<Zmin]=Zmin
    levels = numpy.linspace(Zmin, Zmax, 60)
    F=plt.tricontourf(X,Y,ele,ZZ,vmin=Zmin,vmax=Zmax,cmap=plt.cm.Spectral_r,levels=levels)
    plt.clim(Zmin,Zmax)
    plt.xlabel('Easting (meters)',fontsize = 30)
    plt.ylabel('Northing (meters)',fontsize = 30)

    cbar=plt.colorbar(F,ax=ax)#numpy.linspace(Zmin,Zmax,10))
    cbar.set_label(r"%s" %(params), size=30)
    cbar.ax.tick_params(labelsize=25) 
    plt.draw()
    if bnd is not None: 
 	plotpatch(bnd)

    if lim is not None:
        ax.set_xlim([lim[0], lim[1]])
        ax.set_ylim([lim[2], lim[3]])
    else:
        ax.set_xlim([X.min(), X.max()])
        ax.set_ylim([Y.min(), Y.max()])

    ax.set_axis_bgcolor('black')
           
    if quiver is not None:
        u,v=quiver.split(',')
                                    
        u=ncs.variables[u][ts,:]
        v=ncs.variables[v][ts,:]

        if len(u.shape)>1:
            u=u[level,:]
            v=v[level,:]

        U = griddata(XY, u,(Xreg,Yreg),method='linear')
        V = griddata(XY, v,(Xreg,Yreg),method='linear')

        mag=numpy.sqrt(U**2+V**2)
                    
        U[U==0]=numpy.NaN
        V[V==0]=numpy.NaN
                    
        Q = ax.quiver(Xreg, Yreg, U, V,scale=quiver_scale,zorder=1,color='k')
        qk = ax.quiverkey(Q,0.1,0.1,1,'1m/s')


       
    tt1.set_text(time.strftime("%Y-%m-%d %H:%M"))
    plt.show()
   
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='display.py', usage='%(prog)s filein ts params [options]')
    ## main arguments
    parser.add_argument('filein', type=str,help='name of the netcdf file')
    parser.add_argument('ts', type=int,help='timestep (native)')
    parser.add_argument('params', type=str,help='name of the parameter to plot')
    ## options
    parser.add_argument('-lim', type=float,help='Xlim,Ylim',nargs='+')
    parser.add_argument('-zmin', type=float,help='minimum value')
    parser.add_argument('-zmax', type=float,help='maximum value')
    parser.add_argument('-quiver', type=str,help='name of the quiver variable to plot')
    parser.add_argument('-quiver_res', type=int,help='Quiver resolution (default 10)',default=10)
    parser.add_argument('-quiver_scale', type=int,help='Quiver scale (default 1)',default=1)
    parser.add_argument('-bnd', type=str,help='Blanking file',default=None)
    parser.add_argument('-level',type=int,help='level to plot for 3D data',default=1)
    args = parser.parse_args()



    ### PRINT IT ALL


    process(args.filein,args.ts,args.params,args.quiver,args.quiver_res,args.quiver_scale,args.zmin,args.zmax,args.bnd,args.level-1,args.lim)
