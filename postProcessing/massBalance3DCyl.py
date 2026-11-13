#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: massBalance3DCyl.py
Date: 2026-11-13
Version: 1.0
License: GPL-3.0-or-later
Copyright 2026 Alban Gilletta, Cyrille Bonamy, Julien Chauchat
This file is part of ExnerAnalysis.

ExnerAnalysis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

ExnerAnalysis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with ExnerAnalysis. If not, see <https://www.gnu.org/licenses/>.

Contact: alban.gilletta-de-saint-joseph@univ-grenoble-alpes.fr
"""

#
# Import section
#

import subprocess
import sys
import numpy as np
import fluidfoam
from tqdm import tqdm
import h5py

from func1 import create_point_2Dcyl, save_point, read_point, uns2cylvec, uns2cyl, ddt

class Threedimsimu(object):
    def __init__(self, basepath, case, gC, gZw, gZ1, gZ1b, nzw, nzbw, nz1, nz1b, read_point_nc, save_point_nc, read_mesh, read_result, asint):
        self.basepath = basepath
        self.case = case
        self.sol = self.basepath + self.case + "/"
        self.solsav = self.basepath + self.case + "/constant/"
        self.asint = asint 					# Solid volume fraction criteria
        self.D = 0.1						# Pile diameter
        self.radius = 2*self.D*np.sqrt(2)	# Sandpit radius
        self.read_point_nc = read_point_nc
        self.save_point_nc = save_point_nc
        
        self.createmesh(gC, gZw, gZ1, gZ1b, nzw, nzbw, nz1, nz1b, read_mesh, read_result)
        
    def createmesh(self, gC, gZw, gZ1, gZ1b, nzw, nzbw, nz1, nz1b, read_mesh, read_result):
        if not read_mesh:
            x, y, z = fluidfoam.readmesh(self.sol, precision = 13)
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'w')
            hf.create_dataset('x',data=x)
            hf.create_dataset('y',data=y)
            hf.create_dataset('z',data=z)
            hf.close()
        else :
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'r')
            x = np.array(hf.get('x'))
            y = np.array(hf.get('y'))
            z = np.array(hf.get('z'))
            
        ncell = np.size(x)
        if not self.read_point_nc:
            n2d, nz, pbed, pointer = create_point_2Dcyl(x, y, z)
            if self.save_point_nc:
                save_point(self.solsav, n2d, nz, pbed, x, y, z, pointer)
        else:
            n2d, nz, pbed, pointer = read_point(self.solsav)
            
        ncf = nz*n2d 			# Cell number in the cylindrical part
        theta = np.zeros(n2d)
        r = np.zeros(n2d)
        x = x[pbed]
        y = y[pbed]
        z = z[pointer][0,:]
        
        ### CONVERT INTO CYLINDRICAL COORDINATES ###
        
        for i in range(n2d):
            if np.logical_and(x[i] > 0, y[i] > 0):
                theta[i]=np.arctan(y[i]/x[i])
            elif np.logical_and(x[i] > 0, y[i] < 0):
                theta[i]=np.arctan(y[i]/x[i]) + 2*np.pi
            elif (x[i] < 0):
                theta[i]=np.arctan(y[i]/x[i]) + np.pi
            elif np.logical_and(x[i] == 0, y[i] > 0):
                theta[i]= np.pi/2
            elif np.logical_and(x[i] == 0, y[i] < 0):
                theta[i]= 3*np.pi/2

            r[i]=np.sqrt(x[i]**2+y[i]**2) # radial distance from cell centers to pile center
        
        R = self.D/2															# Pile radius
        nr = np.unique(np.round(r,6)).shape[0]									# Cells number in radial direction inside sandpit
        dr, drN = fluidfoam.meshdesign.getdzs(self.radius-R, gC, nr)			# Get radial cell size
        rval, deltar, gC = fluidfoam.meshdesign.getgz(self.radius-R, dr, nr)	# Get radial edge coordinates
        ntheta = np.unique(np.round(theta,6)).shape[0]							# Cells number in orthoradial direction inside sandpit
        z = np.unique(np.round(z,6))											# Cells vertical coordiantes inside sandpit
        nz = np.unique(np.round(z,6)).shape[0]									# Cells number in vertical direction inside sandpit
        indcyl = np.zeros((nr, ntheta))											# Indice for cylindrical coordinates
        rcyl=np.unique(np.round(r,6))
        thetacyl = np.unique(np.round(theta,6))
        indcart = range(n2d)													# Indice for cartesian coordinates array flattened
        dtheta = 2*np.pi/ntheta
        rcylint = rcyl - deltar/2
        rcylext = rcyl + deltar/2
        
        
        ### CORRESPONDING INDICE OF CARTESIAN ARRAY indcart IN CYLINDRICAL ARRAY indcyl ###
        
        i=-1
        j=-1
        for d in tqdm(np.unique(np.round(r,6))):
            j=-1
            i=i+1
            for t in np.unique(np.round(theta,6)):
                j=j+1
                indcyl[i,j]=int(np.argwhere(np.logical_and(np.round(r,6)==d,np.round(theta,6)==t))[0][0])
        indcyl=indcyl.astype(int)

        print(np.shape(r), np.shape(theta))
        r = np.unique(np.round(r,6))
        theta = np.unique(np.round(theta,6))
        
        #
        # Construct vertical mesh M0123456
        #

        Lzw  = self.D 		# Length upper part (water)
        Lzbw = self.D 		# Length lower part (water)
        Lz1  = self.D/2 	# Length upper part (sediment)
        Lz1b = 0.04 		# Length lower part (sediment)
        
        dzw = np.ones(nz)*Lzw/nzw
        ziw = np.linspace(self.D,2*self.D,nzw+1)
        dzbw, dzbNw = fluidfoam.meshdesign.getdzs(Lzbw, gZw, nzbw)
        zbiw, dzbw, gZw = fluidfoam.meshdesign.getgz(Lzbw, dzbw, nzbw)
        dz1, dz1N = fluidfoam.meshdesign.getdzs(Lz1, gZ1, nz1)
        z1i, dz1, gZ1 = fluidfoam.meshdesign.getgz(Lz1, dz1, nz1)
        dz1b, dz1bN = fluidfoam.meshdesign.getdzs(Lz1b, gZ1b, nz1b)
        z1bi, dz1b, gZ1b = fluidfoam.meshdesign.getgz(Lz1b, dz1b, nz1b)
        gZ1b = 1/gZ1b

        deltaZ = np.concatenate((dz1b, dz1, dzbw, dzw))
        zp = np.concatenate((z1bi-(Lz1+Lz1b), z1i[1:]-Lz1, zbiw[1:],ziw[1:]))
        
        self.x = x
        self.y = y
        self.z = z
        self.ncell = ncell
        self.n2d = n2d
        self.nz = nz
        self.pbed = pbed
        self.pointer = pointer
        self.ncf = ncf
        self.theta = theta
        self.r = r
        self.R = R
        self.nr = nr
        self.dr = dr
        self.deltar = deltar
        self.ntheta = ntheta
        self.indcyl = indcyl
        self.indcart = indcart
        self.dtheta = dtheta
        self.deltaZ = deltaZ
        self.zp = zp
    
        self.createtime(read_result)
        
    def createtime(self, read_result):
        try:
            proc = subprocess.Popen(["foamListTimes", "-case", self.sol], stdout=subprocess.PIPE)
        except FileNotFoundError:
            print("foamListTimes : command not found")
            print("Do you have load OpenFoam environement?")
            sys.exit(0)
        output = proc.stdout.read()
        tread = output.decode().rstrip().split("\n")
        if len(tread) != 1 :
            tread= tread[:480]
            dt = float(tread[1])-float(tread[0])
            nt = len(tread)
        else :
            tread = np.arange(0.25,120.25,0.25)
            tread = np.array(tread).astype(str)
            dt = 0.25
            nt = 480
        k=-1
        for s in tread:
            k=k+1
            s=s.replace(".0","")
            tread[k]=s
            
        time = np.zeros(nt)
        
        self.time = time
        self.dt = dt
        self.nt = nt
        self.tread = tread
        
        self.postProcess(read_result)
          
    def postProcess(self, read_result):
        k = -1
        if not read_result:
            dVuns = np.zeros(self.ncell)
            dV = np.zeros((self.nr, self.ntheta, self.nz))
            print(np.shape(self.r), np.shape(self.theta), np.shape(dV))
            alpha = np.zeros((self.ncell, self.nt))
            alphat, dalphat_sur_dt = np.zeros((2, self.nr, self.ntheta, self.nz, self.nt))
            aU = np.zeros((3, self.ncell, self.nt))
            aUa = np.zeros((3, self.nr, self.ntheta, self.nz, self.nt))
    
            zbed, dzbed_sur_dt, dVolPhi_sur_dt, VoldPhi_sur_dt, VoldPhisusp_sur_dt, VoldPhibed_sur_dt = np.zeros((6, self.nr, self.ntheta, self.nt))
            phi_dzbed_sur_dt, intdVolPhi_sur_dt, intVolPhi = np.zeros((3, self.nt))
            intdVolPhibed_sur_dt, intVoldPhibed_sur_dt, intVolPhibed, intdVolPhisusp_sur_dt, intVoldPhisusp_sur_dt, intVolPhisusp = np.zeros((6, self.nt))
            intS_PindS, intS2_PindS, intS3_PindS, intSlat_PindS, intSout_PindS, intSlat1st_PindS, intSlat1end_PindS, intStop_PindS = np.zeros((8, self.nt))
            intSlatbed_PindS, intSoutbed_PindS, intSlat1stbed_PindS, intSlat1endbed_PindS, intSlatsusp_PindS, intSoutsusp_PindS, intSlat1stsusp_PindS, intSlat1endsusp_PindS = np.zeros((8, self.nt))
    
            dVuns[:] = fluidfoam.readscalar(self.sol, "constant", "V", verbose=False, precision=13)
            varorig = dVuns[:]
            varorig = varorig[self.pointer]
            dV[:,:,:] = uns2cyl(varorig, self.nr, self.ntheta, self.nz, self.indcyl, self.indcart)
            print(np.sum(varorig), np.sum(dV))
            
            k = -1
            for t in tqdm(self.tread):
                print("Reading time: %s s" % t)
                k = k + 1
                alpha[:, k] = fluidfoam.readscalar(self.sol, t, "alpha.a", verbose=False, precision=13)
                aU[:,:,k] = fluidfoam.readvector(self.sol, t, "alphaUa", verbose=False, precision=13)
                self.time[k] = float(t)
                varorig = alpha[:,k]
                varorig = varorig[self.pointer]
                alphat[:,:,:,k] = uns2cyl(varorig, self.nr, self.ntheta, self.nz, self.indcyl, self.indcart)
                varorig = aU[:,:,k]
                varorig = varorig[:,self.pointer]
                print(np.shape(varorig))
                aUa[:,:,:,:,k] = uns2cylvec(varorig, self.nr, self.ntheta, self.nz, self.indcyl, self.indcart)
                        
            thetaZ1st= int(self.ntheta*1/8)  				# Zone 1 : int(self.ntheta*3/8), Zone 2l : int(self.ntheta*1/8), Zone 2r : int(self.ntheta*5/8), Zone 3: int(self.ntheta*7/8) 
            thetaZ1end= int(self.ntheta*3/8) 				# Zone 1 : int(self.ntheta*5/8), Zone 2l : int(self.ntheta*3/8), Zone 2r : int(self.ntheta*7/8), Zone 3: int(self.ntheta*1/8)
            rangetheta = range(thetaZ1st, thetaZ1end)   	# Zone 1, 2l, 2r : range(thetaZ1st, thetaZ1end) Zone 3 : np.concatenate((range(thetaZ1st, self.ntheta), range(thetaZ1end)))
            rangethetam1 = range(thetaZ1st, thetaZ1end-1)	# Zone 1, 2l, 2r : range(thetaZ1st, thetaZ1end-1) Zone 3 : chain(range(int(ntheta*7/8),ntheta), range(int(ntheta*0.125)))
            Nx = self.nr                       
    
            for ri in range(Nx-1):
                for thi in rangetheta:
                    for t in range(self.nt):
                        bedcondi = np.where(alphat[ri,thi,:,t] <= self.asint)
                        suspcondi = np.where(alphat[ri,thi,:,t] <= 0.08)
                        zbed[ri, thi, t] = self.zp[bedcondi[0][0]]
                        intVolPhi[t] = intVolPhi[t] + np.sum((alphat[ri,thi,bedcondi,t])*dV[ri,thi,bedcondi])
                        intVolPhisusp[t] = intVolPhisusp[t] + np.sum((alphat[ri,thi,suspcondi,t])*dV[ri,thi,suspcondi])
                        intVolPhibed[t] = intVolPhibed[t] + np.sum((alphat[ri,thi,bedcondi[0][0]:suspcondi[0][0],t])*dV[ri,thi,bedcondi[0][0]:suspcondi[0][0]])
                    dzbed_sur_dt[ri, thi, :]=ddt('eulerFor', zbed[ri, thi, :], self.dt)
                    
            intdVolPhi_sur_dt = ddt('eulerFor', intVolPhi, self.dt)
            intdVolPhibed_sur_dt = ddt('eulerFor', intVolPhibed, self.dt)
            intdVolPhisusp_sur_dt = ddt('eulerFor', intVolPhisusp, self.dt)
        
            for ri in range(Nx-1):
                for thi in rangetheta:
                    for k in range(self.nz):
                        dalphat_sur_dt[ri, thi, k, :]=ddt('eulerFor', alphat[ri, thi, k, :], self.dt)
    
            for ri in range(Nx-1):
                for thi in rangetheta:
                    for t in range(self.nt):
                        bedcondi = np.where(alphat[ri,thi,:,t] <= self.asint)
                        bedcondi1 = np.where(alphat[ri+1,thi,:,t] <= self.asint)
                        suspcondi = np.where(alphat[ri,thi,:,t] <= 0.08)
                        if bedcondi[0][0] == 0 :
                            phi_dzbed_sur_dt[t] = phi_dzbed_sur_dt[t] + alphat[ri, thi, bedcondi[0][0], t]*dzbed_sur_dt[ri, thi, t]*(self.dtheta*self.r[ri]*self.deltar[ri])
                        else :
                            phi_dzbed_sur_dt[t] = phi_dzbed_sur_dt[t] + (alphat[ri, thi, bedcondi[0][0], t]*self.deltaZ[bedcondi[0][0]-1]+alphat[ri, thi, bedcondi[0][0]-1, t]*self.deltaZ[bedcondi[0][0]])/(self.deltaZ[bedcondi[0][0]]+self.deltaZ[bedcondi[0][0]-1])*(dzbed_sur_dt[ri, thi, t])*(self.dtheta*self.r[ri]*self.deltar[ri])
                            intS_PindS[t] = intS_PindS[t] - (aUa[2, ri, thi, bedcondi[0][0], t]*self.deltaZ[bedcondi[0][0]-1]+aUa[2, ri, thi, bedcondi[0][0]-1, t]*self.deltaZ[bedcondi[0][0]])/(self.deltaZ[bedcondi[0][0]]+self.deltaZ[bedcondi[0][0]-1])*(self.dtheta*self.r[ri]*self.deltar[ri])
                        VoldPhi_sur_dt[ri, thi, t] = np.sum(dalphat_sur_dt[ri,thi,bedcondi,t]*dV[ri,thi,bedcondi])
                        VoldPhisusp_sur_dt[ri, thi, t] = np.sum(dalphat_sur_dt[ri,thi,suspcondi,t]*dV[ri,thi,suspcondi])
                        VoldPhibed_sur_dt[ri, thi, t] = np.sum(dalphat_sur_dt[ri,thi,bedcondi[0][0]:suspcondi[0][0],t]*dV[ri,thi,bedcondi[0][0]:suspcondi[0][0]])
                        intdivPi[t] = intdivPi[t] + np.sum(divPi[ri,thi,bedcondi[0][0]:,t]*dV[ri,thi,bedcondi[0][0]:])
                        intdivPisusp[t] = intdivPisusp[t] + np.sum(divPi[ri,thi,suspcondi[0][0]:,t]*dV[ri,thi,suspcondi[0][0]:])
                        intdivPibed[t] = intdivPibed[t] + np.sum(divPi[ri,thi,bedcondi[0][0]:suspcondi[0][0],t]*dV[ri,thi,bedcondi[0][0]:suspcondi[0][0]])
                        intStop_PindS[t] = intStop_PindS[t] +  aUa[2, ri, thi, -1, t]*(self.dtheta*self.r[ri]*self.deltar[ri])
                        if ri!=Nx-2:
                            deltaZx = zbed[ri+1, thi, t]-zbed[ri, thi, t]
                            match deltaZx:
                                case deltaZx if deltaZx > 0:
                                    S2 = np.where(np.logical_and(alphat[ri,thi,:,t] <= self.asint, alphat[ri+1,thi,:,t] > self.asint))[0]
                                    intS2_PindS[t] = intS2_PindS[t] + np.sum(((aUa[0, ri, thi, S2, t]*self.deltar[ri+1]+aUa[0, ri+1, thi, S2, t]*self.deltar[ri])*np.cos(self.theta[thi])+(aUa[1, ri, thi, S2, t]*self.deltar[ri+1]+aUa[1, ri+1, thi, S2, t]*self.deltar[ri])*np.sin(self.theta[thi]))/(self.deltar[ri]+self.deltar[ri+1])*self.deltaZ[S2]*self.dtheta*(self.r[ri]+self.deltar[ri]/2))
                                case deltaZx if deltaZx < 0:
                                    S2 = np.where(np.logical_and(alphat[ri,thi,:,t] > self.asint, alphat[ri+1,thi,:,t] <= self.asint))[0]
                                    intS2_PindS[t] = intS2_PindS[t] - np.sum(((aUa[0, ri, thi, S2, t]*self.deltar[ri+1]+aUa[0, ri+1, thi, S2, t]*self.deltar[ri])*np.cos(self.theta[thi])+(aUa[1, ri, thi, S2, t]*self.deltar[ri+1]+aUa[1, ri+1, thi, S2, t]*self.deltar[ri])*np.sin(self.theta[thi]))/(self.deltar[ri]+self.deltar[ri+1])*self.deltaZ[S2]*self.dtheta*(self.r[ri]+self.deltar[ri]/2))
                       
            for ri in range(Nx-1):
                for thi in rangetheta:
                    for t in range(self.nt):
                        thip1 = (thi+1)%self.ntheta
                        bedcondi = np.where(alphat[ri,thi,:,t] <= self.asint)
                        bedcondi1 = np.where(alphat[ri,thip1,:,t] <= self.asint)
                        deltaZt = zbed[ri, thip1, t]-zbed[ri, thi, t]
                        match deltaZt:
                            case deltaZt if deltaZt > 0:
                                S3 = np.where(np.logical_and(alphat[ri,thi,:,t] <= self.asint, alphat[ri,thip1,:,t] > self.asint))[0]
                                intS3_PindS[t] = intS3_PindS[t] + np.sum((-(aUa[0, ri, thi, S3, t]+aUa[0, ri+1, thip1, S3, t])*np.sin(self.theta[thi]+self.dtheta*0.5)+(aUa[1, ri, thi, S3, t]+aUa[1, ri+1, thip1, S3, t])*np.cos(self.theta[thi]+self.dtheta*0.5))*0.5*self.deltaZ[S3]*self.deltar[ri])
                            case deltaZt if deltaZt < 0:
                                S3 = np.where(np.logical_and(alphat[ri,thi,:,t] > self.asint, alphat[ri,thip1,:,t] <= self.asint))[0]
                                intS3_PindS[t] = intS3_PindS[t] - np.sum((-(aUa[0, ri, thi, S3, t]+aUa[0, ri+1, thip1, S3, t])*np.sin(self.theta[thi]+self.dtheta*0.5)+(aUa[1, ri, thi, S3, t]+aUa[1, ri+1, thip1, S3, t])*np.cos(self.theta[thi]+self.dtheta*0.5))*0.5*self.deltaZ[S3]*self.deltar[ri])
                                    
            for thi in rangetheta:
                for t in range(self.nt):
                    bedcondi = np.where(alphat[Nx-2,thi,:,t] <= self.asint)
                    suspcondi = np.where(alphat[Nx-2,thi,:,t] <= 0.08)
                    for k in range(bedcondi[0][0], self.nz):
                        aUaX = (aUa[0, Nx-2, thi, k, t]*self.deltar[Nx-1]+aUa[0, Nx-1, thi, k, t]*self.deltar[Nx-2])/(self.deltar[Nx-2]+self.deltar[Nx-1])
                        aUaY = (aUa[1, Nx-2, thi, k, t]*self.deltar[Nx-1]+aUa[1, Nx-1, thi, k, t]*self.deltar[Nx-2])/(self.deltar[Nx-2]+self.deltar[Nx-1])
                        if k >= suspcondi[0][0]:
                            intSoutsusp_PindS[t] = intSoutsusp_PindS[t] + (aUaX*np.cos(self.theta[thi])+aUaY*np.sin(self.theta[thi]))*self.deltaZ[k]*self.dtheta*(self.r[Nx-2]+self.deltar[Nx-2]/2)
                        if k < suspcondi[0][0]:
                            intSoutbed_PindS[t] = intSoutbed_PindS[t] + (aUaX*np.cos(self.theta[thi])+aUaY*np.sin(self.theta[thi]))*self.deltaZ[k]*self.dtheta*(self.r[Nx-2]+self.deltar[Nx-2]/2)
                        intSout_PindS[t] = intSout_PindS[t] + (aUaX*np.cos(self.theta[thi])+aUaY*np.sin(self.theta[thi]))*self.deltaZ[k]*self.dtheta*(self.r[Nx-2]+self.deltar[Nx-2]/2)
    
            for thi in [thetaZ1st, thetaZ1end-1]:
                for t in range(self.nt):
                    for ri in range(Nx-1):
                        bedcondi = np.where(alphat[ri,thi,:,t] <= self.asint)
                        suspcondi = np.where(alphat[Nx-2,thi,:,t] <= 0.08)
                        for k in range(bedcondi[0][0], self.nz):
                            if thi == thetaZ1st:
                                aUaX = 0.5*(aUa[0, ri, thi, k, t]+aUa[0, ri, thi-1, k, t])
                                aUaY = 0.5*(aUa[1, ri, thi, k, t]+aUa[1, ri, thi-1, k, t])
                                if k >= suspcondi[0][0]:
                                    intSlat1stsusp_PindS[t] = intSlat1stsusp_PindS[t] + (aUaX*np.sin(self.theta[thi]-self.dtheta/2)-aUaY*np.cos(self.theta[thi]-self.dtheta/2))*self.deltaZ[k]*self.deltar[ri]
                                if k < suspcondi[0][0]:
                                    intSlat1stbed_PindS[t] = intSlat1stbed_PindS[t] + (aUaX*np.sin(self.theta[thi]-self.dtheta/2)-aUaY*np.cos(self.theta[thi]-self.dtheta/2))*self.deltaZ[k]*self.deltar[ri]
                                intSlat1st_PindS[t] = intSlat1st_PindS[t] + (aUaX*np.sin(self.theta[thi]-self.dtheta/2)-aUaY*np.cos(self.theta[thi]-self.dtheta/2))*self.deltaZ[k]*self.deltar[ri]
                            if thi == thetaZ1end-1:
                                aUaX = 0.5*(aUa[0, ri, thi, k, t]+aUa[0, ri, (thi+1)%self.ntheta, k, t])
                                aUaY = 0.5*(aUa[1, ri, thi, k, t]+aUa[1, ri, (thi+1)%self.ntheta, k, t])
                                if k >= suspcondi[0][0]:
                                    intSlat1endsusp_PindS[t] = intSlat1endsusp_PindS[t] + (-aUaX*np.sin(self.theta[thi]+self.dtheta/2)+aUaY*np.cos(self.theta[thi]+self.dtheta/2))*self.deltaZ[k]*self.deltar[ri]
                                if k < suspcondi[0][0]:
                                    intSlat1endbed_PindS[t] = intSlat1endbed_PindS[t] + (-aUaX*np.sin(self.theta[thi]+self.dtheta/2)+aUaY*np.cos(self.theta[thi]+self.dtheta/2))*self.deltaZ[k]*self.deltar[ri]
                                intSlat1end_PindS[t] = intSlat1end_PindS[t] + (-aUaX*np.sin(self.theta[thi]+self.dtheta/2)+aUaY*np.cos(self.theta[thi]+self.dtheta/2))*self.deltaZ[k]*self.deltar[ri]
    
            for t in range(self.nt):
                intSlat_PindS[t] = np.sum([intSout_PindS[t],intSlat1st_PindS[t],intSlat1end_PindS[t]])
                intSlatbed_PindS[t] = np.sum([intSoutbed_PindS[t],intSlat1stbed_PindS[t],intSlat1endbed_PindS[t]])
                intSlatsusp_PindS[t] = np.sum([intSoutsusp_PindS[t],intSlat1stsusp_PindS[t],intSlat1endsusp_PindS[t]])
            
            hf = h5py.File(self.solsav+'/read_massBalanceZ1'+str(self.asint).replace('.','p')+'_t'+self.tread[-1].replace('.','p')+'.h5', 'w')
            hf.create_dataset('phi_dzbed_sur_dt',data=phi_dzbed_sur_dt)
            hf.create_dataset('intdVolPhi_sur_dt',data=intdVolPhi_sur_dt)
            hf.create_dataset('intdVolPhibed_sur_dt',data=intdVolPhibed_sur_dt)
            hf.create_dataset('intdVolPhisusp_sur_dt',data=intdVolPhisusp_sur_dt)
            hf.create_dataset('intS_PindS',data=intS_PindS)
            hf.create_dataset('intS2_PindS',data=intS2_PindS)
            hf.create_dataset('intS3_PindS',data=intS3_PindS)
            hf.create_dataset('intSlat_PindS',data=intSlat_PindS)
            hf.create_dataset('intSlat1st_PindS',data=intSlat1st_PindS)
            hf.create_dataset('intSlat1end_PindS',data=intSlat1end_PindS)
            hf.create_dataset('intSout_PindS',data=intSout_PindS)
            hf.create_dataset('intSlatbed_PindS',data=intSlatbed_PindS)
            hf.create_dataset('intSlat1stbed_PindS',data=intSlat1stbed_PindS)
            hf.create_dataset('intSlat1endbed_PindS',data=intSlat1endbed_PindS)
            hf.create_dataset('intSoutbed_PindS',data=intSoutbed_PindS)
            hf.create_dataset('intSlatsusp_PindS',data=intSlatsusp_PindS)
            hf.create_dataset('intSlat1stsusp_PindS',data=intSlat1stsusp_PindS)
            hf.create_dataset('intSlat1endsusp_PindS',data=intSlat1endsusp_PindS)
            hf.create_dataset('intSoutsusp_PindS',data=intSoutsusp_PindS)
            hf.create_dataset('intStop_PindS',data=intStop_PindS)
            hf.close()
        else :
            for t in tqdm(self.tread):
                k = k + 1
                self.time[k] = float(t)
            hf = h5py.File(self.solsav+'/read_massBalanceZ1'+str(self.asint).replace('.','p')+'.h5', 'r')
            phi_dzbed_sur_dt = np.array(hf.get('phi_dzbed_sur_dt_final'))
            intdVolPhi_sur_dt = np.array(hf.get('intdVolPhi_sur_dt_final'))
            intdVolPhibed_sur_dt = np.array(hf.get('intdVolPhibed_sur_dt_final'))
            intdVolPhisusp_sur_dt = np.array(hf.get('intdVolPhisusp_sur_dt_final'))
            intS_PindS = np.array(hf.get('intS_PindS_final'))
            intS2_PindS = np.array(hf.get('intS2_PindS_final'))
            intS3_PindS = np.array(hf.get('intS3_PindS_final'))
            intSlat_PindS = np.array(hf.get('intSlat_PindS_final'))
            intSlat1st_PindS = np.array(hf.get('intSlat1st_PindS_final'))
            intSlat1end_PindS = np.array(hf.get('intSlat1end_PindS_final'))
            intSout_PindS = np.array(hf.get('intSout_PindS_final'))
            intSlatbed_PindS = np.array(hf.get('intSlatbed_PindS_final'))
            intSlat1stbed_PindS = np.array(hf.get('intSlat1stbed_PindS_final'))
            intSlat1endbed_PindS = np.array(hf.get('intSlat1endbed_PindS_final'))
            intSoutbed_PindS = np.array(hf.get('intSoutbed_PindS_final'))
            intSlatsusp_PindS = np.array(hf.get('intSlatsusp_PindS_final'))
            intSlat1stsusp_PindS = np.array(hf.get('intSlat1stsusp_PindS_final'))
            intSlat1endsusp_PindS = np.array(hf.get('intSlat1endsusp_PindS_final'))
            intSoutsusp_PindS = np.array(hf.get('intSoutsusp_PindS_final'))
            intStop_PindS = np.array(hf.get('intStop_PindS_final'))
        
        # Forward [:-1] Backward [1:] for temporal derivative
        self.phi_dzbed_sur_dt = phi_dzbed_sur_dt[:-1]				# term 1 : bed evolution
        self.intdVolPhi_sur_dt = intdVolPhi_sur_dt[:-1]				# term 2 : storage evolution
        self.intdVolPhibed_sur_dt = intdVolPhibed_sur_dt[:-1]		# term 2a : storage evolution in bedload layer
        self.intdVolPhisusp_sur_dt = intdVolPhisusp_sur_dt[:-1]		# term 2b : storage evolution in suspended layer
        self.intSlat_PindS = intSlat_PindS[:-1]						# term 3 : lateral flux
        self.intSlat1st_PindS = intSlat1st_PindS[:-1]				# term 3a : lateral flux -> through surface defined by theta = thetaZ1st
        self.intSlat1end_PindS = intSlat1end_PindS[:-1]				# term 3b : lateral flux -> through surface defined by theta = thetaZ1end - 1
        self.intSout_PindS = intSout_PindS[:-1]						# term 3c : lateral flux -> through surface defined by r = Nx - 1
        self.intSlatbed_PindS = intSlatbed_PindS[:-1]				# term 3bed : lateral flux -> bedload layer
        self.intSlat1stbed_PindS = intSlat1stbed_PindS[:-1]			# term 3abed : lateral flux -> through surface defined by theta = thetaZ1st in bedload layer
        self.intSlat1endbed_PindS = intSlat1endbed_PindS[:-1]		# term 3bbed : lateral flux -> through surface defined by theta = thetaZ1end - 1 in bedload layer
        self.intSoutbed_PindS = intSoutbed_PindS[:-1]				# term 3cbed : lateral flux -> through surface defined by r = Nx - 1 in bedload layer
        self.intSlatsusp_PindS = intSlatsusp_PindS[:-1]				# term 3susp : lateral flux -> suspended layer
        self.intSlat1stsusp_PindS = intSlat1stsusp_PindS[:-1]		# term 3asusp : lateral flux -> through surface defined by theta = thetaZ1st in suspended layer
        self.intSlat1endsusp_PindS = intSlat1endsusp_PindS[:-1]		# term 3bsusp : lateral flux -> through surface defined by theta = thetaZ1end - 1 in suspended layer
        self.intSoutsusp_PindS = intSoutsusp_PindS[:-1]        		# term 3csusp : lateral flux -> through surface defined by r = Nx - 1 in suspended layer
        self.intS_PindS = intS_PindS[:-1]							# term 4a : bed flux -> vertical component
        self.intS2_PindS = intS2_PindS[:-1]							# term 4b : bed flux -> radial component
        self.intS3_PindS = intS3_PindS[:-1]							# term 4c : bed flux -> orthoradial component
        self.intStop_PindS = intStop_PindS[:-1]						# term 5 : top flux

