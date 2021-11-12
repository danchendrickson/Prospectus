import numpy as np

class EFIT:
    __ts = 0
    __ds = 0

    def __init__(self, xGrid, yGrid, zGrid, tStep, dStep):
        #Initialize with the size of the grid in X, Y, and Z number of nodes.  The distance step, and the time step


        #Velocity grid with 3 part vector at each node point, for 2 time steps
        self.GridShapeV = (3,xGrid+1,yGrid+1,zGrid+1)
        #Stresses with 3 part vector in each of 3 part phaces, at each nod epoint, for 2 time steps
        self.GridShapeS = (3,3,xGrid,yGrid,zGrid)
        #materials property gird.  Initially 3 properties needed: density, Lame 1, Lame 2
        self.GridShapeP = (3,xGrid+1,yGrid+1,zGrid+1)
     
        self.GridPointsV = (xGrid+1)*(yGrid+1)*(zGrid+1)
        self.GridPointsS = (xGrid)*(yGrid)*(zGrid)
        
        
        #define empty grid for the 3 directions of velocity for 2 time steps
        self.Gv = np.zeros(3*self.GridPointsV,dtype="float32").reshape(*self.GridShapeV)
        #define empty grid for the 3 directions of stress on 3 dimmensions of faces for 2 time steps
        self.Gs = np.zeros(3*3*self.GridPointsS,dtype="float32").reshape(*self.GridShapeS)
        #define empty grid for the 3 scalar material properties at each node point.  Can honly hold scalar properties
        #Assumed properties are density, Lame 1, Lame 2
        self.Gp = np.zeros(3*self.GridPointsS,dtype="float32").reshape(*self.GridShapeP)
        
        self.MaxX = xGrid
        self.MaxY = yGrid
        self.MaxZ = zGrid

        self.ds = dStep
        self.ts = tStep

    def CheckStressBoundary(self,x,y,z,Ds):
        #checks to see if a grid is at a boundary, and if so, adjusts boundary conditions appropriately
        #
        # Inputs: x,y,z coordinates of cube in question
        #
        # Outputs: Updated (if boundary) delta stress matrix

        #at front and back faces, stresses perpendicular to the face are 0:
        if x == 1 or x == self.MaxX-1:
            Ds[0,0]=0
            #Ds[0,1]=0
            #Ds[0,2]=0
        
        #at top face, stresses perpendicular to the face are 0:
        if y == self.MaxY-1:
            #Ds[1,0]=0
            Ds[1,1]=0
            #Ds[1,2]=0
        
        if y == 1:
            Ds[1,1]=-self.Gs[1,1,x,y+1,z]
             
        
        #at side faces, stresses perpendicular to the face are 0:
        if z == 1 or z == self.MaxZ-1:
            #Ds[2,0]=0
            #Ds[2,1]=0
            Ds[2,2]=0
        
        return Ds

    def CheckVelocityBoundary(self,x,y,z,Dv):
        #checks to see if a grid is at a boundary, and if so, adjusts boundary conditions appropriately
        #
        # Inputs: x,y,z coordinates of cube in question
        #
        # Outputs: Updated (if boundary) delta velocity vector
        
        #lower boundary is fixed, velocity is Zero
        if y == 0:
            Dv[1]=0
            if x <= 10 or x >= self.MaxX-1:
                Dv[0] = 0
                Dv[2] = 0

        if x == 0 or y == 0 or z == 0 or x == self.MaxX+1 or y == self.MaxY+1 or z == self.MaxZ+1:
            Dv[:]=0

        return Dv
    
    def DeltaStress(self,x,y,z):
        #Gets the change in the stresses per time at a certain coordinate juncture 
        #
        # Inputs: x,y,z coordinates of the cube in question.  Last time is assumed
        #
        # Outputs: 6 dimmensions of stress.  
        
        #Calculated stresses based on 3.55
        Ds = np.zeros((3,3))
        
        Lame1=self.Gp[1,x,y,z]
        Lame2=self.Gp[2,x,y,z]
        try:
            if x == 0:
                Ds[0,0] = 0
            elif y==0 and z == 0:
                Ds[0,0] = ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]) 
                     )
                    )
            elif y==0:
                Ds[0,0] = ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]) +
                        Lame1*(self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1])
                        )
                    )
            elif z == 0:
                Ds[0,0] = ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]) +
                        Lame1*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z])
                        )
                    )
            else:
                Ds[0,0] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]) +
                        Lame1*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z]+self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1])
                        )
                    )
        except:
            #try:
            #    Ds = self.CheckStressBoundary(x,y,z,Ds)
            #except:
                print('opps Ds00',x,y,z)
                Ds[0,0]=0
        
        try:
            if y == 0:
                Ds[1,1] = 0
            elif x==0 and z == 0:
                Ds[1,1] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z])
                     )
                    )                
            elif x == 0:
                Ds[1,1] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z]) +
                        Lame1*(self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1])
                        )
                    )
            elif z == 0:
                Ds[1,1] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z]) +
                        Lame1*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z])
                        )
                    )
            else:                
                Ds[1,1] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z]) +
                        Lame1*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]+self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1])
                        )
                    )
        except:
            #try:
            #    Ds = self.CheckStressBoundary(x,y,z,Ds)
            #except:
                print('opps Ds11',x,y,z)
                Ds[1,1]=0        

        try:
            if z == 0:
                Ds[2,2]=0
            elif x == 0 and y == 0:
                Ds[2,2] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1])
                     )
                    )
            elif x == 0:
                Ds[2,2] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1]) +
                        Lame1*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z])
                        )
                    )
            elif y == 0:
                Ds[2,2] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1]) +
                        Lame1*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z])
                        )
                    )
            else:
                Ds[2,2] =  ((1/self.ds) *
                    ((Lame1+2*Lame2)*(self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1]) +
                        Lame1*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]+self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z])
                        )
                    )
    
        except:
            #try:
            #    Ds = self.CheckStressBoundary(x,y,z,Ds)
            #except:
                print('opps Ds22',x,y,z)
                Ds[2,2]=0        

        try:
            if x == self.MaxX-1 and y == self.MaxY-1:
                Ds[0,1] =  0
            elif x == self.MaxX-1:
                Ds[0,1] =  (
                    (1/self.ds) *
                    (2/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y+1,z]))) *
                    (self.Gv[0,x,y+1,z]-self.Gv[0,x,y,z] )
                    )
            elif y == self.MaxY-1:
                Ds[0,1] =  (
                    (1/self.ds) *
                    (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y-1,z])+(1/self.Gp[2,x+1,y-1,z]))) *
                    (self.Gv[1,x+1,y,z]-self.Gv[1,x,y,z] )
                    )
            else:
                Ds[0,1] =  (
                    (1/self.ds) *
                    (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x+1,y+1,z]))) *
                    (self.Gv[0,x,y+1,z]-self.Gv[0,x,y,z] +self.Gv[1,x+1,y,z]-self.Gv[1,x,y,z] )
                    )

        except:
            #try:
            #    Ds = self.CheckStressBoundary(x,y,z,Ds)
            #except:
                print('opps Ds01',x,y,z)
                Ds[0,1]=0  
        Ds[1,0]=Ds[0,1]      

        try:
            if x == self.MaxX-1 and z == self.MaxZ-1:
                Ds[0,2] =  0
            elif x == self.MaxX-1:
                Ds[0,2] =  (
                    (1/self.ds) *
                    (2/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y,z+1]))) *
                    (self.Gv[0,x,y,z+1]-self.Gv[0,x,y,z])
                    )
            elif z == self.MaxZ-1:
                Ds[0,2] =  (
                    (1/self.ds) *
                    (2/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z]))) *
                    (self.Gv[2,x+1,y,z]-self.Gv[2,x,y,z] )
                    )
            else:
                Ds[0,2] =  (
                    (1/self.ds) *
                    (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x+1,y,z+1]))) *
                    (self.Gv[0,x,y,z+1]-self.Gv[0,x,y,z] +self.Gv[2,x+1,y,z]-self.Gv[2,x,y,z] )
                    )

        except:
            #try:
            #    Ds = self.CheckStressBoundary(x,y,z,Ds)
            #except:
                print('opps Ds02',x,y,z)
                Ds[0,2]=0        
        Ds[2,0]=Ds[0,2]
        
        try:
            if y == self.MaxY-1  and z == self.MaxZ-1:
                Ds[1,2] =  0
            elif y == self.MaxY-1:
                Ds[1,2] =  (
                    (1/self.ds) *
                    (2/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y,z+1]))) *
                    (self.Gv[1,x,y,z+1]-self.Gv[1,x,y,z])
                    )
            elif z == self.MaxZ-1:
                Ds[1,2] =  (
                    (1/self.ds) *
                    (2/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y+1,z]))) *
                    (self.Gv[2,x,y+1,z]-self.Gv[2,x,y,z] )
                    )
            else:
                Ds[1,2] =  (
                    (1/self.ds) *
                    (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x,y+1,z+1]))) *
                    (self.Gv[1,x,y,z+1]-self.Gv[1,x,y,z] +self.Gv[2,x,y+1,z]-self.Gv[2,x,y,z] )
                    )

        except:
            #try:
            #    Ds = self.CheckStressBoundary(x,y,z,Ds)
            #except:
                print('opps Ds12',x,y,z)
                Ds[1,2]=0        
        Ds[2,1] = Ds[1,2]

        Ds = self.CheckStressBoundary(x,y,z,Ds)

        return Ds

    def DeltaVelocity(self, x,y,z):
        #Gets the change in the velocity per time at a certain coordinate juncture 
        #
        # Inputs: x,y,z coordinates of the cube in question.  Last time is assumed
        #
        # Outputs: 3 dimmensions of velocity.  
        
        #Calculated velocity based on 3.54

        DV = np.zeros(3)

        try:
            if x == self.MaxX and y == 0 and z==0:
                DV[0] = 0
            elif y == 0 and z==0:
                DV[0] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x+1,y,z])) *
                        (self.Gs[0,0,x+1,y,z] - self.Gs[0,0,x,y,z])
                        )
            elif x == self.MaxX and z==0:
                DV[0] = ((1 / self.ds ) *
                        (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x,y-1,z])
                        )
            elif x == self.MaxX and y == 0:
                DV[0] = ((1 / self.ds ) *
                        (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x,y+1,z])
                        )
            elif z==0:
                DV[0] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x+1,y,z])) *
                        (self.Gs[0,0,x+1,y,z] - self.Gs[0,0,x,y,z] + self.Gs[0,1,x,y,z] - self.Gs[0,1,x,y-1,z])
                        )
            elif x == self.MaxX:
                DV[0] = ((1 / self.ds ) *
                        (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x,y-1,z] + self.Gs[0,2,x,y,z] -self.Gs[0,2,x,y,z-1])
                        )
            elif y == 0:
                DV[0] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x+1,y,z])) *
                        (self.Gs[0,0,x+1,y,z] - self.Gs[0,0,x,y,z] + self.Gs[0,2,x,y,z] -self.Gs[0,2,x,y,z-1])
                        )
            else:
                DV[0] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x+1,y,z])) *
                        (self.Gs[0,0,x+1,y,z] - self.Gs[0,0,x,y,z] + self.Gs[0,1,x,y,z] - self.Gs[0,1,x,y-1,z] + self.Gs[0,2,x,y,z] -self.Gs[0,2,x,y,z-1])
                        )

        except:
            print('opps Dv0',x,y,z)
            DV[0]=0
        
        #calculate velocity in y based on 3.54
        try:
            if y==self.MaxY and x==0 and z==0:
                DV[1] = 0
            elif x==0 and z==0:
                DV[1] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y+1,z])) *
                        (self.Gs[1,1,x,y+1,z] - self.Gs[1,1,x,y,z])
                        )
            elif y==self.MaxY and z==0:
                DV[1] = ((1 / self.ds ) *
                        (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x-1,y,z])
                        )
            elif y==self.MaxY and x==0:
                DV[1] = ((1 / self.ds ) 
                        * (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x+1,y,z] + self.Gs[1,1,x,y-1,z] - self.Gs[1,1,x,y,z])
                        )
            elif y==self.MaxY:
                DV[1] = ((1 / self.ds ) *
                        (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x-1,y,z] + self.Gs[1,2,x,y,z] -self.Gs[1,2,x,y,z-1])
                        )
            elif x==0:
                DV[1] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y+1,z])) *
                        (self.Gs[1,1,x,y+1,z] - self.Gs[1,1,x,y,z] + self.Gs[1,2,x,y,z] -self.Gs[1,2,x,y,z-1])
                        )
            elif z==0:
                DV[1] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y+1,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x-1,y,z] + self.Gs[1,1,x,y+1,z] - self.Gs[1,1,x,y,z])
                        )
            else:
                DV[1] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y+1,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x-1,y,z] + self.Gs[1,1,x,y+1,z] - self.Gs[1,1,x,y,z] + self.Gs[1,2,x,y,z] -self.Gs[1,2,x,y,z-1])
                        )
        except:
            print('opps Dv1',x,y,z)
            DV[1]=0

        #calculate velocity in z based on 3.54
        try:
            if z==self.MaxZ and x==0 and y==0:
                DV[2] = 0
            elif x==0 and y==0:
                DV[2] = ( 
                        (1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z] + self.Gp[0,x,y,z+1])) *
                        (self.Gs[2,2,x,y,z+1] -self.Gs[2,2,x,y,z])
                        )
            elif z==self.MaxZ and y==0:
                DV[2] = ((1 / self.ds ) *
                        (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,2,x,y,z] - self.Gs[0,2,x-1,y,z])
                        )
            elif z==self.MaxZ and x==0:
                DV[2] = ( 
                        (1 / self.ds )
                        * (1 / (self.Gp[0,x,y,z])) 
                        * (self.Gs[1,2,x,y,z] - self.Gs[1,2,x,y-1,z])
                        )
            elif z==self.MaxZ:
                DV[2] = ((1 / self.ds ) *
                        (1 / (self.Gp[0,x,y,z])) *
                        (self.Gs[0,2,x,y,z] - self.Gs[0,2,x-1,y,z] + self.Gs[1,2,x,y,z] - self.Gs[1,2,x,y-1,z])
                        )
            elif x==0:
                DV[2] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y,z+1])) *
                        (self.Gs[1,2,x,y,z] - self.Gs[1,2,x,y-1,z] + self.Gs[2,2,x,y,z+1] -self.Gs[2,2,x,y,z])
                        )
            elif y==0:
                DV[2] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y,z+1])) *
                        (self.Gs[0,2,x,y,z] - self.Gs[0,2,x-1,y,z] + self.Gs[2,2,x,y,z+1] -self.Gs[2,2,x,y,z])
                        )
            else:
                DV[2] = ((1 / self.ds ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y,z+1])) *
                        (self.Gs[0,2,x,y,z] - self.Gs[0,2,x-1,y,z] + self.Gs[1,2,x,y,z] - self.Gs[1,2,x,y-1,z] + self.Gs[2,2,x,y,z+1] -self.Gs[2,2,x,y,z])
                        )

        except:
            print('opps Dv2',x,y,z, self.ds, self.MaxZ)
            DV[2]=0
        
        DV = self.CheckVelocityBoundary(x,y,z,DV)

        return DV
    
    def UpdateStresses(self, x,y,z):
        #Updates velocity based off of previous velocities and delta velocity times time
        #
        # Inputs: Coordinates of cube in question.  Assumed last time step
        #
        # Output: updated self.Gs matrix
        
        delS = self.DeltaStress(x,y,z)

        #for i in range(3):
        #        for j in range(3):
        #            self.Gs[i,j,x,y,z] += delS[i,j] * self.ts
        self.Gs[:,:,x,y,z] += delS[:,:] * self.ts

        return self

    def UpdateVelocity(self,x,y,z):
        #Updates velocity based off of previous velocities and delta velocity times time
        #
        # Inputs: Coordinates of cube in question.  Assumed last time step
        #
        # Output: updated self.Gs matrix

        delV = self.DeltaVelocity(x,y,z)

        for i in range(3):
            self.Gv[i,x,y,z] += delV[i] * self.ts

        return self
    
    def ForcingFunctionWave(self, t, Hz = 40000, EP=100):
        # Adds stresses from a force to the stress grid
        # Initially assumed a single force of a small plate sinosoidal ultrasound emitter.  More to be added later
        # 
        # Input:   t is the time
        #          Hz is the frequency of the signal
        #          EP is the force / density
        #
        # Outputs: no direct outputs, last time step stress is updated

        frequency = Hz
        EmitterPreasure = EP
        
        ##run for two periods and then stop:
        #if 2.0 / frequency < t:
            
        EmitterWidth = 0.01 / self.ds
        EmitterWidth = int(EmitterWidth)
        if EmitterWidth <= 3: EmitterWidth = 4

        #emitter placed in middle of top face
        
        StartX = int(self.MaxX / 2 - EmitterWidth / 2)
        StartZ = int(self.MaxZ / 2 - EmitterWidth / 2)

        Temp = np.zeros((EmitterWidth,EmitterWidth))

        Temp[:,:] = np.sin(frequency * t) * EmitterPreasure

        self.Gv[1,StartX:StartX+EmitterWidth,self.MaxY-2,StartZ:StartZ+EmitterWidth] += Temp * self.ts

        return self

    def ForcingFunctionImpulse(self, DeltaV = 100, emitter = 0.01, Dimm=1):
        # Adds stresses from a force to the stress grid
        # Initially assumed a single force of a small plate sinosoidal ultrasound emitter.  More to be added later
        # 
        # Input is the time
        #
        # Outputs: no direct outputs, last time step stress is updated

        EmitterWidth = emitter / self.ds
        EmitterWidth = int(EmitterWidth)
        if EmitterWidth <= 3: EmitterWidth = 4

        #emitter placed in middle of top face
        
        StartX = int((self.MaxX / 2) - (EmitterWidth / 2))
        StartZ = int((self.MaxZ / 2) - (EmitterWidth / 2))

        Temp = np.zeros((EmitterWidth,EmitterWidth))
        #print (EmitterWidth, StartX, StartZ)
        Temp[:,:] = -DeltaV

        self.Gv[Dimm,StartX:StartX+EmitterWidth,self.MaxY-2,StartZ:StartZ+EmitterWidth] += Temp 

        return self

    def VelocityCut(self, Dimm, location = -1):
        #Gets a cut of the magnitude of velocity in a plane defined by the two dimmension, cut at the location given in location
        #
        # Inputs: Dimm is the dimmension that is orthoganl to the plane that will be reutrned.
        #         Location if given will be the cut allong the thrid dimmension that will be returned.  
        #               If not given, -1 is default, it will signal using the centroid by default
        #
        # Outputs: 2 dimesional matrix with the combined velocity in the 2 given dimensions

        Results = []

        if Dimm == 0:
            if location > self.MaxX or location < 0: location = -1
            if location == -1:
                location = int(self.MaxX / 2)
            Component1 = np.matrix((self.MaxY, self.MaxZ))
            Component2 = np.matrix((self.MaxY, self.MaxZ))

            Component1 = self.Gv[1,location,:,:]
            Component2 = self.Gv[2,location,:,:]
            
            Results = Component1 + Component2
            #Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))

        if Dimm == 1:
            if location > self.MaxY or location < 0: location = -1
            if location == -1:
                location = int(self.MaxY / 2)
            Component1 = np.matrix((self.MaxX, self.MaxZ))
            Component2 = np.matrix((self.MaxX, self.MaxZ))

            Component1 = self.Gv[0,:,location,:]
            Component2 = self.Gv[2,:,location,:]
            
            Results = Component1 + Component2
            # Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))
            
        if Dimm == 2:
            if location > self.MaxZ or location < 0: location = -1
            if location == -1:
                location = int(self.MaxZ / 2)
            Component1 = np.matrix((self.MaxX, self.MaxY))
            Component2 = np.matrix((self.MaxX, self.MaxY))

            Component1 = self.Gv[0,:,:,location]
            Component2 = self.Gv[1,:,:,location]
            
            Results = Component1 + Component2
            #Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))
            

        return Results

    def StressCut(self, Dimm, Dimm1, Dimm2, location = -1):
        #Gets a cut of the magnitude of stress on the give face and direction, cut at the location given in location
        #
        # Inputs: Dimm is the dimmension to make the cut across
        #         Dimm1 is the dimmension that is orthoganl to the plane that will be reutrned.
        #         Dimm2 is the vector of the stress
        #         Location if given will be the cut allong the thrid dimmension that will be returned.  
        #               If not given, -1 is default, it will signal using the centroid by default
        #
        # Outputs: 2 dimesional matrix with the combined velocity in the 2 given dimensions

        Results = []

        if Dimm == 0:
            if location > self.MaxX or location < 0: location = -1
            if location == -1:
                location = int(self.MaxX / 2)
            Component1 = np.matrix((self.MaxY, self.MaxZ))
            
            Component1 = self.Gs[Dimm1,Dimm2,location,:,:]
            
            Results = Component1
            #Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))

        if Dimm == 1:
            if location > self.MaxY or location < 0: location = -1
            if location == -1:
                location = int(self.MaxY / 2)
            Component1 = np.matrix((self.MaxX, self.MaxZ))
            
            Component1 = self.Gs[Dimm1, Dimm2,:,location,:]
            
            Results = Component1
            # Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))
            
        if Dimm == 2:
            if location > self.MaxZ or location < 0: location = -1
            if location == -1:
                location = int(self.MaxZ / 2)
            Component1 = np.matrix((self.MaxX, self.MaxY))
            
            Component1 = self.Gs[Dimm1,Dimm2,:,:,location]
            
            Results = Component1
            #Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))
            

        return Results