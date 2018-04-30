import numpy as np
import matplotlib.pyplot as plt


class Model():

    def __init__(self):
        self.t0 = 0.0
        self.dt = 1e-4
        self.Nstep = 10000
        self.Nparticle = 10000
        self.Nparticle = 10
        self.Nstep = 100

        
        self.u = np.zeros((self.Nparticle,self.Nstep),dtype=float)
        self.t = np.arange(self.Nstep)*self.dt

        ## Statistical Quantities
        self.mean = np.zeros(self.Nstep,dtype=float)
        self.std = np.zeros(self.Nstep,dtype=float)
        
    def Initialize(self):
        print("Initialization")
        self.u[:,0] = np.zeros(self.Nparticle)
                
    def Solve(self):

        print("Wienner process")

        dW = np.random.normal(loc=0.0,scale = np.sqrt(self.dt),size=(self.Nparticle,self.Nstep))
        
        for i in range(self.Nstep-1):
            self.u[:,i+1] = self.u[:,i] + dW[:,i]
                            
    def StatisticalQuantities(self):

        self.mean = np.mean(self.u,axis=0)
        self.std = np.std(self.u,axis=0)


    def Display(self):

        
        #plt.plot(self.t,self.u[0,:])
        plt.plot(self.t,self.mean,label="mean")
        for i in range(3):
            plt.plot(self.t,self.u[i,:],label="particle {0}".format(i))
        
        #plt.errorbar(self.t,self.mean,yerr=self.std)

        plt.legend()
        plt.savefig("Path.png")
        plt.show()

    def Exact(self,x):
        var = self.t[-1]
        return np.exp(-(x**2/(2.0*var)))/(np.sqrt(2.0*np.pi*var))
        
    def Histogram(self):
        
        zmax = np.max(np.abs(self.u[:,-1]))
        z = np.linspace(-zmax*1.1,zmax*1.1,1000)

        plt.hist(self.u[:,-1],bins = 100,normed='True',label="simulation")
        plt.plot(z,self.Exact(z),label="Exact")
        plt.legend()
        plt.savefig("Histogram.png")
        plt.show()
    



class OrnsteinUhlenbeck(Model):

    def __init__(self):
        super().__init__()
        self.theta = 10.0        
        self.mu = 0.0
        self.sigma = 0.1
        
    def Initialize(self):
        print("Initialization")
        self.u[:,0] = np.ones(self.Nparticle)
                        
    def Solve(self):
        print("Ornstein-Uhlenbeck process")
        dW = np.random.normal(loc=0.0,scale = np.sqrt(self.dt),size=(self.Nparticle,self.Nstep))
        
        for i in range(self.Nstep-1):
            self.u[:,i+1] = self.u[:,i] + -1.0*self.dt*self.theta*(self.u[:,i]-self.mu)+self.sigma*dW[:,i]
            
    def Exact(self,x):
        D = self.sigma**2/2.0
        return np.exp(-self.theta*(x-self.mu)**2/(2.0*D))*(np.sqrt(self.theta/(2.0*np.pi*D)))



class KPZ(Model):

    def __init__(self):
        super().__init__()
        self.Lx = 1000
        self.dx = 1e-2
        self.nu = 1.0
        self.lambd = 1.0
        self.sigma = 1.0
        self.u = np.zeros((self.Nparticle,self.Nstep,self.Lx),dtype=float)
        self.mean = np.zeros((self.Nstep,self.Lx),dtype=float)
        self.std  = np.zeros((self.Nstep,self.Lx),dtype=float)
        
    def Initialize(self):
        print("Initialization")
        self.u[:,0,:] = np.zeros((self.Nparticle,self.Lx))
                        
    def BoundaryCondition(self,i):

        self.u[:,i,0]  = 0.0
        self.u[:,i,-1] = 0.0

    def Solve(self):
        print("KPZ equation")

        dW = np.random.normal(loc=0.0,scale = np.sqrt(self.dt),size=(self.Nparticle,self.Nstep,self.Lx))
        
        for i in range(self.Nstep-1):
            self.BoundaryCondition(i)
            self.u[:,i+1,1:-1] \
                = self.u[:,i,1:-1] \
                + self.sigma*dW[:,i,1:-1] \
                +self.dt*self.nu*((self.u[:,i,2:]-2*self.u[:,i,1:-1]+self.u[:,i,:-2])/(self.dx**2))\
                +self.dt*self.lambd*0.5**(((self.u[:,i,2:]-self.u[:,i,:-2])/self.dx)**2
                )
            
    def Display(self):

        elem = int(self.Lx/2)
#        plt.plot(self.t,self.mean[:,elem],label="mean")
        for i in range(3):
            for j in range(self.Lx):
                plt.plot(self.t,self.u[i,:,j],label="particle {0}".format(i))
        
        #plt.errorbar(self.t,self.mean,yerr=self.std)

        plt.legend()
        plt.savefig("Path.png")
        plt.show()

    def Exact(self,x):
        return np.zeros_like(x)

        
    def Histogram(self):
        
        zmax = np.max(np.abs(self.u[:,-1]))
        z = np.linspace(-zmax*1.1,zmax*1.1,1000)

        plt.hist(self.u[:,-1,elem],bins = 100,normed='True',label="simulation")
        plt.plot(z,self.Exact(z),label="Exact")
        plt.legend()
        plt.savefig("Histogram.png")
        plt.show()
    


            
if __name__ == "__main__":

    Model = KPZ()
    #Model = OrnsteinUhlenbeck()
    #Model = Model()
    Model.Initialize()
    #Model.Solve()
    Model.Solve()
    Model.StatisticalQuantities()
    Model.Display()
    Model.Histogram()
