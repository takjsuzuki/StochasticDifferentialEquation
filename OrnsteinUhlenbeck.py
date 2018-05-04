import numpy as np
import matplotlib.pyplot as plt

import Model

class OrnsteinUhlenbeck(Model.Model):

    def __init__(self):
        super().__init__()
        self.theta = 5.0        
        self.mu = 0.0
        self.sigma = 0.2
        self.InitialValue = 1.0
        
    def Initialize(self):
        print("Initialization")
        self.u[:,0] = np.ones(self.Nparticle)*self.InitialValue
                        
    def Solve(self):
        print("Ornstein-Uhlenbeck process")
        dW = np.random.normal(loc=0.0,scale = np.sqrt(self.dt),size=(self.Nparticle,self.Nstep))
        
        for i in range(self.Nstep-1):
            self.u[:,i+1] \
                = self.u[:,i] + -1.0*self.dt*self.theta*(self.u[:,i]-self.mu)+self.sigma*dW[:,i]
            
    def Exact(self,x,time):
        D = self.sigma**2/2.0
        return np.exp(-self.theta*(x-self.mu - (self.InitialValue-self.mu)*np.exp(-self.theta*time))**2/(2.0*D*(1.0-np.exp(-2.0*self.theta * time)))) * (np.sqrt(self.theta/(2.0*np.pi*D*(1.0-np.exp(-2.0*self.theta * time)))))



