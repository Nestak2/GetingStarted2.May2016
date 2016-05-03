import matplotlib.pyplot as plt
import numpy as np

def NestorSineWave(lam):
    x = np.arange(-2*np.pi, 2*np.pi, 0.01)
    y = np.sin(2*np.pi/lam * x)
    plt.plot(x, y)
    pl.grid(True)
    print x
    plt.show()
   
    
    
NestorSineWave(0.5*2*np.pi)
