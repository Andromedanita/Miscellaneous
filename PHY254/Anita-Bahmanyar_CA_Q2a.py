import numpy as np
import matplotlib.pylab as plt

time=np.arange(0,60,0.05)
nt=np.len(time)


theta=np.zeros(nt)
v=np.zeros(nt)

theta[0]=1.0
v[0]=1.0

while i < nt-1:


        


#theta vs time plot
plt.plot(time,theta)
#energy vs time plot
plt.plot(time,energy)
plt.title('')
plt.xlabel('')
plt.ylabel('')
plt.show()
