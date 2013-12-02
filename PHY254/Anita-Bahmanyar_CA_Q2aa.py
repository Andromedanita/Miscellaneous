import numpy as np
import matplotlib.pylab as plt

time=np.arange(0,60,0.01)
nt=len(time)


theta=np.zeros(nt)
v=np.zeros(nt)

theta[0]=1.0
v[0]=1.0

i=0
while i<nt-1:
    v[i+1]=v[i]+(0.05*(-np.sin(theta[i])))
    theta[i+1]=theta[i]+(0.05*v[i])
    i+=1

energy=(1-(np.cos(theta)))+(0.5*(v**2))

plt.subplot(2,1,1)
plt.plot(time,theta)
plt.ylabel('theta')
plt.xlabel('time')
plt.grid('on')
plt.subplot(2,1,2)
plt.plot(time,energy)
plt.ylabel('energy')
plt.xlabel('time') 
plt.grid('on')
plt.suptitle('Theta and Energy vs. Time using forward Euler method')
plt.show()
