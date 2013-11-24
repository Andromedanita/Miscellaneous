import numpy as np
import matplotlib.pyplot as plt
import scipy

a=0.1
b=0.1
c=6.5

delta_t=0.01
time_array= np.arange(0,200,delta_t)
nt=len(time_array)

loops=np.arange(4,15,1)

x = np.zeros(nt)
y =np.zeros(nt)
z= np.zeros(nt)

x[0]=10.0
y[0]=10.0
z[0]=10.0

def getmaxes(z,c):
    maxarray=[]
    carray=[]
    a=(np.diff(np.sign(np.diff(z))) < 0).nonzero()[0] + 1
    for q in range(0,np.size(a)):
        maxarray.append(z[a[q]])
        carray.append(c)
    return [maxarray,carray]



maximas=[]
cs=[]

for c in loops:
    i=0
    while i < nt-1:
        x[i+1]=((-y[i]-z[i])*delta_t)+x[i]
        y[i+1]=((x[i]+(a*y[i]))*delta_t)+y[i]
        z[i+1]=((b+(z[i]*(x[i]-c)))*delta_t)+z[i]
        i+=1

    tempmaxes=getmaxes(z,c)
    maximas.extend(tempmaxes[0])
    cs.extend(tempmaxes[1])
    print 'cs values are:', cs, 'z values are:', maximas
'''
plt.subplot(3,1,1)
plt.plot(time_array,x)
plt.ylabel('x(t)')
plt.grid('on')

plt.subplot(3,1,2)
plt.plot(time_array,y)
plt.ylabel('y(t)')
plt.grid('on')

plt.subplot(3,1,3)
plt.plot(time_array,z)
plt.ylabel('z(t)')
plt.grid('on')

plt.suptitle('x(t) , y(t) , z(t)')
plt.xlabel('time')
# savefig('shm_energy_output.pdf')
    '''
plt.plot(cs,maximas,'o',markersize=1.5)
plt.show()
