import numpy as np
import matplotlib.pyplot as plt
import scipy

a=0.1
b=0.1
c=4

a1=0.1
b1=0.1
c1=10

delta_t=0.001
time_array= np.arange(0,200,delta_t)
time_array1= np.arange(0,200,delta_t)
nt=len(time_array)

x = np.zeros(nt)
y =np.zeros(nt)
z= np.zeros(nt)

x1 = np.zeros(nt)
y1 =np.zeros(nt)
z1= np.zeros(nt)

dif=np.zeros(nt)


x[0]=10.0
y[0]=10.0
z[0]=10.0

x1[0]=10.0001
y1[0]=10.0
z1[0]=10.0

dif[0]=0


i=0
while i < nt-1:
    x[i+1]=((-y[i]-z[i])*delta_t)+x[i]
    y[i+1]=((x[i]+(a*y[i]))*delta_t)+y[i]
    z[i+1]=((b+(z[i]*(x[i]-c)))*delta_t)+z[i]
    x1[i+1]=((-y1[i]-z1[i])*delta_t)+x1[i]
    y1[i+1]=((x1[i]+(a1*y1[i]))*delta_t)+y1[i]
    z1[i+1]=((b1+(z1[i]*(x1[i]-c1)))*delta_t)+z1[i]
    dif[i+1]=np.abs(x1[i+1]-x[i+1])
    i+=1


plt.subplot(4,1,1)
plt.plot(time_array,x,'r',label='x(t)')
plt.plot(time_array1,x1,'b',label='x1(t)')
plt.ylabel('x(t)')
plt.grid('on')
plt.legend()

plt.subplot(4,1,2)
plt.plot(time_array,y,'r')
plt.plot(time_array1,y1,'b')
plt.ylabel('y(t)')
plt.grid('on')

plt.subplot(4,1,3)
plt.plot(time_array,z,'r')
plt.plot(time_array1,z1,'b')
plt.ylabel('z(t)')
plt.grid('on')

plt.subplot(4,1,4)
plt.plot(time_array,dif,'g')
plt.ylabel('x(t) difference')
plt.grid('on')

plt.suptitle('x(t) , y(t) , z(t)')
plt.xlabel('time')
# savefig('shm_energy_output.pdf')
plt.show()
