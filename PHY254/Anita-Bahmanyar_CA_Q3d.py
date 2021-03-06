from pylab import *
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from numpy import *

def rhs(xvector,t):
        x,vx,y,vy= xvector
	#x1dot=vx
	#x2dot=(((-1.0/sigma)*x)+((1.0/sigma)*(1-sigma)*x)/(np.sqrt((x**2)-((1-y)**2))))
	#x3dot=vy
	#x4dot=(-1)+((1.0/sigma)*(1-y))-((1.0/sigma)*(1-sigma)*(1-y)*(1.0/(np.sqrt((x**2)-((1-y)**2)))))
        x1dot=vx
	x2dot=(-(1./sigma)*x)+(((1./sigma)-1)*(x/(np.sqrt((x**2)+((1-y)**2)))))
	x3dot=vy
	x4dot= -1 + (1./sigma)*(1-y) - (((1./sigma)-1)*((1-y)/np.sqrt((x**2)+((1-y)**2))))

	return [x1dot,x2dot,x3dot,x4dot]

# Define the starttime, endtime and when we want odeint
# to output our solution
start=0
end=50.0
num=200

# generate a list of times to output the solution at,
# "num" of them spaced equally going from "start", to "end"
t=linspace(start,end,num)

# Our nondimensional parameter
sigma=0.0001

# My IC's
x0=0.8
y0=0.4
vx0=0.0
vy0=0.0

# Call odeint
solution=odeint(rhs,[x0,vx0,y0,vy0],t)

# Extract the solution
x=solution[:,0]
vx=solution[:,1]
y=solution[:,2]
vy=solution[:,3]

# Plot the solution in subplots
subplot(3,1,1)
plot(t,x)
xlabel('t')
ylabel('x')

subplot(3,1,2)
plot(t,y) 
xlabel('t')
ylabel('y')

subplot(3,1,3, aspect='equal')
plot(x,y)
xlabel('x')
ylabel('y')
#savefig("SpringPendulumtest.pdf")
plt.suptitle('x and y values of the pendulum')
show()
