from pylab import *
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from numpy import *

def rhs(xvector,t):

	x1dot=????
	x2dot=????	
	x3dot=????
	x4dot=????

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
sigma=.25

# My IC's
x0=1
y0=1
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
savefig("SpringPendulumtest.pdf")
show()