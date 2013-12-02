from pylab import *
from scipy.integrate import odeint

# This function is the right hand side of the system of equations 
# solved by odeint. Given the array "state" which contains theta
# and thetadot we return the time derivatives of these

def rhs(xvector,t):

	x1dot=xvector[1]
	x2dot=-sin(xvector[0])
	
	return [x1dot, x2dot]
	

start=0.0
end=60.0
num=400

#Get an array of times that we want outputted
t=linspace(start,end,num)

# Our initial conditions
theta0=3.13
thetadot0=0.0

# Call odeint with the name of our right hand side function
# and an array of the initial conditions in the same order
# the rhs function expects them in
solution=odeint(rhs,[theta0,thetadot0],t)

# Get theta and thetadot out of the solution (returned by odeint)
# and calculate the energy
theta=solution[:,0]
thetadot=solution[:,1]
energy=1-cos(theta)+.5*pow(thetadot,2)

# Plot both theta and the energy
figure()
subplot(2,1,1)
plot(t,theta)
xlabel('t')
ylabel('theta')
subplot(2,1,2)
plot(t,energy)
xlabel('t')
ylabel('Energy')

show()
