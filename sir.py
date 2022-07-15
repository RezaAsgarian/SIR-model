import numpy as np
import matplotlib.pyplot as plt

def f(xi,yi,**par):
    si,ii = yi[0],yi[1]
    beta = par["beta"]
    gamma = par["gamma"]
    N = par["N"]
    f_ = [
    - beta * si * ii / N,
    + beta * si * ii / N - gamma * ii
    ]
    return f_
def k1(xi,yi,**par):
    return f(xi,yi,**par)
def k2(xi,yi,**par):
    h = par["h"]
    x_ = xi + h/2
    y_ = [yi_ + k1_*h/2 for yi_,k1_ in zip(yi,k1(xi,yi,**par))]
    return f(x_,y_,**par)
def k3(xi,yi,**par):
    h = par["h"]
    x_ = xi + h/2
    y_ = [yi_ + k2_*h/2 for yi_,k2_ in zip(yi,k2(xi,yi,**par))]
    return f(x_,y_,**par)
def k4(xi,yi,**par):
    h = par["h"]
    x_ = xi + h
    y_ = [yi_ + k3_*h for yi_,k3_ in zip(yi,k3(xi,yi,**par))]
    return f(x_,y_,**par)
def yi1(xi,yi,**par):
    h = par["h"]
    yi1_ = [yi_ + (h/6) * (k1_+2*k2_+2*k3_+k4_)
    for yi_,k1_,k2_,k3_,k4_ in zip(yi,k1(xi,yi,**par),k2(xi,yi,**par),
    k3(xi,yi,**par),k4(xi,yi,**par))
    ]
    return yi1_
N = 80
i0 = (0.01/100.0)*N
s0 = N-i0
gamma = 10.0/100.0
beta = 50.0/100.0
h = 0.1
T = 100
Nsteps = int(T/h)
yy = [[s0,i0]]
tt = [0]
for i in range(1,Nsteps):
    y_ = yy[-1]
    t_ = tt[-1]
    t_new = t_ + h
    y_new = yi1(t_,y_,h=h,gamma=gamma,beta=beta,N=N)
    tt.append(t_new)
    yy.append(y_new)
S = [i[0] for i in yy]
I = [i[1] for i in yy]
R = [N - i[0] - i[1] for i in yy]
t = tt
fig = plt.figure()
ax = fig.add_axes((0.15, 0.15, 0.7, 0.7))
plt.plot(t,S,"b-",label="Susceptible")
plt.plot(t,I,"r-",label="Infectious")
plt.plot(t,R,"g-",label="Recovered")
plt.xlabel('time (days)')
plt.ylabel('population (million persons)')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
plt.legend()
plt.show()
