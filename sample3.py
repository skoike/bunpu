# -*- coding: utf-8 -*-


"""
拡散方程式
水分子間距離:1*10-7mm(0.96Å)
"""
from bunpu import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML
 
nx = 41

dx = 2 / (nx-1)
nt = 25 
nu = 0.3 
alpha =0.3  
dt = alpha* dx**2 /nu
c = nu * dt / dx**2
print('dx=',dx)
print('dt=',dt)
print('c=',c)
div=[30]
x = np.linspace(0,2,nx)
col1=[0]
col2=[0]
u=[]
shw=1
for i in range(nx):
    u.append(bunpu())

for i in range(nx):
    if i<int(nx/4) or i>int(3*nx/4):
        u[i].bunpu_gene([0.99999999],[1.00000001],[1.0],[0.000000002],div,'u'+str(i))
    else:
        u[i].bunpu_gene([1.99999999],[2.00000001],[2.0],[0.000000002],div,'u'+str(i))
 
fig = plt.figure(figsize=(8,4))
ims=[]
for n in range(nt):
    un = []
    umax=[]
    umin=[]
    print(n)
    for i in range(nx):
        un.append(u[i])
        umax.append(u[i].xmax)
        umin.append(u[i].xmin)

    if (nt%1==0):
        im = plt.plot(x,umin, "r",x,umax, "r",)
        ims.append(im)
    for i in range(1, nx - 1):
        uu=[bunpu(),bunpu(),bunpu()]
        uu[0]=un[i+1].bunpu_sub(un[i],div,col1,shw)
        uu[1]=uu[0].bunpu_sub(un[i],div,col1,shw)
        uu[2]=uu[1].bunpu_add(un[i-1],div,col1,shw)
        u[i]=un[i].bunpu_simu(uu[2],c,col2,div,shw)
plt.grid()
plt.ylim([-0.1,3])
plt.xlabel("x")
plt.ylabel("u(m/s)")     
anim = animation.ArtistAnimation(fig, ims)
rc('animation', html='jshtml')
plt.show()
anim
