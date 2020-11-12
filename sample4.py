# -*- coding: utf-8 -*-
"""
Created on Thu Mar. 28 12:15:27 2019
@author: shin
"""
'''
放物線運動
'''

from bunpu import *



v=bunpu()
v.bunpu_gene([1.0,3.0,4.0],[1.2,3.3,4.3],[1.12,3.2,4.2],[0.03,0.03,0.03],[20,20,20],'v')
v.bunpu_graph('v')

x=bunpu()
x.bunpu_gene([0.0,0.0,0.0],[0.05,0.05,0.05],[0.03,0.03,0.03],[0.005,0.005,0.005],[20,20,20],'x')
tend=11
dt=0.02

core1=[0,0,0]
core2=[0,0,0]
core3=[0,0,0]
n=100
vdiv=[5,5,5]
posx=[]
posy1=[]
posy2=[]
posymx=[]
time=[]
vel1=[]
vel2=[]
r=[0.0]
for t in range(tend):
    g=bunpu()
    g.bunpu_gene([0,0,-20],[2,2,-18],[1.2,1.2,-19.2],[0.3,0.3,0.2],[20,20,20],'cons')
    v1=bunpu()
    v1=v.bunpu_simu(g,dt,core2,vdiv)
    v=bunpu()
    v=v1
    x1=bunpu()
    x1=x.bunpu_simu(v,dt,core3,vdiv)
    x=bunpu()
    x=x1
    xlist=(x1.flatten[3]).tolist()
    num=xlist.index(max(xlist))
    posx.append(x1.flatten[1][num])
    posy1.append(max(x1.flatten[2]))
    posy2.append(min(x1.flatten[2]))
    posymx.append(x1.flatten[2][num])
    time.append(t)
    vel1.append(max(v.flatten[2]))
    vel2.append(min(v.flatten[2]))
    q,mod=divmod(t,5)
    if t>=0 and mod==0:
        print('time',t)
        g.bunpu_graph('a'+str(t))
        x.bunpu_graph('x'+str(t))
        if t>2:
            v1.bunpu_graph('v'+str(t))
    #
fig1=plt.figure(figsize=(14,7))
ax1 = fig1.add_subplot(111)
ax1.plot(posx,posy1, color="g")
ax1.plot(posx,posy2, color="b")
ax1.plot(posx,posymx, color="r")
plt.savefig("kiseki")
fig2=plt.figure(figsize=(14,7))
ax2 = fig2.add_subplot(111)
ax2.plot(time,vel1, color="g")
ax2.plot(time,vel2, color="b")
plt.savefig("vel")
            


    
