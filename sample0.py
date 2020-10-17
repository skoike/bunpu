# -*- coding: utf-8 -*-
"""
Created on Thu Mar. 28 12:15:27 2019
@author: shin
"""


from bunpu import *


x1=bunpu()
x1.bunpu_gene([-10.0],[15],[9],[2],[100],'x1')
x1.bunpu_graph('1d x')

x2=bunpu()
x2.bunpu_gene([-1.0,2.0],[0.5,4.0],[-0.1,2.8],[0.2,0.3],[20,20],'x2')
x2.bunpu_graph('2d x')
x2.bunpu_file('2d x')
x3=bunpu()
x3.bunpu_gene([-5.0,-5.0,5.0],[10.0,10.0,25.0],[4,1,13],[2,2,2],[20,20,20],'x3')
x3.bunpu_graph('3d x')
x3.bunpu_file('3d x')

#
y1=bunpu()
y1.bunpu_gene([10],[25],[13],[2],[20],'y1')
y1.bunpu_graph('1d y')

y2=bunpu()
y2.bunpu_gene([1.0,2.0],[3.5,4.0],[2.2,2.8],[0.3,0.2],[20,20],'y2')
y2.bunpu_graph('2d y')
y2.bunpu_file('2d y')

y3=bunpu()
y3.bunpu_gene([30,10,20],[55,35,40],[39,22,28],[2,3,2],[20,20,20],'y3')
y3.bunpu_graph('3d y')
y3.bunpu_file('3d y')

z1=bunpu()
z1=x1+y1
z1.bunpu_graph('1dx+y')

z1=bunpu()
z1=x1-y1
z1.bunpu_graph('1d x-y')

z1=bunpu()
z1=x1*y1
z1.bunpu_graph('1d x＊y')

z1=bunpu()
z1=x1/y1
z1.bunpu_graph('1d x／y')

z2=bunpu()
z2=x2+y2
z2.bunpu_graph('2d x+y')

z2=bunpu()
z2=x2-y2
z2.bunpu_graph('2d x-y')

z2=bunpu()
z2=x2*y1
z2.bunpu_graph('2d x＊y')
z2.bunpu_file('2d x＊y')

z2=bunpu()
z2=x2/y1
z2.bunpu_graph('2d x／y')
z2.bunpu_file('2d x／y')

z3=bunpu()
z3=x3+y3
z3.bunpu_graph('3d x+y')

z3=bunpu()
z3=x3-y3
z3.bunpu_graph('3d x-y')

z3=bunpu()
z3=x3*y1
z3.bunpu_graph('3d x＊y',1)
z3.bunpu_file('3d x＊y')

z3=bunpu()
z3=x3/y1
z3.bunpu_graph('3d x／y',1)
z3.bunpu_file('3d x／y')

