# -*- coding: utf-8 -*-
"""
Created on Thu Mar. 28 12:15:27 2019
@author: shin
ファイルデータのヒストグラムから分布生成
ヒストグラムから生成した分布で四則演算
"""


from bunpu import *



lists0=['nikkei']
lists1=['toyota']
lists2=['toyota_sony']
lists3=['toyota_sony_nikkei']

x1=[bunpu(),bunpu(),bunpu(),bunpu()]
x2=[bunpu(),bunpu(),bunpu(),bunpu()]
x3=[bunpu(),bunpu(),bunpu(),bunpu()]

infile2=[lists2[0]+'.txt']
infile3=[lists3[0]+'.txt']
outfile0=lists0[0]
outfile1=lists1[0]
outfile2=lists2[0]
outfile3=lists3[0]

nik=bunpu()
nik.bunpu_data(infile3,outfile0,2,[17],[20],0,0)
toyo=bunpu()
toyo.bunpu_data(infile3,outfile1,2,[11],[20],0,0)
toyoso=bunpu()
toyoso.bunpu_data(infile2,outfile2,2,[5,11],[20,20],0,0)
toyosonk=bunpu()
toyosonk.bunpu_data(infile3,outfile3,2,[5,11,17],[20,20,20],0,0)

x1[0]=nik+toyo
x1[0].bunpu_graph('c+t')

x1[1]=nik-toyo
x1[1].bunpu_graph('c-t')

x1[2]=nik*toyo
x1[2].bunpu_graph('c＊t')

x1[3]=nik/toyo
x1[3].bunpu_graph('c／t')

x2[0]=toyoso+toyoso
x2[0].bunpu_graph('ts+ts')

x2[1]=toyoso-toyoso
x2[1].bunpu_graph('ts-ts')

x2[2]=toyoso*toyo
x2[2].bunpu_graph('ts＊t')

x2[3]=toyoso/toyo
x2[3].bunpu_graph('ts／t')

x3[0]=toyosonk+toyosonk
#x3[0].bunpu_graph('tsn+tsn')
print('matplotlib2.2.2より新しいとsavefigを行うと表示が壊れるので表示のみとします、2.2.2ならコメントを外してください')
x3[0].bunpu_graph()
x3[0].bunpu_file('tsn+tsn.csv')

x3[1]=toyosonk-toyosonk
#x3[1].bunpu_graph('tsn-tsn')
x3[0].bunpu_graph()
x3[1].bunpu_file('tsn-tsn.csv')

x3[2]=toyosonk*toyo
#x3[2].bunpu_graph('tsn＊t')
x3[0].bunpu_graph()
x3[2].bunpu_file('tsn＊t.csv')

x3[0]=toyosonk/toyo
#x3[0].bunpu_graph('tsn／t')
x3[0].bunpu_graph()
x3[0].bunpu_file('tsn／t.csv')
        
