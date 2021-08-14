# -*- coding: utf-8 -*-
"""
© 2020 Shin Koike  bunpu@a1.rim.or.jp

このソフトウェアをそのままの複製を学習や研究を目的として、利用する場合、本ソフトウェアおよび今後作成されるものを含めたそのブランチの利用を無償で許可します。

このソフトウェアは未完成で、改善の提案や機能拡張の協力を求めています、このソフトの改善や協力の為にに、変更、追加、結合、移植を含む派生を、
利用可能な情報とともに、公開を前提として、前記アドレスにその情報提供をお願いします。
その内容は公共性に基づいて本ソフトまたはそのブランチに反映させていきます。

このソフトを利用・参考にする場合は、このソフトの著作権と特許出願（PCT/JP2020/034566とその分割、関連出願）およびその協力者における権利を尊重ください。
このソフトウェアの一部分を利用または参考にして、変更、追加、結合、継承や移植を含む派生を、配布または商用利用する場合は前記アドレスに相談してください。

ソフトウェアは、未完成で、何らの保証もなく提供されます。
ここでいう保証とは、商品性、特定の目的への適合性、および権利非侵害についての保証も含みますが、それに限定されるものではありません。 
このソフト作者または著作権者は、契約行為、不法行為、またはそれ以外であろうと、ソフトウェアに起因または関連し、あるいはソフトウェアの使用または
その他の扱いによって生じる一切の請求、損害、その他の義務について何らの責任も負わないものとします。

以上の表示および本許諾表示を、ソフトウェアのすべての複製または部分の利用または分布処理を参考とする場合に、作成される著作物に記載するものとします。

© 2020 Shin Koike  bunpu@a1.rim.or.jp

Permission is hereby granted, free of charge, to any person obtaining a exact copy of this software,
its branches and associated documentation files (the "Software"),for learning or research purposes, to deal in the Software with restriction.

This software is incomplete and we are seeking suggestions for improvement and cooperation in enhancements.
For the improvement and cooperation of this software, please provide the derivation
including modification, addition, mergers,combination, translation with available information to above address
 the assumption that it will be published.
The contents will be reflected in this software and its branches based on public nature and my leeway.

When using or referring to this software, please correspond the copyright of this software
and the rights in patent applications(PCT/JP2020/034566 and divisional other).
Please contact with above address if you want to use or refer to a part of this software and distribute it privately or use it for commercial purposes.

The software is incomplete and is provided without warranty.Warranties here include, but are not limited to, warranties of merchantability, 
fitness for a particular purpose, and non-infringement.
The author or copyright holder of this software, whether contractual, tort, or otherwise, is due to or related to the software, or uses or uses the software.
We shall not be liable for any claims, damages or other obligations arising from any other dealings.

The above copyright notice and this permission notice shall be included in all copies, portions or reference of the software related to dstribution .


機能:分布ベクトル演算処理

・ファイルから読み取ったデータのヒストグラムから分布を作成
・範囲や平均値、標準偏差などを指定して分布を生成する
・1次元から3次元までの分布ベクトルの四則演算（相関係数含む）
・分布ベクトルの微小積分
・分布のグラフ表示
・分布パラメータのファイル出力

"""
#from scipy import stats
from scipy.stats import multivariate_normal,matrix_normal
from scipy.stats import powerlognorm,norminvgauss
import numpy as np
#import cupy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal, interpolate
from matplotlib import cm
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
import glob
import re
import csv
import os
import pandas as pd
import copy
import matplotlib.colors as colors
#from sklearn.neighbors.kde import KernelDensity#vr.0.22
from sklearn.neighbors import KernelDensity#vr.0.24
from scipy.stats import gaussian_kde
from matplotlib.font_manager import FontProperties



def genes(xmin,xmax,xmean,xdev,div=200,filename=''):#
    fp = FontProperties(fname=r'C:\Windows\Fonts\msgothic.ttc',size=24)
    nm=0
    ratio=1/(xmax-xmin)
    mean=(xmean-xmin)/(xmax-xmin)

    s=xdev*ratio
    dx=(xmax-xmin)/(div-1)

    m0=[[0.973736201,0.934747739,0.78264355,0.659429535,0.593065765,0.5,0.406934235,0.340570465,0.21735645,0.158253011,0.026263799],
        [0.96548985,0.82750139,0.76916612,0.65155935,0.588782244,0.5,0.411217756,0.34844065,0.23083388,0.17249861,0.03451015],
        [0.923366197,0.735090976,0.686755062,0.602139917,0.5595659,0.5,0.4404341,0.397860083,0.313244938,0.264909024,0.076633803],
        [0.868392295,0.666750228,0.630089404,0.570106766,0.540730027,0.5,0.459269973,0.429893234,0.369910596,0.333249772,0.131607705],
        [0.847765623,0.740479408,0.697816791,0.614424967,0.569021281,0.5,0.430978719,0.385575033,0.302183209,0.259520592,0.152234377],
        [0.8944794,0.75360135,0.707325691,0.621018408,0.574224197,0.5,0.425775803,0.378981592,0.292674309,0.24639865,0.1055206],
        [0.851113717,0.683406566,0.648475292,0.586926827,0.553495272,0.5,0.446504728,0.413073173,0.351524708,0.316593434,0.148886283],
        [0.795705603,0.634927311,0.608449192,0.562767704,0.538215521,0.5,0.461784479,0.437232296,0.391550808,0.365072689,0.204294397]]

    pk0=[[0.984924623,0.949748744,0.804020101,0.673366834,0.603015075,0.502512563,0.396984925,0.326633166,0.195979899,0.135678392,0.015075377],
         [0.994974874,0.904522613,0.844221106,0.698492462,0.618090452,0.497487437,0.381909548,0.301507538,0.155778894,0.095477387,0.005025126],
         [0.979899497,0.753768844,0.698492462,0.603015075,0.557788945,0.502512563,0.442211055,0.396984925,0.301507538,0.246231156,0.020100503],
         [0.91959799,0.809045226,0.75879397,0.648241206,0.587939698,0.497487437,0.412060302,0.351758794,0.24120603,0.190954774,0.08040201],
         [0.984924623,0.884422111,0.824120603,0.683417085,0.608040201,0.497487437,0.391959799,0.316582915,0.175879397,0.115577889,0.015075377],
         [0.984924623,0.743718593,0.688442211,0.59798995,0.557788945,0.497487437,0.442211055,0.40201005,0.311557789,0.256281407,0.015075377]]

    s0= [[0.048706732,0.06493023,0.09245454,0.102939157,0.1053216,0.10357657,0.1053216,0.102939157,0.09245454,0.084491614,0.048706732],
         [0.080432844,0.135852903,0.14445847,0.153933144,0.155733199,0.153522145,0.155733199,0.153933144,0.14445847,0.135852903,0.080432844],
         [0.117606404,0.172616385,0.177265418,0.180906322,0.180824964,0.177982922,0.180824964,0.180906322,0.177265418,0.172616385,0.117606404],
         [0.144828966,0.185096647,0.187140374,0.187919868,0.187152342,0.184536877,0.187152342,0.187919868,0.187140374,0.185096647,0.144828966],
         [0.136676644,0.166987943,0.175448815,0.185569018,0.186535589,0.179118666,0.186535589,0.185569018,0.175448815,0.166987943,0.136676644],
         [0.138952735,0.186754717,0.195625697,0.205495282,0.206341314,0.197554193,0.206341314,0.205495282,0.195625697,0.186754717,0.138952735],
         [0.16699984,0.212740576,0.216777868,0.219333118,0.217690489,0.208681783,0.217690489,0.219333118,0.216777868,0.212740576,0.16699984],
         [0.188564579,0.220684598,0.222088376,0.221227419,0.218695369,0.211008781,0.218695369,0.221227419,0.222088376,0.220684598,0.188564579]]

    a1=[0.1,1.0,5.0,13.51]
    a2=[0.22,1.0,5.0,13.0]
    b1=[-0.999,-0.9,-0.8,-0.5,-0.3,0.0,0.3,0.5,0.8,0.9,0.999]
    b2=[-0.99,-0.9,-0.8,-0.5,-0.3,0.0,0.3,0.5,0.8,0.9,0.99]
    c0=[0.014,0.1]
    mean0=mean
    if mean < m0[0][10]:
        mean0=m0[0][10]
    elif mean > m0[0][0]:
        mean0=m0[0][0]
    f1 = interpolate.interp1d(m0[0],s0[0])
    w1 = f1(mean0)
    g1 = interpolate.interp1d(m0[0],b1)
    v1 = g1(mean0)
    if mean < m0[1][10]:
        mean0=m0[1][10]
    elif mean > m0[1][0]:
        mean0=m0[1][0]
    f2 = interpolate.interp1d(m0[1],s0[1])
    w2 = f2(mean0)
    g2 = interpolate.interp1d(m0[1],b1)
    v2 = g2(mean0)
    if mean < m0[2][10]:
        mean0=m0[2][10]
    elif mean > m0[2][0]:
        mean0=m0[2][0]
    f3 = interpolate.interp1d(m0[2],s0[2])
    w3 = f3(mean0)
    g3 = interpolate.interp1d(m0[2],b1)
    v3 = g3(mean0)
    if mean < m0[3][10]:
        mean0=m0[3][10]
    elif mean > m0[3][0]:
        mean0=m0[3][0]
    f4 = interpolate.interp1d(m0[3],s0[3])
    w4 = f4(mean0)
    g4 = interpolate.interp1d(m0[3],b1)
    v4 = g4(mean0)
    if mean < m0[4][10]:
        mean0=m0[4][10]
    elif mean > m0[4][0]:
        mean0=m0[4][0]
    f5 = interpolate.interp1d(m0[4],s0[4])
    w5 = f5(mean0)
    g5 = interpolate.interp1d(m0[4],b2)
    v5 = g5(mean0)
    if mean < m0[5][10]:
        mean0=m0[5][10]
    elif mean > m0[5][0]:
        mean0=m0[5][0]
    f6 = interpolate.interp1d(m0[5],s0[5])
    w6 = f6(mean0)
    g6 = interpolate.interp1d(m0[5],b2)
    v6 = g6(mean0)
    if mean < m0[6][10]:
        mean0=m0[6][10]
    elif mean > m0[6][0]:
        mean0=m0[6][0]
    f7 = interpolate.interp1d(m0[6],s0[6])
    w7 = f7(mean0)
    g7 = interpolate.interp1d(m0[6],b2)
    v7 = g7(mean0)
    if mean < m0[7][10]:
        mean0=m0[7][10]
    elif mean > m0[7][0]:
        mean0=m0[7][0]
    f8 = interpolate.interp1d(m0[7],s0[7])
    w8 = f8(mean0)
    g8 = interpolate.interp1d(m0[7],b2)
    v8 = g8(mean0)

    if s <= w4:
        if s < w1:
            a=a1[0]
            b=v1
        else:
            f0 = interpolate.interp1d([w1,w2,w3,w4],a1)
            a = f0(s)
            g0 = interpolate.interp1d([w1,w2,w3,w4],[v1,v2,v3,v4])
            b = g0(s)
        c=c0[0]

    else:
        if s > w8:
            a=a2[2]
            b=v8
        else:
            f0 = interpolate.interp1d([w5,w6,w7,w8],a2)
            a = f0(s)
            g0 = interpolate.interp1d([w5,w6,w7,w8],[v5,v6,v7,v8])
            b = g0(s)
        c=c0[1]
    bx=a*b

    x0=np.linspace(norminvgauss.ppf(c, a, bx),norminvgauss.ppf(1-c, a, bx), div)
    y0=norminvgauss.pdf(x0,a,bx)
    if y0[0]<y0[div-1]:
        dy=y0[0]
    else:
        dy=y0[div-1]
        #
    area=0
    i0=0
    yx=0
    y1=[]
    x=[]

    x=np.linspace(xmin,xmax,len(x0))
    for i in range(len(x0)):
        y1.append(y0[i])
        area += y1[i]*dx
    y=y1/np.array([area]*len(x0))

    return x,y


def file_data(infilename,outdataname,skipline,colum=[[0,1],[5,6],[10,11]],judg=[1,2], limitcol=[0,0], limitvalue=[-1000,1000],cntmx=0):

    if type(infilename) is list:
        if len(infilename)==2:
            dirname=[]
            infile=infilename[0]
            sht=infilename[1]
            if 'csv' in sht:
                dirname=infilename[0]
                infile=infilename[1]
        elif len(infilename)==3:
            dirname=infilename[0]
            infile=infilename[1]
            sht=infilename[2]
            
    else:
        infile=infilename
        dirname=[]
    filename=os.path.splitext(infile)
    outfilename=os.path.splitext(outdataname)
    if filename[1]=='.csv':
        if dirname==[]:
            df = pd.read_csv(infile,skiprows=skipline,na_values=999)
            dt = df.values
        else:
            dt=[]
            allFiles = glob.glob(dirname+infile)
            for file_ in allFiles:
                df = pd.read_csv(file_,skiprows=skipline,na_values=999)
                dt.extend(df.values)
            
    elif filename[1]=='.xlsx':
        if dirname==[]:
            df = pd.read_excel(infile,skiprows=skipline,sheet_name=sht,na_values=999)
            dt = df.values
        else:
            dt=[]
            allFiles = glob.glob(dirname+infile)
            for file_ in allFiles:
                df = pd.read_excel(file_,skiprows=skipline,sheet_name=sht,na_values=999)
                dt.extend(df.values)

    data=dt
    pu_data=[]
    num_data=len(colum)
    tmp_data=[]
    lastcol1=[]
    lastcol2=[]
    cnt=[]
    flagy=[]
    cntflag=[]
    fstvale=[]
    tmp=[]
    for i in range(num_data):
        tmp_data.append([])
        lastcol1.append([0])
        lastcol2.append([0])
        cnt.append([0])
        flagy.append([0])
        cntflag.append([0])
        fstvale.append([])
        tmp.append([])
    ncase=6
    case=[]
    vcase=[]
    for i in range(ncase):
        case.append([])
        vcase.append([])
    njudg=len(judg)
    flag1=[]
    flag2=[]
    for i in range(njudg):
        case[judg[i]].append(limitcol[i])
        vcase[judg[i]].append(limitvalue[i])
        if i<=3:
            flag1.append([])
        else:
            flag2.append([])

    for i in range(len(data)):
        xx=data[i]
                       
        for j in range(len(colum)):
            for k in range(len(flag1)):
                flag1[k]=[]
            
            flagx=[]
            for k in range(len(case[0])):
                if xx[colum[j][case[0][k]]]>=vcase[0][k]:
                    
                    flagx.append(1)
                else:
                    flagx.append(0)
            if len(case[0])==0:
                flag1[0].append(1)
            elif all(n==1 for n in flagx):
                
                flag1[0].append(1)
            else:
                flag1[0].append(0)
                
            flagx=[]
            for k in range(len(case[1])):
                if xx[colum[j][case[1][k]]]<=vcase[1][k]:
                    flagx.append(1)
                else:
                    flagx.append(0)
            if len(case[1])==0:
                flag1[1].append(1)
            elif all(n==1 for n in flagx):
                flag1[1].append(1)
            else:
                flag1[1].append(0)
            
            if len(case[2])>=1 and i>=2:
                if (xx[colum[j][case[2][0]]]-lastcol1[j]>=vcase[2][0]):
                    flag1[2].append(1)
                else:
                    flag1[2].append(0)
            else:
                flag1[2].append(1)

            if len(case[3])>=1 and i>=2:
                if (lastcol2[j]-xx[colum[j][case[3][0]]]>=vcase[3][0]):
                    flag1[3].append(1)
                else:
                    flag1[3].append(0)
            else:
                flag1[3].append(1)
            
            if all(n[0]==1 for n in flag1):
                flagy[j][0]=1
            else:
                flagy[j][0]=0
            
            if len(case[4])>=1:
                if xx[colum[j][case[4][0]]]>=vcase[4][0] and flagy[j][0]==1:
                    cntflag[j][0]=1
                    cnt[j][0]+=1
                else:
                    cntflag[j][0]=0
                if xx[colum[j][case[4][0]]]>=vcase[4][0] and (flagy[j][0]==1 or cnt[j][0]>0):
                    cnt[j][0]+=1
                else:
                    cnt[j][0]=0
                                                                          
                if cnt[j][0]>=cntmx[0]:
                    cntflag[j][0]=2
                    cnt[j][0]=0
            else:
                if flagy[j][0]==1:
                    cntflag[j][0]=-1
                else:
                    cntflag[j][0]=0
                    
            flag6=0

            if cntflag[j][0]==1:
                for k in range(len(colum[j])):
                    tmp[j].append(xx[colum[j][k]])
            if cntflag[j][0]==2:
                fstvale[j].extend(tmp[j])
                tmp[j]=[]
                for k in range(len(colum[j])):
                    fstvale[j].append(xx[colum[j][k]])
                
            if cntflag[j][0]==-1:
                for k in range(len(colum[j])):
                    fstvale[j].append(xx[colum[j][k]])
            if fstvale[j]!=[]:
                pu_data.append(fstvale[j])
                fstvale[j]=[]
            if case[2]!=[]:
                lastcol1[j]=xx[colum[j][case[2][0]]]
            if case[3]!=[]:
                lastcol2[j]=xx[colum[j][case[3][0]]]
    data_np=np.array(pu_data)
    dim=len(colum)-1
    df=pd.DataFrame(pu_data)
    df.to_csv(outdataname)

def predistrb(tmin,tmax,div,dim):

    dt=[]
    if dim == 1:
        divt=div[0]
        t=np.linspace(tmin[0],tmax[0],divt)#11/10
        t1=t
        dt.append((tmax[0]-tmin[0])/divt)
        t0=t1
        p1=np.zeros(divt)
        p0=p1
    elif dim ==2:
        divtt=div[0]*div[1]
        t=[np.linspace(tmin[0],tmax[0],div[0]),np.linspace(tmin[1],tmax[1],div[1])]
        tx,ty=np.meshgrid(t[0],t[1])#
        dt.append((tmax[0]-tmin[0])/div[0])
        dt.append((tmax[1]-tmin[1])/div[1])
        p1=np.zeros(divtt)
        p0=np.zeros((div[1],div[0]))
        t1=[tx.flatten(),ty.flatten()]
        t0=[tx,ty]
    elif dim==3:
        divtt=div[0]*div[1]*div[2]
        t=[np.linspace(tmin[0],tmax[0],div[0]),np.linspace(tmin[1],tmax[1],div[1]),np.linspace(tmin[2],tmax[2],div[2])]
        tx,ty,tz=np.meshgrid(t[0],t[1],t[2])#
        dt.append((tmax[0]-tmin[0])/div[0])
        dt.append((tmax[1]-tmin[1])/div[1])
        dt.append((tmax[2]-tmin[2])/div[2])
        p1=np.zeros(divtt)
        p0=np.zeros((div[2],div[1],div[0]))
        t1=[tx.flatten(),ty.flatten(),tz.flatten()]
        t0=[tx,ty,tz]
    return t,t0,p0,t1,p1,dt


def file2cor(infilename,skipline,colum1,colum2):

    if type(infilename) is list:
        if len(infilename)==2:
            dirname=[]
            infile=infilename[0]
            sht=infilename[1]
            if 'csv' in sht:
                dirname=infilename[0]
                infile=infilename[1]
        elif len(infilename)==3:
            dirname=infilename[0]
            infile=infilename[1]
            sht=infilename[2]
            
    else:
        infile=infilename
        dirname=[]
    filename=os.path.splitext(infile)
    if filename[1]=='.csv':
        if dirname==[]:
            df = pd.read_csv(infile,usecols=[colum1,colum2],skiprows=skipline,na_values=999)
            dt = df.values
        else:
            dt=[]
            allFiles = glob.glob(dirname+infile)
            for file_ in allFiles:
                df = pd.read_csv(file_,usecols=[colum1,colum2],skiprows=skipline,na_values=999)
                dt.extend(df.values)
        
    elif filename[1]=='.xlsx':
        if dirname==[]:
            df = pd.read_excel(infile,usecols=[colum1,colum2],skiprows=skipline,sheet_name=sht,na_values=999)
            dt = df.values
        else:
            dt=[]
            allFiles = glob.glob(dirname+infile)
            for file_ in allFiles:
                df = pd.read_excel(file_,usecols=[colum1,colum2],skiprows=skipline,sheet_name=sht,na_values=999)
                dt.extend(df.values)

    data=np.array(dt)
    x=data[:,0]
    y=data[:,1]
    meanx=sum(x)/len(x)
    meany=sum(y)/len(y)
    disco=0
    disx=0
    disy=0
    for i in range(len(data)):
        dx=x[i]-meanx
        dy=y[i]-meany
        disco+=dx*dy
        disx+=dx*dx
        disy+=dy*dy
    corel=disco/(disx**0.5*disy**0.5)
    return corel

def cor2ran(x,xmin,xmax,y,ymin,ymax,r=0):
    #
    if y>=ymin and y<=ymax and x>=xmin and x<=xmax:
        kx=(x-xmin)/(xmax-xmin)
        ky=(y-ymin)/(ymax-ymin)
        
        if r==0:
            flag=1
        elif r<0:
            r1=-1*r
            if kx<=ky and kx+ky<=1:
                if kx==0:
                    r0=0
                else:
                    r0=2/((1-ky)/kx+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
            elif kx>ky and kx+ky<=1:
                #
                if ky==0:
                    r0=0
                else:
                    r0=2/((1-kx)/ky+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
            elif kx>ky and kx+ky>1:
                #
                if kx==1:
                    r0=0
                else:
                    r0=2/(ky/(1-kx)+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
            elif kx<=ky and kx+ky>1:
                #
                if ky==1:
                    r0=0
                else:
                    r0=2/(kx/(1-ky)+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
        elif r>0:
            r1=r
            if kx<=ky and kx+ky<=1:
                #
                if kx==0:
                    r0=0
                else:
                    r0=2/(ky/kx+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
            elif kx>ky and kx+ky<=1:
                #
                if ky==0:
                    r0=0
                else:
                    r0=2/(kx/ky+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
            elif kx>ky and kx+ky>1:
                #
                if kx==1:
                    r0=0
                else:
                    r0=2/((1-ky)/(1-kx)+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
            elif kx<=ky and kx+ky>1:
                if ky==1:
                    r0=0
                else:
                    r0=2/((1-kx)/(1-ky)+1)
                if r1<=r0:
                    flag=1
                else:
                    flag=0
    else:
        flag=0
    #
    return flag


def slice2d_bunpu(ti,t1,tmin,tmax,dt,divs,xmin,xmax,ymin,ymax):
    plmax=9
    plmin=9
    if xmin[0]>=0 and xmin[1]>=0:
        if xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            plmax=0
        elif xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmax=1
        if xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            plmin=3
        elif xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmin[1]
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            plmin=2
    elif xmin[0]<0 and xmax[0]>=0 and xmin[1]>0:#
        if xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti]:#
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            plmax=0
        elif xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmax=1
        elif xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti]:#
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[0][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            plmax=3
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmin[1]
            min_tx=tmin[0]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            plmin=2
    elif xmax[0]<=0 and xmin[1]>=0:###原点が右下
        if xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_tx=tmin[0]*t1[1][ti]/t1[0][ti]
            max_ty=tmin[0]
            plmax=3
        elif xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]
            max_sy=xmax[1]
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            plmax=0
        if xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]>=xmin[1]*t1[0][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmin[1]
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            plmin=2
        elif xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmin=1
    elif xmin[0]>=0 and xmin[1]<=0 and xmax[1]>=0:
        if xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti]:#
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            plmax=0
        elif xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmax=1
        elif xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            max_sy=xmin[1]
            max_tx=tmin[1]
            max_ty=tmin[1]*t1[0][ti]/t1[1][ti]
            plmax=2
        if xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            plmin=3
    elif xmin[0]<0 and xmax[0]>0 and xmin[1]<0 and xmax[1]>0:
        if xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti] and t1[1][ti]>=0:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            plmax=0
        elif xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti] and t1[0][ti]>=0:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[0][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti] and t1[1][ti]<0:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            max_sy=xmin[1]
            max_tx=tmin[1]
            max_ty=tmin[1]*t1[0][ti]/t1[1][ti]
            plmax=2
        elif xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti] and t1[0][ti]<0:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[0][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            plmax=3
        if xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti] and t1[1][ti]<0:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmax[1]
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            plmin=0
        elif xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti] and t1[0][ti]<0:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmin=1
        elif xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti] and t1[1][ti]>=0:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmin[1]
            min_tx=tmin[0]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            plmin=2
        elif xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti] and t1[0][ti]>=0:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            plmin=3
    elif xmax[0]<=0 and xmin[1]<=0 and xmax[1]>=0:
        if xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_tx=tmin[0]*t1[1][ti]/t1[0][ti]
            max_ty=tmin[0]
            plmax=3
        elif xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]
            max_sy=xmax[1]
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            plmax=0
        elif xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmin[1]
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            plmax=2
        if xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmin=1
    elif xmin[0]>=0 and xmax[1]<=0:
        if xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            max_sy=xmin[1]
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            plmax=2
        if xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            plmin=3
        elif xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmax[1]
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            plmin=0
    elif xmin[0]<0 and xmax[0]>=0 and xmax[1]<0:
        if xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            plmax=2
        elif xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmax=1
        elif xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]
            max_tx=tmax[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            plmax=3
        if xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmax[1]
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            plmin=0
    elif xmax[0]<=0 and xmax[1]<=0:
        if xmin[1]<=xmin[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_tx=tmin[0]*t1[1][ti]/t1[0][ti]
            max_ty=tmin[0]
            plmax=3
        elif xmin[0]<=xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]<=xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]
            max_sy=xmin[1]
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            plmax=2
        if xmin[0]<=xmax[1]*t1[0][ti]/t1[1][ti] and xmax[0]>=xmax[1]*t1[0][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]
            min_sy=xmax[1]
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            plmin=0
        elif xmin[1]<=xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>=xmax[0]*t1[1][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            plmin=1
    if plmax!=9 and plmin!=9:
        max_s=(max_sx**2+max_sy**2)**0.5
        max_t=(max_tx**2+max_ty**2)**0.5
        min_s=(min_sx**2+min_sy**2)**0.5
        min_t=(min_tx**2+min_ty**2)**0.5
        dt_t=(dt[0]**2+dt[1]**2)**0.5
        divst=-1*((-1*(max_t-min_t))//dt_t)
    else:
        min_sx=0
        max_sx=0
        min_sy=0
        max_sy=0
        max_s=0
        max_t=0
        min_s=0
        min_t=0
        divst=1
    return divst, min_t, max_t, min_s, max_s, min_sx, max_sx, min_sy, max_sy, plmin, plmax



def slice3d_bunpu(ti,t1,tmin,tmax,dt,divs,xmin,xmax,ymin,ymax):
    plmax=9
    plmin=9

    if xmin[0]>0 and xmin[1]>0 and xmin[2]>0:
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:#
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmin[0]>0 and xmin[1]<0 and xmax[1]>0 and xmin[2]>0:
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:#
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmin[0]>0 and xmax[1]<0 and xmin[2]>0:

        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmin[0]>0 and xmin[1]>0 and xmin[2]<0 and xmax[2]>0:
        #plmax0145
        #plmin23
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:#
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
        
    elif xmin[0]>0 and xmin[1]<0 and xmax[1]>0 and xmin[2]<0 and xmax[2]>0:
        #plmax01245
        #plmin3
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:#
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:#
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3

    elif xmin[0]>0 and xmin[1]>0 and xmax[2]<0:
        #plmax1245
        #plmin03
        #
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
    elif xmin[0]>0 and xmin[1]<0 and xmax[1]>0 and xmax[2]<0:
        #plmax015
        #plmin234
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    elif xmin[0]>0 and xmax[1]<0 and xmax[2]<0:
        #plmax125
        #plmin034
        #
        #
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    elif xmin[0]<0 and xmax[0]>0 and xmin[1]>0 and xmin[2]>0:
        #plmax0134
        #plmin25
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmin[0]<0 and xmax[0]>0 and xmin[1]<0 and xmax[1]>0 and xmin[2]>0:
        #plmax01234
        #plmin5
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmin[0]<0 and xmax[0]>0 and xmax[1]<0 and xmin[2]>0:
        #plmax1234
        #plmin05
        #
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmin[0]<0 and xmax[0]>0 and xmin[1]>0 and xmin[2]<0 and xmin[2]>0:
        #plmax01345
        #plmin2
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
    elif xmin[0]<=0 and xmax[0]>=0 and xmin[1]<=0 and xmax[1]>=0 and xmin[2]<=0 and xmax[2]>=0:
        #原点が範囲内にある場合
        #plmax012345
        #plmin012345
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and ti[1][ti]>=0:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and ti[0][ti]>=0:#
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and ti[1][ti]<0:#
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and ti[2][ti]>=0:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti] and ti[2][ti]<0:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and ti[1][ti]<0:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and ti[0][ti]<0:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and ti[1][ti]>=0:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti] and xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and ti[0][ti]>=0:
            min_sx=xmin[0]
            min_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmin[0]
            min_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmin=3
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and ti[2][ti]<0:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti] and ti[2][ti]>=0:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
        
    elif xmin[0]<=0 and xmax[0]>=0 and xmax[1]<0 and xmin[2]<=0 and xmax[2]>0:
        #plmax12345
        #plmin0
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        #
    elif xmin[0]<=0 and xmax[0]>=0 and xmin[1]>0 and xmax[2]<0:
        #plmax0135
        #plmin24
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    elif xmin[0]<=0 and xmax[0]>=0 and xmin[1]<0 and xmax[1]>0 and xmax[2]<0:
        #plmax01235
        #plmin4
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    elif xmin[0]<=0 and xmax[0]>=0 and xmax[1]<0 and xmax[2]<0:
        #plmax1235
        #plmin04
        #
        #
        if xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmax[0]
            max_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmax[0]
            max_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmax=1
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    elif xmax[0]<0 and xmin[1]>0 and xmin[2]>0:
        #plmax034
        #plmin125
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        elif xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmax[0]<0 and xmin[1]<0 and xmax[1]>0 and xmin[2]>0:
        #plmax0234
        #plmin15
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmax[0]<0 and xmax[1]<0 and xmin[2]>0:
        #plmax234
        #plmin015
        #
        if xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        if xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:
            min_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmin[2]
            min_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmin[2]
            plmin=5
    elif xmax[0]<0 and xmin[1]>0 and xmin[2]<0 and xmax[2]>0:
        #plmax0345
        #plmin12
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        elif xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
    elif xmax[0]<0 and xmin[1]<0 and xmax[1]>0 and xmin[2]<0 and xmax>0:
        #plmax02345
        #plmin1
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
    elif xmax[0]<0 and xmax[1]<0 and xmin[2]<0 and xmax>0:
        #plmax2345
        #plmin01
        #
        if xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti] and xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti]:
            max_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmax[2]
            max_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmax[2]
            plmax=4
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
    elif xmax[0]<0 and xmin[1]>0 and xmax<0:
        #plmax035
        #plmin124
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        elif xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[0]>xin[1]*t1[0][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmin[1]
            min_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmin[1]
            min_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmin=2
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    elif xmax[0]<0 and xmin[1]<0 and xmax[1]>0 and xmax<0:
        #plmax0235
        #plmin14
        #
        if xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti] and xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti]:
            max_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmax[1]
            max_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmax[1]
            max_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmax=0
        elif xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        elif xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    elif xmax[0]<0 and xmax[1]<0 and xmax<0:
        #plmax235
        #plmin014
        if xmax[0]>xmin[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmin[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmin[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmin[1]*t1[2][ti]/t1[1][ti]:
            max_sx=xmin[1]*t1[0][ti]/t1[1][ti]#
            max_sy=xmin[1]
            max_sz=xmin[1]*t1[2][ti]/t1[1][ti]#
            max_tx=tmin[1]*t1[0][ti]/t1[1][ti]
            max_ty=tmin[1]
            max_tz=tmin[1]*t1[2][ti]/t1[1][ti]
            plmax=2
        elif xmax[2]>xmin[0]*t1[2][ti]/t1[0][ti] and xmin[2]<xmin[0]*t1[2][ti]/t1[0][ti] and xmax[1]>xmin[0]*t1[1][ti]/t1[0][ti] and xmin[1]<xmin[0]*t1[1][ti]/t1[0][ti]:
            max_sx=xmin[0]
            max_sy=xmin[0]*t1[1][ti]/t1[0][ti]#
            max_sz=xmin[0]*t1[2][ti]/t1[0][ti]#
            max_tx=tmin[0]
            max_ty=tmin[0]*t1[1][ti]/t1[0][ti]
            max_tz=tmin[0]*t1[2][ti]/t1[0][ti]
            plmax=3
        elif xmax[0]>xmin[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmin[2]*t1[0][ti]/t1[2][ti] and xmax[1]>xmin[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmin[2]*t1[1][ti]/t1[2][ti]:#
            max_sx=xmin[2]*t1[0][ti]/t1[2][ti]#
            max_sy=xmin[2]*t1[1][ti]/t1[2][ti]#
            max_sz=xmin[2]
            max_tx=tmin[2]*t1[0][ti]/t1[2][ti]
            max_ty=tmin[2]*t1[1][ti]/t1[2][ti]
            max_tz=tmin[2]
            plmax=5
        if xmax[0]>xmax[1]*t1[0][ti]/t1[1][ti] and xmin[0]<xmax[1]*t1[0][ti]/t1[1][ti] and xmax[2]>xmax[1]*t1[2][ti]/t1[1][ti] and xmin[2]<xmax[1]*t1[2][ti]/t1[1][ti]:
            min_sx=xmax[1]*t1[0][ti]/t1[1][ti]#
            min_sy=xmax[1]
            min_sz=xmax[1]*t1[2][ti]/t1[1][ti]#
            min_tx=tmax[1]*t1[0][ti]/t1[1][ti]
            min_ty=tmax[1]
            min_tz=tmax[1]*t1[2][ti]/t1[1][ti]
            plmin=0
        elif xmin[1]<xmax[0]*t1[1][ti]/t1[0][ti] and xmax[1]>xmax[0]*t1[1][ti]/t1[0][ti] and xmin[2]<xmax[0]*t1[2][ti]/t1[0][ti] and xmax[2]>xmax[0]*t1[2][ti]/t1[0][ti]:
            min_sx=xmax[0]
            min_sy=xmax[0]*t1[1][ti]/t1[0][ti]#
            min_sz=xmax[0]*t1[2][ti]/t1[0][ti]#
            min_tx=tmax[0]
            min_ty=tmax[0]*t1[1][ti]/t1[0][ti]
            min_tz=tmax[0]*t1[2][ti]/t1[0][ti]
            plmin=1
        elif xmax[1]>xmax[2]*t1[1][ti]/t1[2][ti] and xmin[1]<xmax[2]*t1[1][ti]/t1[2][ti] and xmax[0]>xmax[2]*t1[0][ti]/t1[2][ti] and xmin[0]<xmax[2]*t1[0][ti]/t1[2][ti]:
            min_sx=xmax[2]*t1[0][ti]/t1[2][ti]#
            min_sy=xmax[2]*t1[1][ti]/t1[2][ti]#
            min_sz=xmax[2]
            min_tx=tmax[2]*t1[0][ti]/t1[2][ti]
            min_ty=tmax[2]*t1[1][ti]/t1[2][ti]
            min_tz=tmax[2]
            plmin=4
    if plmax!=9 and plmin!=9:
        max_s=(max_sx**2+max_sy**2+max_sz**2)**0.5
        max_t=(max_tx**2+max_ty**2+max_tz**2)**0.5
        min_s=(min_sx**2+min_sy**2+min_sz**2)**0.5
        min_t=(min_tx**2+min_ty**2+min_tz**2)**0.5
        dt_t=(dt[0]**2+dt[1]**2+dt[2]**2)**0.5
        divst=-1*((-1*(max_t-min_t))//dt_t)
    if plmax==9 or plmin==9:# or abs(max_s-min_s)<min()
        max_sx=0
        max_sy=0
        max_sz=0
        max_s=0
        max_t=0
        min_sx=0
        min_sy=0
        min_sz=0
        min_s=0
        min_t=0
        divst=1
    return divst, min_t, max_t, min_s, max_s, min_sx, max_sx, min_sy, max_sy, min_sz, max_sz, plmin, plmax






class bunpu(object):
    """
    確率分布の定義:bunpu()
    """
    def __init__(self):
        
        d=0
        x=[]
        r=0
        o=0
        tgt=''
        nm=''
        self.name=[]
        self.dim=d
        self.para=[]
        self.mesh=[]
        self.flatten=[]
        self.xmin=[]
        self.xmax=[]
        self.xmean=[]
        self.sdev=[]
        self.skew=[]
        self.div=[]
        self.dx=[]
        self.hist=0
        self.kh=0
        self.pmax=1
        self.cmap=0

    def __enter__(self):
        print('START')
        return self

    def __add__(self,other):

        divt=self.div

        WW=bunpu()
        WW=self.bunpu_add(other,divt)
        return WW
    
    def __sub__(self,other):

        divt=self.div

        WW=bunpu()
        WW=self.bunpu_sub(other,divt)
        return WW
        
    def __mul__(self,other):

        divt=self.div
        WW=bunpu()
        WW=self.bunpu_product(other,divt)
        return WW

    def __truediv__(self,other):

        divt=self.div

        WW=bunpu()
        WW=self.bunpu_division(other,divt)
        return WW
    
    def bunpu_data(self,infilename,outdataname,skipline,colum,divn=100, limit=0,kh=0, fontsize=24):
        fp = FontProperties(fname=r'C:\Windows\Fonts\msgothic.ttc',size=fontsize)
        plt.rcParams["font.size"] = fontsize
        if kh!=0:
            bw=kh
        if type(infilename) is list:
            if len(infilename)==2:
                dirname=[]
                infile=infilename[0]
                sht=infilename[1]
                if 'csv' or 'txt' in sht:
                    dirname=infilename[0]
                    infile=infilename[1]
            elif len(infilename)==3:
                dirname=infilename[0]
                infile=infilename[1]
                sht=infilename[2]
                
        else:
            infile=infilename
            dirname=[]
        filename=os.path.splitext(infile)
        outfilename=os.path.splitext(outdataname)
        fg=0
        if filename[1]=='.csv':
            if dirname==[]:
                df = pd.read_csv(infile,usecols=colum,skiprows=skipline,na_values=999,encoding='cp932')
                
                dt = df.values
                    
            else:
                dt=[]
                allFiles = glob.glob(dirname+infile)
                for file_ in allFiles:
                    df = pd.read_csv(file_,usecols=colum,skiprows=skipline,na_values=999,encoding='cp932')
                    dt.extend(df.values)
            
        elif filename[1]=='.xlsx':
            if dirname==[]:
                df = pd.read_excel(infile,usecols=colum,skiprows=skipline,sheet_name=sht,na_values=999)
                dt = df.values
            else:
                dt=[]
                allFiles = glob.glob(dirname+infile)
                for file_ in allFiles:
                    df = pd.read_excel(file_,usecols=colum,skiprows=skipline,sheet_name=sht,na_values=999)
                    dt.extend(df.values)
        elif filename[1]=='.txt':
            if dirname==[]:
                df = pd.read_table(infile,usecols=colum,skiprows=skipline,na_values=999,encoding='cp932')
                dt = df.values
            else:
                fg=1
                dt=[]
                allFiles = glob.glob(dirname+infile)
                for file_ in allFiles:
                    df = pd.read_table(file_,usecols=colum,skiprows=skipline,na_values=999,encoding='cp932')
                    dt.extend(df.values)
        data=dt
        pu_data=[]
        flag=0
        if limit == 0:
            flag=0
        else:
            flag=1
            limitcol=limit[0]
            lowlimit=limit[1]
            uplimit=limit[2]
        col=np.array(colum)
        tmp=col.argsort()
        for i in range(len(data)):
            xx=[]
            for j in range(len(tmp)):
                xx.append(data[i][tmp[j]])
            if flag==0:
                pu_data.append(xx)
            elif xx[limitcol]<=uplimit and xx[limitcol]>=lowlimit:# or np.isnan(xx[limitcol]):
                pu_data.append(xx)
        data_np=np.array(pu_data)
        #
        dim=len(colum)
        if dim==1:
            tmin = [np.min(data_np)]
            tmax = [np.max(data_np)]
            if kh==0:
                bw=2*(tmax[0]-tmin[0])/divn[0]
            x_grid = np.linspace(tmin[0], tmax[0], divn[0])
            wt=0
            if wt==0.0:
                weights = np.ones_like(data_np)/float(len(data_np))
            else:
                weights = np.ones_like(data_np)/(wt*float(len(data_np)))
            #
            kde_mode_bw = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_np)#
            score_bw = kde_mode_bw.score_samples(x_grid[:, np.newaxis])

            fg=plt.figure(figsize=(14,7))
            ax1 = fg.add_subplot(111)
            
            ln1=ax1.plot(x_grid, np.exp(score_bw), label="origin")
            ax2 = ax1.twinx()
            ln2 = ax2.hist(data_np, alpha=0.3, bins=divn[0], weights=weights)
            
            h1, l1 = ax1.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax1.legend(h1+h2, l1+l2, loc='lower right')
            
            ax1.set_xlabel(u''+'x軸', fontproperties=fp)
            ax1.set_ylabel(r'kakuritsubunpu')
            ax1.grid(True)
            ax2.set_ylabel(r'hindo')
            
            
            plt.title(u''+'分布', fontproperties=fp) 
            plt.xlabel(u''+'x軸', fontproperties=fp) 
            outgraphname=outfilename[0]+'.png'
            plt.savefig(outgraphname)
            
            y=[x_grid,np.exp(score_bw)]
            adarray=pd.DataFrame(y)
            outdata=adarray.T
            outfile=outfilename[0]+'.csv'
            outdata.to_csv(outfile)
            self.pmax=max(y[1])
            self.div=divn
            self.para=y
            self.mesh=y
            self.flatten=y
            self.dim=1
            self.xmax=tmax
            self.xmin=tmin
            self.dx=[(tmax[0]-tmin[0])/(divn[0]-1)]
        elif dim==2:
            tmin=[data_np[:,0].min(),data_np[:,1].min()]
            tmax=[data_np[:,0].max(),data_np[:,1].max()]
            if kh==0:
                bw0=[]
                for i in range(dim):
                    bw0.append(2*(tmax[i]-tmin[i])/divn[i])
                bw=(bw0[0]+bw0[1])/2

            fig = plt.figure(figsize=(14,7))
            ax = fig.gca()
            ax.set_xlim(tmin[0],tmax[0])
            ax.set_ylim(tmin[1],tmax[1])

            H = ax.hist2d(data_np[:,0], data_np[:,1], bins=divn, cmap=cm.gray)
            fig.colorbar(H[3],ax=ax)
            outgraphname=outfilename[0]+'.png'
            plt.savefig(outgraphname)

            X = np.vstack((data_np[:,0], data_np[:,1])).T

            kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(X)
            np.exp(kde.score_samples(X))

            f0 = np.linspace(tmin[0],tmax[0], divn[0])
            f1 = np.linspace(tmin[1],tmax[1], divn[1])
            xx, yy = np.meshgrid(f0, f1)
            positions = np.vstack([xx.ravel(), yy.ravel()]).T
            scores = np.exp(kde.score_samples(positions)).T
            f = np.reshape(scores, xx.shape)
            gr=1
            if gr==1:
                fig = plt.figure(figsize=(14,7))
                ax = fig.gca()
                ax.set_xlim(tmin[0],tmax[0])
                ax.set_ylim(tmin[1],tmax[1])
                if fontsize>15:
                    cset = ax.contour(xx, yy, f, linewidths=4, colors='r')
                else:
                    cset = ax.contour(xx, yy, f, colors='r')
                ax.clabel(cset, inline=1, fontsize=fontsize)
                
                H = ax.hist2d(data_np[:,0], data_np[:,1], bins=divn, cmap=cm.gray)
                fig.colorbar(H[3],ax=ax)
                outgraphname=outfilename[0]+'2.png'
                plt.savefig(outgraphname)
                plt.close()

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(xx, yy, f,cmap='viridis',linewidth=0)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            outgraphname=outfilename[0]+'.png'
            plt.savefig(outgraphname)
            plt.close()
            self.pmax=(f.flatten()).max()
            self.div=divn
            self.para=[f0,f1]
            self.mesh=[xx,yy,f]
            self.flatten=[xx.flatten(),yy.flatten(),f.flatten()]
            self.dim=2
            self.xmin=[tmin[0],tmin[1]]
            self.xmax=[tmax[0],tmax[1]]
            self.dx=[(tmax[0]-tmin[0])/(divn[0]-1),(tmax[1]-tmin[1])/(divn[1]-1)]
        elif dim==3:
            tmin=[data_np[:,0].min(),data_np[:,1].min(),data_np[:,2].min()]
            tmax=[data_np[:,0].max(),data_np[:,1].max(),data_np[:,2].max()]
            if kh==0:
                bw0=[]
                for i in range(dim):
                    bw0.append(2*(tmax[i]-tmin[i])/divn[i])
                bw=(bw0[0]+bw0[1]+bw0[2])/3

            fig = plt.figure(figsize=(14,7))
            ax = fig.add_subplot(111, projection='3d')
            plt.figaspect(1)
            cmap = cm.binary
            cmap_data = cmap(np.arange(cmap.N))
            lnmap=len(cmap_data)
            for i in range(lnmap):
                trn=i/lnmap
                cmap_data[i, 3] = trn 
            customized_binary = colors.ListedColormap(cmap_data)
            gr=3
            if gr==0 or gr==2:
                density,edges=np.histogramdd(data_np,bins=divn)
                f0=[(edges[0][i]+edges[0][i+1])/2 for i in range(len(edges[0])-1)]
                f1=[(edges[1][i]+edges[1][i+1])/2 for i in range(len(edges[1])-1)]
                f2=[(edges[2][i]+edges[2][i+1])/2 for i in range(len(edges[2])-1)]
                edgx,edgy,edgz=np.meshgrid(f0, f1, f2)
                sc = ax.scatter(edgx, edgy, edgz, c=density.ravel(), cmap=customized_binary)
                fig.colorbar(sc)
            
            elif gr==1 or gr==3:
                fig = plt.figure(figsize=(14,7))
                ax = fig.add_subplot(111, projection='3d')
                plt.figaspect(1)
                X = data_np.T
                kde = gaussian_kde(X)
                f0 = np.linspace(tmin[0],tmax[0], divn[0])
                f1 = np.linspace(tmin[1],tmax[1], divn[1])
                f2 = np.linspace(tmin[2],tmax[2], divn[2])
                xx, yy, zz = np.mgrid[tmin[0]:tmax[0]:(tmax[0]-tmin[0])/divn[0],tmin[1]:tmax[1]:(tmax[1]-tmin[1])/divn[1],\
                                      tmin[2]:tmax[2]:(tmax[2]-tmin[2])/divn[2]]
                positions = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()])
                
                scores = kde(positions)

                density = np.reshape(scores, xx.shape)
                sc = ax.scatter(xx, yy, zz, c=density.ravel(), cmap=customized_binary)
                fig.colorbar(sc)
            if gr==2 or gr==3:
                ax.set_xlim((tmin[0],tmax[0]))
                ax.set_ylim((tmin[1],tmax[1]))
                ax.set_zlim((tmin[2],tmax[2]))
                plotdat = np.sum(density, axis=2)
                plotdat = plotdat / np.max(plotdat)
                plotx, ploty = np.mgrid[tmin[0]:tmax[0]:(tmax[0]-tmin[0])/divn[0],tmin[1]:tmax[1]:(tmax[1]-tmin[1])/divn[1]]
                ax.contour(plotx, ploty, plotdat, offset=tmin[2], zdir='z')
                plotdat = np.sum(density, axis=1)
                plotdat = plotdat / np.max(plotdat)
                plotx, plotz = np.mgrid[tmin[0]:tmax[0]:(tmax[0]-tmin[0])/divn[0],tmin[2]:tmax[2]:(tmax[2]-tmin[2])/divn[2]]
                ax.contour(plotx, plotdat, plotz, offset=tmax[1], zdir='y')
                plotdat = np.sum(density, axis=0)
                plotdat = plotdat / np.max(plotdat)
                ploty, plotz = np.mgrid[tmin[1]:tmax[1]:(tmax[1]-tmin[1])/divn[1],tmin[2]:tmax[2]:(tmax[2]-tmin[2])/divn[2]]
                ax.contour(plotdat, ploty, plotz, offset=tmin[0], zdir='x')
            outgraphname=outfilename[0]+'.png'
            #plt.savefig(outgraphname)
            #plt.show()
            plt.close()
            self.pmax=(scores.flatten()).max()
            self.div=divn
            self.para=[f0,f1,f2]
            self.mesh=[xx,yy,zz,density]
            self.flatten=[xx.flatten(),yy.flatten(),zz.flatten(),density.flatten()]
            self.dim=3
            self.xmin=[tmin[0],tmin[1],tmin[2]]
            self.xmax=[tmax[0],tmax[1],tmax[2]]
            self.dx=[(tmax[0]-tmin[0])/(divn[0]-1),(tmax[1]-tmin[1])/(divn[1]-1),(tmax[2]-tmin[2])/(divn[2]-1)]

            
    def bunpu_gene(self,xmin,xmax,xmean,xdev,div=[200],filename=''):
        dim=len(xmean)
        self.xmin=xmin
        self.xmax=xmax
        self.div=div
        if dim==1:
            s1,t1=genes(xmin[0],xmax[0],xmean[0],xdev[0],div[0],filename)
            self.para=[s1,t1]
            self.dim=1
            self.dx=[(xmax[0]-xmin[0])/(div[0]-1)]
            self.mesh=[s1,t1]
            self.flatten=[s1,t1]
            self.dim=1
            data=[s1,t1]
            self.pmax=max(t1)
        if dim==2:
            ZZ=[]
            #Z=[]
            s1,t1=genes(xmin[0],xmax[0],xmean[0],xdev[0],div[0])
            s2,t2=genes(xmin[1],xmax[1],xmean[1],xdev[1],div[1])
            base=(xmax[0]-xmin[0])*(xmax[1]-xmin[1])/((div[0]-1)*(div[1]-1))
            x = np.linspace(xmin[0],xmax[0],div[0])
            y = np.linspace(xmin[1],xmax[1],div[1])
            seki=0
            X,Y=np.meshgrid(x,y)
            case=3
            if case==0:
                for j in range(div[1]):
                    Z0=np.array(t1)*t2[j]
                    seki+=sum(Z0)*base
                    ZZ.append(Z0)
                    Z1=np.array(ZZ)
            else:
                haba=case
                Z0=np.zeros((div[1],div[0]))
                ZZ=np.zeros((div[1],div[0]))
                hasi=[t1[0]*t2[0],t1[div[0]-1]*t2[0],t1[div[0]-1]*t2[div[1]-1],t1[0]*t2[div[1]-1]]
                for i in range(div[1]):
                    for j in range(div[0]):
                        Z0[i][j]=t1[j]*t2[i]
                        k=0
                        ph=0
                        if i<=haba-1:
                            ph=((div[0]-1-j)*hasi[0]+j*hasi[3])/(div[0]-1)
                            k=haba-i-1
                        elif j<=haba-1:
                            ph=((div[1]-1-i)*hasi[0]+i*hasi[1])/(div[1]-1)
                            k=haba-j-1
                        elif j>=div[0]-haba:
                            ph=((div[1]-1-i)*hasi[1]+i*hasi[2])/(div[1]-1)
                            k=j+haba-div[0]
                        elif i>=div[1]-haba:
                            ph=((div[0]-1-j)*hasi[3]+j*hasi[2])/(div[0]-1)
                            k=i+haba-div[1]
                        if k!=0 or ph!=0:
                            ZZ[i][j]=((haba-k)*Z0[i][j]+k*ph)/haba
                        else:
                            ZZ[i][j]=Z0[i][j]
                        seki+=ZZ[i][j]*base
                Z1=np.array(ZZ)
                            
            if seki != 0:
                Z=Z1/seki
                print('menseki',seki)
            else:
                print('seki==0')
            self.pmax=(Z.flatten()).max()
            self.para=[x,y]
            self.mesh=[X,Y,Z]
            self.flatten=[X.flatten(),Y.flatten(),Z.flatten()]
            self.dim=2
            self.dx=[(xmax[0]-xmin[0])/(div[0]-1),(xmax[1]-xmin[1])/(div[1]-1)]
            data=[X.flatten(),Y.flatten(),Z.flatten()]
        elif dim==3:
            WWW=[]
            seki=0
            div2=div
            s1,t1=genes(xmin[0],xmax[0],xmean[0],xdev[0],div[0])
            s2,t2=genes(xmin[1],xmax[1],xmean[1],xdev[1],div[1])
            s3,t3=genes(xmin[2],xmax[2],xmean[2],xdev[2],div[2])
            base=(xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])/((div[0]-1)*(div[1]-1)*(div[2]-1))
            x = np.linspace(xmin[0],xmax[0],div[0])
            y = np.linspace(xmin[1],xmax[1],div[1])
            z = np.linspace(xmin[2],xmax[2],div[2])
            X,Y,Z=np.meshgrid(x,y,z)
            case=3
            if case==0:
                for j in range(div[1]):
                    WW=[]
                    for k in range(div[2]):
                        W0=np.array(t1)*t2[j]*t3[k]
                        seki+=sum(W0)*base
                        WW.append(W0)
                    W1=np.array(WW)
                    WWW.append(W1)
                W2=np.array(WWW)
            else:
                haba=case
                Z0=np.zeros((div[2],div[1],div[0]))
                ZZ=np.zeros((div[2],div[1],div[0]))
                hasi=[t1[0]*t2[0]*t3[0],t1[div[0]-1]*t2[0]*t3[0],t1[0]*t2[div[1]-1]*t3[0],t1[0]*t2[0]*t3[div[2]-1],t1[div[0]-1]*t2[div[1]-1]*t3[0],\
                      t1[0]*t2[div[1]-1]*t3[div[2]-1],t1[div[0]-1]*t2[0]*t3[div[2]-1],t1[div[0]-1]*t2[div[1]-1]*t3[div[2]-1]]
                for j in range(div[1]):
                    WW=[]
                    for k in range(div[2]):
                        for i in range(div[0]):
                            Z0[k][j][i]=t1[i]*t2[j]*t3[k]
                            ph=0
                            g=0
                            h=0
                            if i<=haba-1 and k<=haba-1:
                                ph=((div[0]-1-j)*hasi[0]+j*hasi[2])/(div[0]-1)
                                g=haba-k-1
                                h=haba-i-1
                            elif j<=haba-1 and k<=haba-1:
                                ph=((div[1]-1-i)*hasi[0]+i*hasi[1])/(div[1]-1)
                                g=haba-k-1
                                h=haba-j-1
                            elif j<=haba-1 and i<=haba-1:
                                ph=((div[1]-1-k)*hasi[0]+k*hasi[3])/(div[2]-1)
                                g=haba-i-1
                                h=haba-j-1
                                
                            elif j<=haba-1 and i>=div[1]-haba:
                                ph=((div[2]-1-k)*hasi[1]+k*hasi[6])/(div[2]-1)
                                g=haba-j-1
                                h=i+haba-div[0]
                            elif j<=haba-1 and k>=div[2]-haba:
                                ph=((div[2]-1-i)*hasi[3]+i*hasi[6])/(div[0]-1)
                                g=haba-j-1
                                h=k+haba-div[2]
                            elif i<=haba-1 and k>=div[2]-haba:
                                ph=((div[1]-1-j)*hasi[3]+j*hasi[5])/(div[1]-1)
                                g=haba-i-1
                                h=k+haba-div[2]
                            elif i<=haba-1 and j>=div[1]-haba:
                                ph=((div[2]-1-k)*hasi[2]+k*hasi[5])/(div[2]-1)
                                g=haba-i-1
                                h=j+haba-div[1]
                            elif k<=haba-1 and j>=div[1]-haba:
                                ph=((div[0]-1-i)*hasi[2]+i*hasi[4])/(div[0]-1)
                                g=haba-k-1
                                h=j+haba-div[1]
                            elif k<=haba-1 and i>=div[1]-haba:
                                ph=((div[1]-1-j)*hasi[1]+j*hasi[4])/(div[1]-1)
                                g=haba-k-1
                                h=i+haba-div[0]
                                
                            elif j>=div[1]-haba and k>=div[2]-haba:
                                ph=((div[1]-1-i)*hasi[2]+i*hasi[4])/(div[1]-1)
                                g=j+haba-div[1]
                                h=k+haba-div[2]
                            elif i>=div[0]-haba and k>=div[2]-haba:
                                ph=((div[1]-1-j)*hasi[1]+j*hasi[4])/(div[1]-1)
                                g=i+haba-div[0]
                                h=k+haba-div[2]
                            elif i>=div[0]-haba and j>=div[1]-haba:
                                ph=((div[1]-1-k)*hasi[4]+k*hasi[7])/(div[1]-1)
                                g=i+haba-div[0]
                                h=j+haba-div[1]
                            if h!=0 or ph!=0:
                                if g>=h:
                                    ZZ[k][j][i]=((haba-h)*Z0[k][j][i]+h*ph)/haba
                                else:
                                    ZZ[k][j][i]=((haba-g)*Z0[k][j][i]+g*ph)/haba
                            else:
                                ZZ[k][j][i]=Z0[k][j][i]
                            seki+=ZZ[k][j][i]*base
                W2=np.array(ZZ)
            W=W2/seki
            print('menseki',seki)

            self.pmax=(Z.flatten()).max()
            self.para=[x,y,z]
            self.mesh=[X,Y,Z,W]
            self.flatten=[X.flatten(),Y.flatten(),Z.flatten(),W.flatten()]
            self.div=div
            self.dim=3
            self.dx=[(xmax[0]-xmin[0])/(div[0]-1),(xmax[1]-xmin[1])/(div[1]-1),(xmax[2]-xmin[2])/(div[2]-1)]
            data=[X.flatten(),Y.flatten(),Z.flatten(),W.flatten()]
        else:
            z=0
    
        
    def bunpu_add(self,other,divt0=[],corel=[0],shw=0,cormap=0):
            
        r=corel
        max_x=0
        max_y=0
        max_p=0
        cmap=0
        X,X1,X0,x0min,x0max,dim,dx,divx,m,y0min,y0max,divt0,flag0=self.bunpu_filter(other,0,divt0)#today,X(para)だけが3重カッコになっている
        if flag0!=0:
            if flag0!=3:
                Y=m
                divt=divx
                flagx=0
                flagy=0
                if dim==1:
                    y0=np.linspace(y0min[0],y0max[0],len(X0[0]))
                    t=[X0[0][i]+y0[i] for i in range(len(X0[0]))]
                    t1=t
                    t0=t
                    tmin=[x0min[0]+y0min[0]]
                    tmax=[x0max[0]+y0max[0]]
                    p1=[X1[1][i]*(x0max[0]-x0min[0])/(tmax[0]-tmin[0]) for i in range(len(X1[1]))]
                    p0=p1

                    divtt=divt[0]
                    max_p=max(p1)
                    dt=np.array(dx)*(tmax[0]-tmin[0])/(x0max[0]-x0min[0])
                elif dim==2:
                    t=[]
                    t1=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    y0=[np.linspace(y0min[0],y0max[0],len(X0[0])),np.linspace(y0min[1],y0max[1],len(X0[1]))]
                    t.append([X[0][i]+y0[0][i] for i in range(len(X[0]))])
                    t.append([X[1][i]+y0[1][i] for i in range(len(X[1]))])
                    tx,ty=np.meshgrid(t[0],t[1])
                    t0=[tx,ty]
                    t1=[tx.flatten(),ty.flatten()]
                    tmin.append(x0min[0]+y0min[0])
                    tmin.append(x0min[1]+y0min[1])
                    tmax.append(x0max[0]+y0max[0])
                    tmax.append(x0max[1]+y0max[1])
                    tmp0=(x0max[0]-x0min[0])*(x0max[1]-x0min[1])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1]))
                    p1=X1[2]*tmp0
                    p0=X0[2]*tmp0
                    max_p=max(p1)
                    divtt=divt[0]*divt[1]
                    dt=np.array(dx)/tmp0
                elif dim==3:
                    t=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    y0=[np.linspace(y0min[0],y0max[0],len(X0[0])),np.linspace(y0min[1],y0max[1],len(X0[1])),np.linspace(y0min[2],y0max[2],len(X0[2]))]
                    t.append([X[0][i]+y0[0][i] for i in range(len(X[0]))])
                    t.append([X[1][i]+y0[1][i] for i in range(len(X[1]))])
                    t.append([X[2][i]+y0[2][i] for i in range(len(X[2]))])
                    tx,ty,tz=np.meshgrid(t[0],t[1],t[2])
                    t0=[tx,ty,tz]
                    t1=[tx.flatten(),ty.flatten(),tz.flatten()]
                    tmin.append(x0min[0]+y0min[0])
                    tmin.append(x0min[1]+y0min[1])
                    tmin.append(x0min[2]+y0min[2])
                    tmax.append(x0max[0]+y0max[0])
                    tmax.append(x0max[1]+y0max[1])
                    tmax.append(x0max[2]+y0max[2])
                    tmp0=(x0max[0]-x0min[0])*(x0max[1]-x0min[1])*(x0max[2]-x0min[2])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1])*(tmax[2]-tmin[2]))
                    p1=X1[3]*tmp0
                    p0=X0[3]*tmp0
                    divtt=divt[0]*divt[1]*divt[2]
                    max_p=max(p1)
                    dt=np.array(dx)/tmp0
            else:#yがベクトル
                Y=other
                divt=divx
                flagx=0
                flagy=0
                if dim==1:
                    
                    t=[X0[0][i]+Y[0] for i in range(len(X0[0]))]
                    t1=t
                    t0=t
                    p1=X1[1]
                    p0=X1[1]
                    tmin=[x0min[0]+Y[0]]
                    tmax=[x0max[0]+Y[0]]
                    divtt=divt[0]
                    max_p=max(p1)
                elif dim==2:
                    t=[]
                    t1=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    t.append([X[0][i]+Y[0] for i in range(len(X[0]))])
                    t.append([X[1][i]+Y[1] for i in range(len(X[1]))])
                    t1.append([X1[0][i]+Y[0] for i in range(len(X1[0]))])
                    t1.append([X1[1][i]+Y[1] for i in range(len(X1[1]))])
                    t0.append([[X0[0][i][j]+Y[0] for j in range(len(X0[0][0]))] for i in range(len(X0[0]))])
                    t0.append([[X0[1][i][j]+Y[1] for j in range(len(X0[1][0]))] for i in range(len(X0[1]))])
                    tmin.append(x0min[0]+Y[0])
                    tmin.append(x0min[1]+Y[1])
                    tmax.append(x0max[0]+Y[0])
                    tmax.append(x0max[1]+Y[1])
                    p1=X1[2]
                    p0=X0[2]
                    max_p=max(p1)
                    divtt=divt[0]*divt[1]
                elif dim==3:
                    t=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    t.append([X[0][i]+Y[0] for i in range(len(X[0]))])
                    t.append([X[1][i]+Y[1] for i in range(len(X[1]))])
                    t.append([X[2][i]+Y[2] for i in range(len(X[2]))])
                    t1.append([X1[0][i]+Y[0] for i in range(len(X1[0]))])
                    t1.append([X1[1][i]+Y[1] for i in range(len(X1[1]))])
                    t1.append([X1[2][i]+Y[2] for i in range(len(X1[2]))])
                    t0.append([[[X0[0][i][j][k]+Y[0] for k in range(len(X0[0][0][0]))] for j in range(len(X0[0][0]))] for i in range(len(X0[0]))])
                    t0.append([[[X0[1][i][j][k]+Y[1] for k in range(len(X0[1][0][0]))] for j in range(len(X0[1][0]))] for i in range(len(X0[1]))])
                    t0.append([[[X0[2][i][j][k]+Y[2] for k in range(len(X0[2][0][0]))] for j in range(len(X0[2][0]))] for i in range(len(X0[2]))])
                    tmin.append(x0min[0]+Y[0])
                    tmin.append(x0min[1]+Y[1])
                    tmin.append(x0min[2]+Y[2])
                    tmax.append(x0max[0]+Y[0])
                    tmax.append(x0max[1]+Y[1])
                    tmax.append(x0max[2]+Y[2])
                    p1=X1[3]
                    p0=X0[3]
                    divtt=divt[0]*divt[1]*divt[2]
                    max_p=max(p1)
                dt=dx
            menseki=0
            for ti in range(divtt):
                if dim==1:
                    menseki+=p1[ti]*dt[0]
                if dim==2:
                    menseki+=p1[ti]*(dt[0]*dt[1])
                elif dim==3:
                    menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
        else:
            #y
            Y=other.flatten
            yf,py0,y0,dim,divy,dy,y0min,y0max=other.kaiseki()
            if other.xmin==[] or other.xmax==[] or other.dim==0 or other.div==[]:
                
                other.xmin=y0min
                other.xmax=y0max
                other.dim=dim
                other.div=divy
                other.dx=dy
            else:
                y0min=other.xmin
                y0max=other.xmax
                dim=other.dim
                divy=other.div
                dy=other.dx
            #
            #確率値が小さい場合、10000倍しておいてmensekiも10000倍する>*1
            #chosei=1000
            chosei=1
            flagx=0
            flagy=0
            if self.pmax<=0.001:
                flagx=1
            if other.pmax<=0.001:
                flagy=1
            #
            #
            div=self.div
            if self.div==[]:
                div=len(x0[0])
            xmin=[]
            xmax=[]
            ymin=[]
            ymax=[]
            x=[]
            y=[]
            
            filename='x+y_log'
            logname=filename+'.csv'
            log = open(logname,'w' , encoding='shift_jis' )
            info=['tパラメータ','p値']
            writer = csv.writer(log, lineterminator='\n')
            writer.writerow(info)
            #各次元のmin,maxを抽出
            for i in range(dim):
                xmin.append(min(X[i]))
                xmax.append(max(X[i]))
                ymin.append(min(Y[i]))
                ymax.append(max(Y[i]))
            #
            divt=[]
            dt=[]
            tmin=np.array(xmin)+np.array(ymin)
            tmax=np.array(xmax)+np.array(ymax)
            for i in range(dim):
                if dx[i]<=dy[i]:
                    divt.append(divx[i])
                else:
                    divt.append(divy[i])
                dt.append((tmax[i]-tmin[i])/(divt[i]-1))
            x2,divx2,dx2=self.divide(dt)
            y2,divy2,dy2=other.divide(dt)
            t=[]
            t1=[]#flatten
            p1=[]
            t0=[]#
            p0=[]
            lastp=[]
            #dx=[]
            #dy=[]
            x=[]#X:演算前のパラメータX
            y=[]#Y:演算前のパラメータY
            menseki=0
            divxt=1
            if dim==1:
                divtt=divt[0]
                t=np.linspace(tmin[0],tmax[0],divt[0])#11/10
                t1=t
                t0=t
                p1=np.zeros(divt[0])
                p0=p1
                lastp=np.zeros(divt[0])
                #>*1
                if flagx==1:
                    for i in range(len(x2[1])):
                        x2[1][i]=x2[1][i]*chosei
                if flagy==1:
                    for i in range(len(y2[1])):
                        y2[1][i]=y2[1][i]*chosei
                zfun_smooth_rbf = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0)  # default smooth=0 for interpolation
            elif dim ==2:
                divtt=divt[0]*divt[1]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1])]
                tx,ty=np.meshgrid(t[0],t[1])#11/10
                t1=[tx.flatten(),ty.flatten()]
                t0=[tx,ty]
                p1=np.zeros(divtt)
                p0=np.zeros((divt[1],divt[0]))#7/2
                lastp=np.zeros(divtt)
                #>*1
                if flagx==1:
                    x2[2]=x2[2]*chosei
                if flagy==1:
                    y2[2]=y2[2]*chosei
                zfun_smooth_rbf = interpolate.Rbf(y2[0], y2[1], y2[2], function='linear', smooth=0)
            elif dim==3:
                divtt=divt[0]*divt[1]*divt[2]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1]),np.linspace(tmin[2],tmax[2],divt[2])]
                tx,ty,tz=np.meshgrid(t[0],t[1],t[2])#11/10
                t1=[tx.flatten(),ty.flatten(),tz.flatten()]
                t0=[tx,ty,tz]
                p1=[0.]*divtt
                p0=np.zeros((divt[2],divt[1],divt[0]))
                lastp=np.zeros(divtt)
                if flagx==1:
                    x2[3]=x2[3]*chosei
                if flagy==1:
                    y2[3]=y2[3]*chosei
                zfun_smooth_rbf = interpolate.Rbf(y2[0], y2[1], y2[2], y2[3], function='linear', smooth=0)
            for i in range(dim):
                x.append(np.arange(xmin[i],xmax[i],dt[i]))
                if max(x[i])<xmax[i]:
                    np.append(x,xmax[i])
                y.append(np.arange(ymin[i],ymax[i],dt[i]))
                if max(y[i])<ymax[i]:
                    np.append(y,ymax[i])
                divxt*=divx2[i]
            if len(corel)==1 and dim==2:
                corel=[corel[0],corel[0]]
            elif len(corel)==1 and dim==3:
                corel=[corel[0],corel[0],corel[0]]
            if type(cormap)==list:
                cmap=np.zeros((divtt,divxt))
            elif type(cormap)==int and cormap==1:
                cmap=np.zeros((divtt,divxt))
            else:
                cmap=0
            for ti in range(divtt):
                if dim>=2:
                    tq,tmod=divmod(ti,divt[0])
                if dim==3:
                    tq,tmod=divmod(ti,divt[0])
                    tq1,tq2=divmod(tq,divt[1])
                #

                for xi in range(divxt):
                    flag=np.zeros(dim)
                    yh=np.zeros(dim)
                    if dim>=2:
                        xq,xmod=divmod(xi,divx2[0])
                    if dim==3:
                        xq,xmod=divmod(xi,divx2[0])
                        xq1,xq2=divmod(tq,divx2[1])
                    #
                    for i in range(dim):
                        if dim==1:
                            #yh[i]=round(t1[ti]-x2[i][xi],10)
                            yh[i]=np.round(t1[ti]-x2[i][xi],decimals=9)
                            flag[i]=cor2ran(x2[i][xi],xmin[i],xmax[i],yh[i],ymin[i],ymax[i],corel[i])
                        else:
                            yh[i]=t1[i][ti]-(x2[i].flatten())[xi]
                            flag[i]=cor2ran((x2[i].flatten())[xi],xmin[i],xmax[i],yh[i],ymin[i],ymax[i],corel[i])
                    #
                    if all(flag)==1.0:
                        if dim==1:
                            pyh = zfun_smooth_rbf(yh[0])
                            ptmp2=pyh*dt[0]*x2[1][xi]
                            if type(cormap)==list:
                                ptmp2=ptmp2*cormap[ti][xi]
                                cmap[ti][xi]=ptmp2
                            elif type(cormap)==int and cormap==1:
                                cmap[ti][xi]=ptmp2
                            p1[ti]+=ptmp2
                            p0[ti]=p1[ti]
                            #
                        elif dim==2:
                            ptmp = zfun_smooth_rbf(yh[0], yh[1])
                            if ptmp>0:
                                pyh = ptmp 
                            else:
                                pyh = 0
                            ptmp2=pyh*dt[0]*dt[1]*(x2[2].flatten())[xi]
                            if type(cormap)==list:
                                ptmp2=ptmp2*cormap[ti][xi]
                                cmap[ti][xi]=ptmp2
                            elif type(cormap)==int and cormap==1:
                                cmap[ti][xi]=ptmp2
                            p1[ti]+=ptmp2
                            p0[tq][tmod]+= ptmp2
                            if max_x <= (x2[2].flatten())[xi]:
                                max_x=(x2[2].flatten())[xi]
                            if max_y <= pyh:
                                max_y=pyh
                            if max_p <= p1[ti]:
                                max_p=p1[ti]
                        elif dim==3:
                            ptmp = zfun_smooth_rbf(yh[0], yh[1], yh[2])
                            if ptmp>0:
                                pyh = ptmp 
                            else:
                                pyh = 0
                            ptmp2=pyh*dt[0]*dt[1]*dt[2]*(x2[3].flatten())[xi]
                            if type(cormap)==list:
                                ptmp2=ptmp2*cormap[ti][xi]
                                cmap[ti][xi]=ptmp2
                            elif type(cormap)==int and cormap==1:
                                cmap[ti][xi]=ptmp2
                            p1[ti]+= ptmp2
                            p0[tq1][tq2][tmod]+= ptmp2
                            if max_x <= (x2[3].flatten())[xi]:
                                max_x=(x2[3].flatten())[xi]
                            if max_y <= pyh:
                                max_y=pyh
                            if max_p <= p1[ti]:
                                max_p=p1[ti]
                        else:
                            print('error')

                if dim==1:
                    menseki+=p1[ti]*dt[0]
                if dim==2:
                    menseki+=p1[ti]*(dt[0]*dt[1])
                elif dim==3:
                    menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
                info=[ti,max_x,max_y,max_p]
                writer = csv.writer(log, lineterminator='\n')
                writer.writerow(info)
            
            log.close()
        
        if menseki>0:
            if dim==1:
                p00=p1/menseki
                p=p0/menseki
            else:
                p=p1/menseki
                p00=p0/menseki
        else:
            if dim==1:
                p00=p1
                p=p0
            else:
                p=p1
                p00=p0
        #
        if flagx==1:
            menseki=menseki/chosei
        if flagy==1:
            menseki=menseki/chosei
        if shw==0:
            #flag:0 分布、1 otherが狭い、2　selfが狭い、3　otherがベクトル
            flaginf=['bunpu+bunpu','bunpu+lean','lean+bunpu','bunpu+vector']
            print(flaginf[flag0],menseki)
        #
        if dim==1:
            W=[t]
            W0=[t0,p00]
            W1=[t1,p]
        elif dim==2:
            W=[t[0],t[1]]
            W0=[t0[0],t0[1],p00]#mesh
            W1=[t1[0],t1[1],p]#flatten
        elif dim==3:
            W=[t[0],t[1],t[2]]
            W0=[t0[0],t0[1],t0[2],p00]
            W1=[t1[0],t1[1],t1[2],p]

        WW=bunpu()
        WW.cmap=cmap
        WW.para=W
        WW.flatten=W1
        WW.mesh=W0
        WW.div=divt
        WW.dx=dt
        WW.dim=dim
        WW.xmin=tmin
        WW.xmax=tmax
        WW.pmax=p.max()
        return WW
    
    def bunpu_sub(self,other,divt0=[],corel=[0],shw=0):
        r=corel        
        max_x=0
        max_y=0
        max_p=0
        X,X1,X0,x0min,x0max,dim,dx,divx,m,y0min,y0max,divt0,flag0=self.bunpu_filter(other,0,divt0)
        if flag0!=0:
            if flag0!=3:
                Y=m
                divt=divx
                flagx=0
                flagy=0
                if dim==1:
                    y0=np.linspace(y0min[0],y0max[0],len(X0[0]))
                    t=[X0[0][i]-y0[i] for i in range(len(X0[0]))]
                    t1=t
                    t0=t
                    tmin=[x0min[0]-y0min[0]]
                    tmax=[x0max[0]-y0max[0]]
                    p1=[X1[1][i]*(x0max[0]-x0min[0])/(tmax[0]-tmin[0]) for i in range(len(X1[1]))]
                    p0=p1
                    divtt=divt[0]
                    max_p=max(p1)
                    dt=np.array(dx)*(tmax[0]-tmin[0])/(x0max[0]-x0min[0])
                elif dim==2:
                    t=[]
                    t1=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    y0=[np.linspace(y0min[0],y0max[0],len(X0[0])),np.linspace(y0min[1],y0max[1],len(X0[1]))]
                    t.append([X[0][i]-y0[0][i] for i in range(len(X[0]))])# :para
                    t.append([X[1][i]-y0[1][i] for i in range(len(X[1]))])
                    tx,ty=np.meshgrid(t[0],t[1])
                    t0=[tx,ty]
                    t1=[tx.flatten(),ty.flatten()]
                    tmin.append(x0min[0]-y0min[0])
                    tmin.append(x0min[1]-y0min[1])
                    tmax.append(x0max[0]-y0max[0])
                    tmax.append(x0max[1]-y0max[1])
                    tmp0=(x0max[0]-x0min[0])*(x0max[1]-x0min[1])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1]))
                    p1=X1[2]*tmp0
                    p0=X0[2]*tmp0
                    max_p=max(p1)
                    divtt=divt[0]*divt[1]
                    dt=np.array(dx)/tmp0
                elif dim==3:
                    t=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    y0=[np.linspace(y0min[0],y0max[0],len(X0[0])),np.linspace(y0min[1],y0max[1],len(X0[1])),np.linspace(y0min[2],y0max[2],len(X0[2]))]
                    t.append([X[0][i]-y0[0][i] for i in range(len(X[0]))])
                    t.append([X[1][i]-y0[1][i] for i in range(len(X[1]))])
                    t.append([X[2][i]-y0[2][i] for i in range(len(X[2]))])
                    tx,ty,tz=np.meshgrid(t[0],t[1],t[2])
                    t0=[tx,ty,tz]
                    t1=[tx.flatten(),ty.flatten(),tz.flatten()]
                    tmin.append(x0min[0]-y0min[0])
                    tmin.append(x0min[1]-y0min[1])
                    tmin.append(x0min[2]-y0min[2])
                    tmax.append(x0max[0]-y0max[0])
                    tmax.append(x0max[1]-y0max[1])
                    tmax.append(x0max[2]-y0max[2])
                    tmp0=(x0max[0]-x0min[0])*(x0max[1]-x0min[1])*(x0max[2]-x0min[2])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1])*(tmax[2]-tmin[2]))
                    p1=X1[3]*tmp0
                    p0=X0[3]*tmp0
                    divtt=divt[0]*divt[1]*divt[2]
                    max_p=max(p1)
                    dt=np.array(dx)/tmp0
            else:
                Y=other
                divt=divx
                flagx=0
                flagy=0
                if dim==1:
                    t=[X0[0][i]-Y[0] for i in range(len(X0[0]))]
                    t1=t
                    t0=t
                    p1=X1[1]
                    p0=X0[1]
                    tmin=[x0min[0]-Y[0]]
                    tmax=[x0max[0]-Y[0]]
                    divtt=divt[0]
                    max_p=max(p1)
                elif dim==2:
                    t=[]
                    t1=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    t.append([X[0][i]-Y[0] for i in range(len(X[0]))])
                    t.append([X[1][i]-Y[1] for i in range(len(X[1]))])
                    t1.append([X1[0][i]-Y[0] for i in range(len(X1[0]))])
                    t1.append([X1[1][i]-Y[1] for i in range(len(X1[1]))])
                    t0.append([[X0[0][i][j]-Y[0] for j in range(len(X0[0][0]))] for i in range(len(X0[0]))])
                    t0.append([[X0[1][i][j]-Y[1] for j in range(len(X0[1][0]))] for i in range(len(X0[1]))])
                    tmin.append(x0min[0]-Y[0])
                    tmin.append(x0min[1]-Y[1])
                    tmax.append(x0max[0]-Y[0])
                    tmax.append(x0max[1]-Y[1])
                    p1=X1[2]
                    p0=X0[2]
                    divtt=divt[0]*divt[1]
                    max_p=max(p1)
                elif dim==3:
                    t=[]
                    t0=[]
                    tmin=[]
                    tmax=[]
                    t.append([X[0][i]-Y[0] for i in range(len(X[0]))])
                    t.append([X[1][i]-Y[1] for i in range(len(X[1]))])
                    t.append([X[2][i]-Y[2] for i in range(len(X[2]))])
                    t1.append([X1[0][i]-Y[0] for i in range(len(X1[0]))])
                    t1.append([X1[1][i]-Y[1] for i in range(len(X1[1]))])
                    t1.append([X1[2][i]-Y[2] for i in range(len(X1[2]))])
                    t0.append([[[X0[0][i][j][k]-Y[0] for k in range(len(X0[0][0]))] for j in range(len(X0[0]))] for i in range(len(X0[0]))])
                    t0.append([[[X0[1][i][j][k]-Y[1] for k in range(len(X0[1][0]))] for j in range(len(X0[1]))] for i in range(len(X0[1]))])
                    t0.append([[[X0[2][i][j][k]-Y[2] for k in range(len(X0[2][0]))] for j in range(len(X0[2]))] for i in range(len(X0[2]))])
                    tmin.append(x0min[0]-Y[0])
                    tmin.append(x0min[1]-Y[1])
                    tmin.append(x0min[2]-Y[2])
                    tmax.append(x0max[0]-Y[0])
                    tmax.append(x0max[1]-Y[1])
                    tmax.append(x0max[2]-Y[2])
                    p1=X1[3]
                    p0=X0[3]
                    max_p=max(p1)
                    divtt=divt[0]*divt[1]*divt[2]
                #divt=divx
                dt=dx
            menseki=0
            for ti in range(divtt):
                if dim==1:
                    menseki+=p1[ti]*dt[0]
                if dim==2:
                    menseki+=p1[ti]*(dt[0]*dt[1])
                elif dim==3:
                    menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
        else:
            Y=other.flatten
            yf,py0,y0,dim,divy,dy,y0min,y0max=other.kaiseki()
            if other.xmin==[] or other.xmax==[] or other.dim==0 or other.div==[]:
                
                other.xmin=y0min
                other.xmax=y0max
                other.dim=dim
                other.div=divy
                other.dx=dy
            else:
                y0min=other.xmin
                y0max=other.xmax
                dim=other.dim
                divy=other.div
                dy=other.dx
                
            chosei=1
            flagx=0
            flagy=0
            if self.pmax<=0.001:
                flagx=1
            if other.pmax<=0.001:
                flagy=1
            #
            xmin=[]
            xmax=[]
            ymin=[]
            ymax=[]
            x=[]
            y=[]
            
            filename='x-y_log'
            logname=filename+'.csv'
            log = open(logname,'w' , encoding='shift_jis' )
            info=['tパラメータ','p値']
            writer = csv.writer(log, lineterminator='\n')
            writer.writerow(info)
            for i in range(dim):
                xmin.append(min(X[i]))
                xmax.append(max(X[i]))
                ymin.append(min(Y[i]))
                ymax.append(max(Y[i]))
            #
            divt=[]
            dt=[]
            tmin=np.array(xmin)-np.array(ymax)
            tmax=np.array(xmax)-np.array(ymin)
            for i in range(dim):
                if dx[i]<=dy[i]:
                    divt.append(divx[i])
                else:
                    divt.append(divy[i])
                dt.append((tmax[i]-tmin[i])/(divt[i]-1))
            x2,divx2,dx2=self.divide(dt)
            y2,divy2,dy2=other.divide(dt)

            t=[]#
            t1=[]
            p1=[]
            t0=[]#
            p0=[]
            lastp=[]
            x=[]
            y=[]
            menseki=0
            divxt=1
            if dim==1:
                divtt=divt[0]
                t=np.linspace(tmin[0],tmax[0],divt[0])
                t1=t
                t0=t
                p1=np.zeros(divt[0])
                p0=p1
                lastp=np.zeros(divt[0])
                if flagx==1:
                    x2[1]=x2[1]*chosei
                if flagy==1:
                    y2[1]=y2[1]*chosei
                divxt=divx2[0]
                zfun_smooth_rbf = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0)
            elif dim ==2:
                divtt=divt[0]*divt[1]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1])]
                tx,ty=np.meshgrid(t[0],t[1])
                t1=[tx.flatten(),ty.flatten()]
                t0=[tx,ty]
                p1=np.zeros(divtt)
                p0=np.zeros((divt[1],divt[0]))
                lastp=np.zeros(divtt)
                if flagx==1:
                    x2[2]=x2[2]*chosei
                if flagy==1:
                    y2[2]=y2[2]*chosei
                divxt=divx2[0]*divx2[1]
                zfun_smooth_rbf = interpolate.Rbf(y2[0], y2[1], y2[2], function='linear', smooth=0)  
                
            elif dim==3:
                divtt=divt[0]*divt[1]*divt[2]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1]),np.linspace(tmin[2],tmax[2],divt[2])]
                tx,ty,tz=np.meshgrid(t[0],t[1],t[2])
                t1=[tx.flatten(),ty.flatten(),tz.flatten()]
                t0=[tx,ty,tz]
                p1=[0.]*divtt
                p0=np.zeros((divt[2],divt[1],divt[0]))
                lastp=np.zeros(divtt)
                if flagx==1:
                    x2[3]=x2[3]*chosei
                if flagy==1:
                    y2[3]=y2[3]*chosei
                divxt=divx2[0]*divx2[1]*divx2[2]
                zfun_smooth_rbf = interpolate.Rbf(y2[0], y2[1], y2[2], y2[3], function='linear', smooth=0)
            for i in range(dim):
                x.append(np.arange(xmin[i],xmax[i],dt[i]))
                if max(x[i])<xmax[i]:
                    np.append(x,xmax[i])
                y.append(np.arange(ymin[i],ymax[i],dt[i]))
                if max(y[i])<ymax[i]:
                    np.append(y,ymax[i])
            if len(corel)==1 and dim==2:
                corel=[corel[0],corel[0]]
            elif len(corel)==1 and dim==3:
                corel=[corel[0],corel[0],corel[0]]

            for ti in range(divtt):#

                if dim>=2:
                    tq,tmod=divmod(ti,divt[0])
                if dim==3:
                    tq,tmod=divmod(ti,divt[0])    
                    tq1,tq2=divmod(tq,divt[1])
                #
                for xi in range(divxt):#
                    flag=np.zeros(dim)
                    yh=np.zeros(dim)
                    if dim>=2:
                        xq,xmod=divmod(xi,divx2[0])
                    if dim==3:
                        xq,xmod=divmod(xi,divx2[0])
                        xq1,xq2=divmod(tq,divx2[1])
                    #
                    for i in range(dim):
                        if dim==1:
                            yh[i]=round(x2[i][xi]-t1[ti],10)
                            flag[i]=cor2ran(x2[i][xi],xmin[i],xmax[i],yh[i],ymin[i],ymax[i],corel[i])
                        else:
                            if corel==[0]:
                                for j in range(dim-1):
                                    corel.append(0)
                            yh[i]=(x2[i].flatten())[xi]-t1[i][ti]
                            flag[i]=cor2ran((x2[i].flatten())[xi],xmin[i],xmax[i],yh[i],ymin[i],ymax[i],corel[i])
                    #
                    if all(flag)==1.0:
                        if dim==1:
                            pyh = zfun_smooth_rbf(yh[0])
                            p1[ti]+=pyh*dt[0]*x2[1][xi]
                            p0[ti]=p1[ti]
                            #
                        elif dim==2:
                            ptmp = zfun_smooth_rbf(yh[0], yh[1])
                            if ptmp>0:
                                pyh = ptmp 
                            else:
                                pyh = 0
                            p1[ti]+=pyh*dt[0]*dt[1]*(x2[2].flatten())[xi]
                            p0[tq][tmod]+= pyh*dt[0]*dt[1]*(x2[2].flatten())[xi]
                            if max_x <= (x2[2].flatten())[xi]:
                                max_x=(x2[2].flatten())[xi]
                            if max_y <= pyh:
                                max_y=pyh
                            if max_p <= p1[ti]:
                                max_p=p1[ti]
                        elif dim==3:
                            ptmp = zfun_smooth_rbf(yh[0], yh[1], yh[2])  
                            if ptmp>0:
                                pyh = ptmp  
                            else:
                                pyh = 0
                            p1[ti]+= pyh*dt[0]*dt[1]*dt[2]*(x2[3].flatten())[xi]
                            p0[tq1][tq2][tmod]+= pyh*dt[0]*dt[1]*dt[2]*(x2[3].flatten())[xi]
                            if max_x <= (x2[3].flatten())[xi]:
                                max_x=(x2[3].flatten())[xi]
                            if max_y <= pyh:
                                max_y=pyh
                            if max_p <= p1[ti]:
                                max_p=p1[ti]
                        else:
                            print('error')
                if dim==1:
                    menseki+=p1[ti]*dt[0]
                elif dim==2:
                    menseki+=p1[ti]*(dt[0]*dt[1])
                elif dim==3:
                    menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
                info=[ti,max_x,max_y,max_p]
                writer = csv.writer(log, lineterminator='\n')
                writer.writerow(info)
            log.close()

        if menseki>0:
            if dim==1:
                p00=p1/menseki
                p=p0/menseki
            else:
                p=p1/menseki
                p00=p0/menseki
        else:
            if dim==1:
                p00=p1
                p=p0
            else:
                p=p1
                p00=p0
        #>*1
        if flagx==1:
            menseki=menseki/chosei
        if flagy==1:
            menseki=menseki/chosei
        if shw==0:
            flaginf=['bunpu-bunpu','bunpu-lean','lean-bunpu','bunpu-vector']
            print(flaginf[flag0],menseki)
        #
        if dim==1:
            W=[t]
            W0=[t0,p00]
            W1=[t1,p]
        elif dim==2:
            W=[t[0],t[1]]
            W0=[t0[0],t0[1],p00]
            W1=[t1[0],t1[1],p]
        elif dim==3:
            W=[t[0],t[1],t[2]]
            W0=[t0[0],t0[1],t0[2],p00]
            W1=[t1[0],t1[1],t1[2],p]

        WW=bunpu()
        WW.para=W
        WW.flatten=W1
        WW.mesh=W0
        WW.div=divt
        WW.dx=dt
        WW.dim=dim
        WW.xmin=tmin
        WW.xmax=tmax
        WW.pmax=max_p
        return WW



    def bunpu_simu(self,other,dt,corel=0,vdiv=[10],shw=0,cormap=0):
        X0=self.mesh
        X1=self.flatten
        dim=self.dim

        div=[]
        if dim==1:
            div.append(len(X0[0]))
        elif dim==2:
            div.append(len(X0[0]))
            div.append(len(X0[0][0]))
        elif dim==3:
            div.append(len(X0[0]))
            div.append(len(X0[0][0]))
            div.append(len(X0[0][0][0]))
            
        
        #
        Z=bunpu()
        Z=self.bunpu_add(other,div,corel,shw,cormap)
        zmin=Z.xmin
        zmax=Z.xmax
        Z0=Z.mesh
        Z1=Z.flatten
        tmax=[]
        tmin=[]
        divt=1#
        divz=1
        ndt=1/dt
        for i in range(dim):
            
            tmax.append((self.xmax[i]*(ndt-1)+zmax[i])/ndt)
            tmin.append((self.xmin[i]*(ndt-1)+zmin[i])/ndt)
            divt*=div[i]
            divz*=(2*vdiv[i]-1)

        tw,tw0,pw0,tw1,pw1,dt=predistrb(tmin,tmax,div,dim)
        tw1tmp=tw1
        flag=1
        if dim==1 and flag==0:
            xh,xmx_value,xmx_index,xh_value,xh_div,xlayer=self.bunpu_contour(vdiv[0])#
            zh,zmx_value,zmx_index,zh_value,zh_div,zlayer=Z.bunpu_contour(vdiv[0])
        else:
            xmx_value2,xmx_index2,xdiv_flatten,xdiv_value,xdiv_flattens,xdiv_valus=self.bunpu_contour2(vdiv[0])
            zmx_value2,zmx_index2,zdiv_flatten,zdiv_value,zdiv_flattens,zdiv_valus=Z.bunpu_contour2(vdiv[0])

        nmin=0
        nmax=divt-1

        if dim==1:
            if flag==0:
                zptmp=np.zeros(divt)
                ztmp=np.zeros(divt)
                ztmp[xmx_index]=(X1[0][xmx_index]*(ndt-1)+Z1[0][zmx_index])/ndt
                zptmp[xmx_index]=(X1[1][xmx_index]*(ndt-1)+Z1[1][zmx_index])/ndt
            else:
                zptmp=np.zeros(divz)
                ztmp=[np.zeros(divz)]
        elif dim==2:
            zptmp=np.zeros(divz)
            ztmp=[np.zeros(divz),np.zeros(divz)]#
        elif dim==3:
            zptmp=np.zeros(divz)
            ztmp=[np.zeros(divz),np.zeros(divz),np.zeros(divz)]#
        ztmptest=[]
        
        zptmptest=[]
        if dim==1 and flag==0:
            for i in range(vdiv[0]+1):
            #
            
                nz0=len(zh_div[0][i])
                nz0h=max(3,divt*0.1)
                if nz0==0:
                    zin0tmp=zh_div[0][i-1][len(zh_div[0][i-1])-1]
                    zindx0=[zin0tmp,zin0tmp]
                    zh_div[0][i].append(zin0tmp)
                elif nz0==1:
                    zin0tmp=zh_div[0][i][0]
                    zindx0=[zin0tmp,zin0tmp]
                elif nz0>=2:
                    zindx0=[zh_div[0][i][0],zh_div[0][i][len(zh_div[0][i])-1]]#
                nx0=len(xh_div[0][i])
                if nz0<nz0h:
                    zpos0=np.linspace(Z1[0][zindx0[0]],Z1[0][zindx0[1]],nx0)
                    zpp0=np.linspace(Z1[1][zindx0[0]],Z1[1][zindx0[1]],nx0)
                elif nz0>=nz0h:
                    zdivt=[Z1[0][j] for j in zh_div[0][i]]
                    zdivp=[Z1[1][j] for j in zh_div[0][i]]
                    zfun_smooth_rbf = interpolate.Rbf(zdivt, zdivp, function='linear', smooth=0)
                    zpos0=np.linspace(Z1[0][zindx0[0]],Z1[0][zindx0[1]],nx0)
                    zpp0=[zfun_smooth_rbf(j) for j in zpos0]
                for j in range(nx0):
                    if dim==1:
                        ztmp[xh_div[0][i][j]]=(X1[0][xh_div[0][i][j]]*(ndt-1)+zpos0[j])/ndt
                        zptmp[xh_div[0][i][j]]=(X1[1][xh_div[0][i][j]]*(ndt-1)+zpp0[j])/ndt
                #
                nz1=len(zh_div[1][i])
                nz1h=max(3,divt*0.1)
                if nz1==0:
                    zin1tmp=zh_div[1][i-1][0]
                    zindx1=[zin1tmp,zin1tmp]
                    zh_div[1][i].append(zin1tmp)
                elif nz1==1:
                    zin1tmp=zh_div[1][i][0]
                    zindx1=[zin1tmp,zin1tmp]
                elif nz1>=2:
                    zindx1=[zh_div[1][i][0],zh_div[1][i][len(zh_div[1][i])-1]]#
                nx1=len(xh_div[1][i])
                if nz1<nz1h:
                    zpos1=np.linspace(Z1[0][zindx1[0]],Z1[0][zindx1[1]],nx1)#
                    zpp1=np.linspace(Z1[1][zindx1[0]],Z1[1][zindx1[1]],nx1)
                elif nz1>=nz1h:
                    zdivt=[Z1[0][j] for j in zh_div[1][i]]
                    zdivp=[Z1[1][j] for j in zh_div[1][i]]
                    zfun_smooth_rbf = interpolate.Rbf(zdivt, zdivp, function='linear', smooth=0)
                    zpos1=np.linspace(Z1[0][zindx1[0]],Z1[0][zindx1[1]],nx1)
                    zpp1=[zfun_smooth_rbf(j) for j in zpos1]
                    
                for j in range(nx1):
                    if dim==1:
                        ztmp[xh_div[1][i][j]]=(X1[0][xh_div[1][i][j]]*(ndt-1)+zpos1[j])/ndt
                        zptmp[xh_div[1][i][j]]=(X1[1][xh_div[1][i][j]]*(ndt-1)+zpp1[j])/ndt
        else:
            for k in range(divz):
                for j in range(dim):
                    ztmp[j][k]=(xdiv_flattens[j][k]*(ndt-1)+zdiv_flattens[j][k])/ndt
                zptmp[k]=(xdiv_valus[k]*(ndt-1)+zdiv_valus[k])/ndt
        #
        #
        if dim==1:
            zfun_smooth_rbf = interpolate.Rbf(ztmp, zptmp, function='linear', smooth=0)
        elif dim==2:
            zptmpmsh0=zptmp.reshape(2*vdiv[0]-1,2*vdiv[1]-1)
            zptmpmsh=np.array(zptmpmsh0)
            ztmpmsh=[]
            for i in range(dim):
                ztmpmsh0=(ztmp[i].reshape(2*vdiv[0]-1,2*vdiv[1]-1)).tolist()
                ztmpmsh.append(np.array(ztmpmsh0))
            zfun_smooth_rbf = interpolate.Rbf(ztmpmsh[0],ztmpmsh[1], zptmpmsh, function='linear', smooth=0)#

        elif dim==3:
            zptmpmsh=zptmp.reshape(2*vdiv[0]-1,2*vdiv[1]-1,2*vdiv[2]-1)
            ztmpmsh=[]
            for i in range(dim):
                ztmpmsh.append(ztmp[i].reshape(2*vdiv[0]-1,2*vdiv[1]-1,2*vdiv[2]-1))
            zfun_smooth_rbf = interpolate.Rbf(ztmpmsh[0],ztmpmsh[1],ztmpmsh[2], zptmpmsh, function='linear', smooth=0)
        menseki=0

        for i in range(divt):#
            if dim==2:
                tq,tmod=divmod(i,div[0])
            if dim==3:
                tq,tmod=divmod(i,div[0])
                tq1,tq2=divmod(tq,div[1])
            if dim==1:
                pw0[i]=zfun_smooth_rbf(tw0[i])
                pw1[i]=zfun_smooth_rbf(tw1[i])
                menseki+=pw0[i]*dt[0]
            elif dim==2:
                pw1[i]=zfun_smooth_rbf(tw1[0][i],tw1[1][i])
                pw0[tq][tmod]=pw1[i]
                menseki+=pw1[i]*dt[0]*dt[1]
            elif dim==3:
                pw1[i]=zfun_smooth_rbf(tw1[0][i],tw1[1][i],tw1[2][i])
                pw0[tq1][tq2][tmod]=pw1[i]
                menseki+=pw1[i]*dt[0]*dt[1]*dt[2]
        W=bunpu()
        if shw==0:
            print('simu:menseki',menseki)
        if dim==1:
            W.para=[tw]
            W.flatten=[tw1,pw1/menseki]
            W.mesh=[tw0,pw0/menseki]
        elif dim==2:
            W.para=[tw[0],tw[1]]
            W.flatten=[tw1[0],tw1[1],pw1/menseki]
            W.mesh=[tw0[0],tw0[1],pw0/menseki]
        elif dim==3:
            W.para=[tw[0],tw[1],tw[2]]
            W.flatten=[tw1[0],tw1[1],tw1[2],pw1/menseki]
            W.mesh=[tw0[0],tw0[1],tw0[2],pw0/menseki]
        W.cmap=Z.cmap
        W.dim=dim
        W.xmin=tmin
        W.xmax=tmax
        W.div=div
        W.dx=dt
        W.pmax=(pw1/menseki).max()
        return W

    def bunpu_product(self,other,divt=[],divs=[],corel=[0],shw=0):
        filename='x・y_log'
        logname=filename+'.csv'
        log = open(logname,'w' , encoding='shift_jis' )
        info=['tパラメータ','p値']
        writer = csv.writer(log, lineterminator='\n')
        writer.writerow(info)
        x2=self.flatten
        x20=self.mesh
        x,x2,x20,xmin,xmax,dim,dx,divx,m,y0min,y0max,divt0,flag0=self.bunpu_filter(other,1,divt)#
        if flag0!=0:
            if flag0!=3:
                Y=m
                if divt==[]:
                    divt=divx
                if dim==1:
                    y0=np.linspace(y0min[0],y0max[0],len(x2[0]))
                    t10=[x2[0][i]*y0[i] for i in range(len(x2[0]))]
                    t00=t10
                    tmin=[xmin[0]*y0min[0]]
                    tmax=[xmax[0]*y0max[0]]
                    p1=x2[1]
                    p0=x2[1]
                    divtt=divt[0]
                    zfun_smooth_rbf = interpolate.Rbf(t00, x20[1], function='linear', smooth=0)#
                elif dim==2:

                    y0=[np.linspace(y0min[0],y0max[0],len(x2[0])),np.linspace(y0min[0],y0max[0],len(x2[1]))]
                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]*y0[0][i] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]*y0[1][i] for i in range(len(x2[1]))])
                    t00.append([[x20[0][i][j]*y0[0][j] for j in range(len(x20[0,:]))] for i in range(len(x20[:,0]))])
                    t00.append([[x20[1][i][j]*y0[1][i] for j in range(len(x20[0,:]))] for i in range(len(x20[:,0]))])
                    
                    tmin.append(xmin[0]*y0min[0])
                    tmin.append(xmin[1]*y0min[0])
                    tmax.append(xmax[0]*y0max[0])
                    tmax.append(xmax[1]*y0max[0])
                    p1=x2[2]
                    p0=x20[2]
                    divtt=divt[0]*divt[1]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], x20[2], function='linear', smooth=0)#

                elif dim==3:
                    y0=[np.linspace(y0min[0],y0max[0],len(x2[0])),np.linspace(y0min[0],y0max[0],len(x2[1])),np.linspace(y0min[0],y0max[0],len(x2[2]))]
                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]*y0[0][i] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]*y0[1][i] for i in range(len(x2[1]))])
                    t10.append([x2[2][i]*y0[2][i] for i in range(len(x2[2]))])
                    t00.append([[[x20[0][i][j][k]*y0[0][j] for k in range(len(x20[0][0,0,:]))] for j in range(len(x20[0][0,:,0]))] for i in range(len(x20[0][:,0,0]))])
                    t00.append([[[x20[1][i][j][k]*y0[1][i] for k in range(len(x20[1][0,0,:]))] for j in range(len(x20[1][0,:,0]))] for i in range(len(x20[1][:,0,0]))])
                    t00.append([[[x20[2][i][j][k]*y0[2][k] for k in range(len(x20[2][0,0,:]))] for j in range(len(x20[2][0,:,0]))] for i in range(len(x20[2][:,0,0]))])
                    tmin.append(xmin[0]*y0min[0])
                    tmin.append(xmin[1]*y0min[0])
                    tmin.append(xmin[2]*y0min[0])
                    tmax.append(xmax[0]*y0max[0])
                    tmax.append(xmax[1]*y0max[0])
                    tmax.append(xmax[2]*y0max[0])
                    p1=x2[3]
                    p0=x20[3]
                    divtt=divt[0]*divt[1]*divt[2]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], t00[2], x20[3], function='linear', smooth=0)#
                tt=[]
                dt=[]
                for i in range(dim):
                    tt.append(np.linspace(tmin[i],tmax[i],divx[i]))
                    dt.append((tmax[i]-tmin[i])/(divt[i]-1))
                menseki=0
                if dim == 1:
                    t=[t00]
                    t1=t
                    t0=t
                    p0=np.zeros(divt[0])
                    p1=np.zeros(divt[0])
                elif dim == 2:
                    t=[tt[0],tt[1]]
                    ta,tb=np.meshgrid(tt[0],tt[1])
                    t1=[ta.flatten(),tb.flatten()]
                    t0=[ta,tb]
                    p0=np.zeros((divt[1],divt[0]))
                    p1=np.zeros(divtt)
                elif dim == 3:
                    t=[tt[0],tt[1],tt[2]]
                    ta,tb,tc=np.meshgrid(tt[0],tt[1],tt[2])
                    t1=[ta.flatten(),tb.flatten(),tc.flatten()]
                    t0=[ta,tb,tc]
                    p0=np.zeros((divt[2],divt[1],divt[0]))
                    p1=np.zeros(divtt)
                    
                for ti in range(divtt):
                    if dim==1:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti])*(xmax[0]-xmin[0])/(tmax[0]-tmin[0])
                        p0[ti]=p1[ti]
                        menseki+=p1[ti]*dt[0]
                    elif dim==2:
                        
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1]))
                        tq,tmod=divmod(ti,divt[0])

                        p0[tq][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1])
                    elif dim==3:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti],t1[2][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])\
                            /((tmax[0]-tmin[0])*(tmax[1]-tmin[1])*(tmax[2]-tmin[2]))
                        tq,tmod=divmod(ti,divt[0])
                        tq1,tq2=divmod(tq,divt[1])
                        p0[tq1][tq2][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
            #
            else:#
                Y=other
                if divt==[]:
                    divt=divx
                if dim==1:
                
                    t10=[x2[0][i]*Y[0] for i in range(len(x2[0]))]
                    t00=[x2[0][i]*Y[0] for i in range(len(x2[0]))]
                    p1=x2[1]
                    p0=x2[1]
                    tmin=[xmin[0]*Y[0]]
                    tmax=[xmax[0]*Y[0]]
                    divtt=divt[0]
                    zfun_smooth_rbf = interpolate.Rbf(t00, x20[1], function='linear', smooth=0)#
                elif dim==2:

                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]*Y[0] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]*Y[0] for i in range(len(x2[1]))])
                    t00.append([[x20[0][i][j]*Y[0] for j in range(len(x20[0][0]))] for i in range(len(x20[0]))])
                    t00.append([[x20[1][i][j]*Y[0] for j in range(len(x20[1][0]))] for i in range(len(x20[1]))])
                    
                    tmin.append(xmin[0]*Y[0])
                    tmin.append(xmin[1]*Y[0])
                    tmax.append(xmax[0]*Y[0])
                    tmax.append(xmax[1]*Y[0])
                    p1=x2[2]
                    p0=x20[2]
                    divtt=divt[0]*divt[1]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], x20[2], function='linear', smooth=0)#

                elif dim==3:
                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]*Y[0] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]*Y[0] for i in range(len(x2[1]))])
                    t10.append([x2[2][i]*Y[0] for i in range(len(x2[2]))])
                    t00.append([[[x20[0][i][j][k]*Y[0] for k in range(len(x20[0][0][0]))] for j in range(len(x20[0][0]))] for i in range(len(x20[0]))])
                    t00.append([[[x20[1][i][j][k]*Y[0] for k in range(len(x20[1][0][0]))] for j in range(len(x20[1][0]))] for i in range(len(x20[1]))])
                    t00.append([[[x20[2][i][j][k]*Y[0] for k in range(len(x20[2][0][0]))] for j in range(len(x20[2][0]))] for i in range(len(x20[2]))])
                    tmin.append(xmin[0]*Y[0])
                    tmin.append(xmin[1]*Y[0])
                    tmin.append(xmin[2]*Y[0])
                    tmax.append(xmax[0]*Y[0])
                    tmax.append(xmax[1]*Y[0])
                    tmax.append(xmax[2]*Y[0])
                    p1=x2[3]
                    p0=x20[3]
                    divtt=divt[0]*divt[1]*divt[2]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], t00[2], x20[3], function='linear', smooth=0)#
                tt=[]
                dt=[]
                for i in range(dim):
                    tt.append(np.linspace(tmin[i],tmax[i],divx[i]))
                    dt.append((tmax[i]-tmin[i])/(divt[i]-1))
                menseki=0
                if dim == 1:
                    t=[t00]
                    t1=t
                    t0=t
                    p0=np.zeros(divt[0])
                    p1=np.zeros(divt[0])
                elif dim == 2:
                    t=[tt[0],tt[1]]
                    ta,tb=np.meshgrid(tt[0],tt[1])
                    t1=[ta.flatten(),tb.flatten()]
                    t0=[ta,tb]
                    p0=np.zeros((divt[1],divt[0]))
                    p1=np.zeros(divtt)
                elif dim == 3:
                    t=[tt[0],tt[1],tt[2]]
                    ta,tb,tc=np.meshgrid(tt[0],tt[1],tt[2])
                    t1=[ta.flatten(),tb.flatten(),tc.flatten()]
                    t0=[ta,tb,tc]
                    p0=np.zeros((divt[2],divt[1],divt[0]))
                    p1=np.zeros(divtt)
                    
                for ti in range(divtt):
                    if dim==1:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti])*(xmax[0]-xmin[0])/(tmax[0]-tmin[0])
                        p0[ti]=p1[ti]
                        menseki+=p1[ti]*dt[0]
                    elif dim==2:
                        
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1]))
                        tq,tmod=divmod(ti,divt[0])
                        p0[tq][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1])
                    elif dim==3:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti],t1[2][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])\
                            /((tmax[0]-tmin[0])*(tmax[1]-tmin[1])*(tmax[2]-tmin[2]))
                        tq,tmod=divmod(ti,divt[0])
                        tq1,tq2=divmod(tq,divt[1])
                        p0[tq1][tq2][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
            #
        else:
            y2=other.flatten
            if other.xmin==[] or other.xmax==[] or other.dim==0 or other.div==[]:
                yf,py0,y0,dimy,divy,dy,ymin,ymax=other.kaiseki()
                other.xmin=ymin
                other.xmax=ymax
                other.dim=dimy
                other.div=divy
                other.dx=dy
            else:
                ymin=other.xmin
                ymax=other.xmax
                dimy=other.dim
                divy=(other.div)[0]
                dy=other.dx
            if dimy==2:
                print('no anser')
                return None
            divxt=1
            for i in range(dim):
                divxt*=divx[i]
    
            corner1=np.array(xmin)*ymin[0]
            corner2=np.array(xmin)*ymax[0]
            corner3=np.array(xmax)*ymin[0]
            corner4=np.array(xmax)*ymax[0]
            tmin=[]
            tmax=[]
            dt=[]
            da=[]
            if divt==[]:
                divt=divx
            if divs==[]:
                divs=copy.copy(divt)
            for i in range(dim):
                tmin.append(min(corner1[i],corner2[i],corner3[i],corner4[i]))
                tmax.append(max(corner1[i],corner2[i],corner3[i],corner4[i]))
                dt.append((tmax[i]-tmin[i])/(divt[i]))
            if dim==1:
                divtt=divt[0]
                t=[np.linspace(tmin[0],tmax[0],divt[0])]
                t0=t
                t1=t
                p1=np.zeros(divt[0])
                lastp=np.zeros(divt[0])
                base_x = interpolate.Rbf(x2[0], x2[1], function='linear', smooth=0)  #
                base_y = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0)  #
                da=(ymax[0]-ymin[0]+xmax[0]-xmin[0])/(divs[0]-1)
                signt=[0 if i < divtt-1 and np.sign(t1[0][i])!=np.sign(t1[0][i+1]) else np.sign(t1[0][i]) for i in range(divtt)]
            elif dim ==2:
                divtt=divt[0]*divt[1]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1])]
                tx,ty=np.meshgrid(t[0],t[1])#
                t1=[tx.flatten(),ty.flatten()]
                xx,xy=np.meshgrid(np.linspace(xmin[0],xmax[0],divt[0]),np.linspace(xmin[1],xmax[1],divt[1]))#
                t0=[tx,ty]
                p1=np.zeros(divtt)
                p0=np.zeros((divt[1],divt[0]))#7/2
                lastp=np.zeros(divtt)
                base_x = interpolate.Rbf(x20[0], x20[1], x20[2], function='linear', smooth=0)  #
                base_y = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0) 
                
                
                
            elif dim==3:
                divtt=divt[0]*divt[1]*divt[2]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1]),np.linspace(tmin[2],tmax[2],divt[2])]
                tx,ty,tz=np.meshgrid(t[0],t[1],t[2])#
                t1=[tx.flatten(),ty.flatten(),tz.flatten()]
                xx,xy,xz=np.meshgrid(np.linspace(xmin[0],xmax[0],divt[0]),np.linspace(xmin[1],xmax[1],divt[1]),np.linspace(xmin[2],xmax[2],divt[2]))#
                t0=[tx,ty,tz]
                p1=[0.]*divtt
                p0=np.zeros((divt[2],divt[1],divt[0]))
                lastp=np.zeros(divtt)
                base_x = interpolate.Rbf(x20[0], x20[1], x20[2], x20[3], function='linear', smooth=0)  #
                base_y = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0)  #
            menseki=0
            test1x=[[],[]]
            test1y=[[],[]]
            test2x=[[],[]]
            test2y=[[],[]]
            for ti in range(divtt):#
                max_x=0
                max_y=0
                max_p=0
                if dim==1:
                    dts=copy.copy(dt[0])
                    max_s=xmax[0]
                    min_s=xmin[0]
                    min_t=tmin[0]
                    max_t=tmax[0]
                    tti=ti
                if dim==2:
                    
                    fg=0
                    tq,tmod=divmod(ti,divt[0])

                    divst, min_t, max_t, min_s, max_s, min_sx, max_sx, min_sy, max_sy, plmin,\
                        plmax=slice2d_bunpu(ti,t1,tmin,tmax,dt,divs,xmin,xmax,ymin,ymax)
                    if plmin!=9 and plmax!=9  and max_t>min_t:
                        if divst<=1:
                            divst=1
                        dts=(max_t-min_t)/divst
                        tspan=(t1[0][ti]**2+t1[1][ti]**2)**0.5
                        tti=(tspan-min_t)//dts#

                        divs[0]=divst
                        divss=int(divs[0])
                        s=np.linspace(min_s,max_s,divss)
                        sx=np.linspace(min_sx,max_sx,divss)
                        sy=np.linspace(min_sy,max_sy,divss)
                        if divss<=1:
                            ds=(max_s-min_s)
                        else:
                            ds=(max_s-min_s)/(divss-1)
                        ps=[]
                        for i in(range(divss)):
                            ps.append(base_x(sx[i],sy[i]))
                        divs[0]=divss
                        if divs[0]<=1:
                            da=(ymax[0]-ymin[0]+max_s-min_s)
                        else:
                            da=(ymax[0]-ymin[0]+max_s-min_s)/(divs[0]-1)#
                        if len(s)>1 and max_s>min_s:
                            base_s = interpolate.Rbf(s, ps, function='linear', smooth=0)
                        else:
                            base_s=interpolate.Rbf([s[0]-ds,s[0],s[0]+ds], [ps[0],ps[0],ps[0]], function='linear', smooth=0)
                        signt=[0 if i < divss-1 and np.sign(s[i])!=np.sign(s[i+1]) else np.sign(s[i]) for i in range(divss)]#
                    else:
                        divs[0]=0
                if dim==3:
                    tq,tmod=divmod(ti,divt[0])
                    tq1,tq2=divmod(tq,divt[1])
                    divst, min_t, max_t, min_s, max_s, min_sx, max_sx, min_sy, max_sy, min_sz, max_sz, \
                        plmin, plmax=slice3d_bunpu(ti,t1,tmin,tmax,dt,divs,xmin,xmax,ymin,ymax)
                    if plmin!=9 and plmax!=9 and abs(max_s-min_s)>min((xmax[0]-xmin[0])/divx[0],(xmax[1]-xmin[1])/divx[1],(xmax[2]-xmin[2])/divx[2]):# and max_t>min_t:
                        
                        if divst<=1:
                            divst=1
                        dts=(max_t-min_t)/divst
                        tspan=(t1[0][ti]**2+t1[1][ti]**2+t1[2][ti]**2)**0.5

                        tti=-1*(-1*(tspan-min_t)//dts)
                        divs[0]=divst
                        divss=int(divs[0])
                        s=np.linspace(min_s,max_s,divss)
                        sx=np.linspace(min_sx,max_sx,divss)
                        sy=np.linspace(min_sy,max_sy,divss)
                        sz=np.linspace(min_sz,max_sz,divss)
                        if divss<=1:
                            ds=(max_s-min_s)
                        else:
                            ds=(max_s-min_s)/(divss-1)
                        ps=[]
                        for i in(range(divss)):
                            ps.append(max(0,base_x(sx[i],sy[i],sz[i])))
                        divs[0]=divss
                        if divs[0]<=1:
                            da=(ymax[0]-ymin[0]+max_s-min_s)
                        else:
                            da=(ymax[0]-ymin[0]+max_s-min_s)/(divs[0]-1)#
                        if len(s)>1 and max_s>min_s:
                            base_s = interpolate.Rbf(s, ps, function='linear', smooth=0)

                        else:

                            base_s=interpolate.Rbf([s[0]-ds,s[0],s[0]+ds], [ps[0],ps[0],ps[0]], function='linear', smooth=0)
                        signt=[0 if i < divss-1 and np.sign(s[i])!=np.sign(s[i+1]) else np.sign(s[i]) for i in range(divss)]
                    else:
                        divs[0]=0

                #
                test1x=[]
                test1y=[]
                for js in range(divs[0]):
                    xh=np.zeros(dim)
                    i=0

                    if signt[js]>0:
                        if (ymax[i]-min_s-js*da)**2+4*(min_t+tti*dts)>=0 and (ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts)>=0 \
                           and (ymax[i]-min_s-(js+1)*da)**2+4*(min_t+tti*dts)>=0 and (ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts)>=0:
                            xij=((-(ymax[i]-min_s-js*da)+((ymax[i]-min_s-js*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            xi1j=((-(ymax[i]-min_s-js*da)+((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)
                            xij1=((-(ymax[i]-min_s-(js+1)*da)+((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            yij=(((ymax[i]-min_s-js*da)+((ymax[i]-min_s-js*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            xi1j1=((-(ymax[i]-min_s-(js+1)*da)+((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#
                            yij1=(((ymax[i]-min_s-(js+1)*da)+((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            yi1j=(((ymax[i]-min_s-js*da)+((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)
                            yi1j1=(((ymax[i]-min_s-(js+1)*da)+((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)
                            #対
                            xijp=((-(ymax[i]-min_s-js*da)-((ymax[i]-min_s-js*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            xi1jp=((-(ymax[i]-min_s-js*da)-((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)
                            xij1p=((-(ymax[i]-min_s-(js+1)*da)-((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            yijp=(((ymax[i]-min_s-js*da)-((ymax[i]-min_s-js*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            xi1j1p=((-(ymax[i]-min_s-(js+1)*da)-((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)
                            yij1p=(((ymax[i]-min_s-(js+1)*da)-((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+tti*dts))**0.5)/2)
                            yi1jp=(((ymax[i]-min_s-js*da)-((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)
                            yi1j1p=(((ymax[i]-min_s-(js+1)*da)-((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)
                            spxij=max(xi1j1,xij)-min(xi1j1,xij)
                            spyij=max(yi1j,yij1)-min(yi1j,yij1)
                            spxijp=max(xi1j1p,xijp)-min(xi1j1p,xijp)
                            spyijp=max(yi1jp,yij1p)-min(yi1jp,yij1p)
                            xa=xij
                            ya=yij
                            xb=xijp
                            yb=yijp
                            dareaa=(spxij*spyij)/2#
                            dareab=(spxijp*spyijp)/2
                            flag2a=cor2ran(xij,min_s,max_s,yij,ymin[i],ymax[i],corel[i])
                            flag2b=cor2ran(xijp,min_s,max_s,yijp,ymin[i],ymax[i],corel[i])
                            flag2=flag2a or flag2b
                            darea=flag2a*dareaa+flag2b*dareab
                        else:
                            flag2=0

                    elif signt[js]<=0:
                        if (ymax[i]+max_s-js*da)**2-4*(min_t+(tti-1)*dts)>=0 and (ymax[i]+max_s-js*da)**2-4*(min_t+(tti)*dts)>=0 \
                           and (ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti)*dts)>=0 and (ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti-1)*dts)>=0:
                            xij=(((ymax[i]+max_s-js*da)+((ymax[i]+max_s-js*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#js==16
                            xi1j=(((ymax[i]+max_s-js*da)+((ymax[i]+max_s-js*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#16
                            xij1=(((ymax[i]+max_s-(js+1)*da)+((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#15
                            yij=(((ymax[i]+max_s-js*da)+((ymax[i]+max_s-js*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#16
                            xi1j1=(((ymax[i]+max_s-(js+1)*da)+((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#15
                            yij1=(((ymax[i]+max_s-(js+1)*da)+((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#15
                            yi1j=(((ymax[i]+max_s-js*da)+((ymax[i]+max_s-js*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#16
                            yi1j1=(((ymax[i]+max_s-(js+1)*da)+((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#15
                            #対
                            xijp=(((ymax[i]+max_s-js*da)-((ymax[i]+max_s-js*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#16
                            xi1jp=(((ymax[i]+max_s-js*da)-((ymax[i]+max_s-js*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#16
                            xij1p=(((ymax[i]+max_s-(js+1)*da)-((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#15
                            yijp=(((ymax[i]+max_s-js*da)-((ymax[i]+max_s-js*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#16
                            xi1j1p=(((ymax[i]+max_s-(js+1)*da)-((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#15
                            yij1p=(((ymax[i]+max_s-(js+1)*da)-((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti-1)*dts))**0.5)/2)#15
                            yi1jp=(((ymax[i]+max_s-js*da)-((ymax[i]+max_s-js*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#16
                            yi1j1p=(((ymax[i]+max_s-(js+1)*da)-((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+(tti)*dts))**0.5)/2)#15
                            spxija=max(min(xi1j,xi1jp),min(xij1,xij1p))-min(min(xi1j,xi1jp),min(xij1,xij1p))#
                            spyija=max(max(yij,yijp),max(yi1j1,yi1j1p))-min(max(yij,yijp),max(yi1j1,yi1j1p))
                            spxijb=max(max(xij,xijp),max(xi1j1,xi1j1p))-min(max(xij,xijp),max(xi1j1,xi1j1p))#
                            spyijb=max(min(yi1j,yi1jp),min(yij1,yij1p))-min(min(yi1j,yi1jp),min(yij1,yij1p))
                            xa=min(xij,xijp)
                            ya=max(yij,yijp)
                            xb=max(xij,xijp)
                            yb=min(yij,yijp)
                            dareaa=(spxija*spyija)/2
                            dareab=(spxijb*spyijb)/2
                            flag2a=cor2ran(min(xi1j,xi1jp),min_s,max_s,max(yij,yijp),ymin[i],ymax[i],corel[i])
                            flag2b=cor2ran(max(xi1j,xi1jp),min_s,max_s,min(yij,yijp),ymin[i],ymax[i],corel[i])
                            flag2=flag2a or flag2b
                            darea=flag2a*dareaa+flag2b*dareab
                        else:
                            flag2=0

                    else:

                        xijs=(((ymax[i]+max_s-js*da)+((ymax[i]+max_s-js*da)**2-4*(min_t+tti*dts))**0.5)/2)#
                        xij1s=(((ymax[i]+max_s-(js+1)*da)+((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+tti*dts))**0.5)/2)
                        yijs=(((ymax[i]+max_s-js*da)+((ymax[i]+max_s-js*da)**2-4*(min_t+tti*dts))**0.5)/2)#
                        yij1s=(((ymax[i]+max_s-(js+1)*da)+((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+tti*dts))**0.5)/2)

                        xijsp=(((ymax[i]+max_s-js*da)-((ymax[i]+max_s-js*da)**2-4*(min_t+tti*dts))**0.5)/2)#
                        xij1sp=(((ymax[i]+max_s-(js+1)*da)-((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+tti*dts))**0.5)/2)
                        yijsp=(((ymax[i]+max_s-js*da)-((ymax[i]+max_s-js*da)**2-4*(min_t+tti*dts))**0.5)/2)#
                        yij1sp=(((ymax[i]+max_s-(js+1)*da)-((ymax[i]+max_s-(js+1)*da)**2-4*(min_t+tti*dts))**0.5)/2)

                        xi1jt=((-(ymax[i]-min_s-js*da)+((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#16
                        xi1j1t=((-(ymax[i]-min_s-(js+1)*da)+((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#15
                        yi1jt=(((ymax[i]-min_s-js*da)+((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#16
                        yi1j1t=(((ymax[i]-min_s-(js+1)*da)+((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#15
                        #対
                        xi1jtp=((-(ymax[i]-min_s-js*da)-((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#16
                        xi1j1tp=((-(ymax[i]-min_s-(js+1)*da)-((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#15
                        yi1jtp=(((ymax[i]-min_s-js*da)-((ymax[i]-min_s-js*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#16
                        yi1j1tp=(((ymax[i]-min_s-(js+1)*da)-((ymax[i]-min_s-(js+1)*da)**2+4*(min_t+(tti+1)*dts))**0.5)/2)#15
                       
                        spxija=max(max(xi1j1t,xi1j1tp),min(xij1s,xij1sp))-min(max(xi1j1t,xi1j1tp),min(xij1s,xij1sp))#左上
                        spyija=max(max(yi1jt,yi1jtp),max(yijs,yijsp))-min(max(yi1j1t,yi1j1tp),max(yij1s,yij1sp))
                        spxijb=max(max(xijs,xijsp),min(xi1jt,xi1jtp))-min(max(xijs,xijsp),min(xi1jt,xi1jtp))#右下
                        spyija=max(min(yijs,yijsp),min(yi1jt,yi1jtp))-min(min(yij1s,yij1sp),min(yi1j1t,yi1j1tp))
                        xa=min(xijs,xijsp)
                        ya=max(yijs,yijsp)
                        xb=max(xijs,xijsp)
                        yb=min(yijs,yijsp)
                        dareaa=(spxija*spyija)#
                        dareab=(spxijb*spyijb)#
                        flag2a=cor2ran(min(xijs,xijsp),min_s,max_s,max(yijs,yijsp),ymin[i],ymax[i],corel[i])
                        flag2b=cor2ran(max(xi1jt,xi1jtp),min_s,max_s,min(yi1jt,yi1jtp),ymin[i],ymax[i],corel[i])
                        flag2=flag2a or flag2b
                        darea=flag2a*dareaa+flag2b*dareab
                    test1x.append(xij)
                    test1y.append(yij)

                    if flag2==1.0:
                        if dim==1:
                            pxa = base_x(xa)
                            pya = base_y(ya)
                            pxb = base_x(xb)
                            pyb = base_y(yb)
                            p1[ti]+=(pxa*pya*dareaa*flag2a+pxb*pyb*dareab*flag2b)/dt[0]
                            p0=p1
                            #
                        elif dim==2:
                            pxa = base_s(xa)
                            pya = base_y(ya)
                            pxb = base_s(xb)
                            pyb = base_y(yb)
                            pp=(xa*pxa*pya*dareaa*flag2a/tspan+xb*pxb*pyb*dareab*flag2b/tspan)/dts#

                            p1[ti]+=max(0,pp)
                            p0[tq][tmod]= p1[ti]
                            
                        elif dim==3:
                            pxa = max(0,base_s(xa))
                            pya = base_y(ya)
                            pxb = max(0,base_s(xb))
                            pyb = base_y(yb)
                            p1[ti]+=((xa/tspan)**2*pxa*pya*dareaa*flag2a+(xb/tspan)**2*pxb*pyb*dareab*flag2b)/dts
                            p0[tq1][tq1][tmod]= p1[ti]
                        else:
                            print('error')
                if dim==1:
                    menseki+=p1[ti]*dt[0]
                if dim==2:
                    menseki+=p1[ti]*(dt[0]*dt[1])
                elif dim==3:
                    menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
                gp=0
                if gp==1:
                    fig=plt.figure()
                    flmx=[min_s,min_s,max_s,max_s,min_s]
                    flmy=[ymin[0],ymax[0],ymax[0],ymin[0],ymin[0]]
                    plt.plot(flmx,flmy)
                    plt.plot(test1x[0],test1y[0])

                    plt.savefig("test"+str(ti))
        #
        if shw==0:
            flaginf=['bunpu*bunpu','bunpu*lean','lean*bunpu','bunpu*vector']
            print(flaginf[flag0],menseki)
        if menseki>0:
            if dim==1:
                p00=p0/menseki
                p=p1/menseki
            else:
                p=p1/menseki
                p00=p0/menseki

        else:
            if dim==1:
                p00=p0
                p=p1
            else:
                p=p1
                p00=p0
        #
        if dim==1:
            W=[t[0]]
            W0=[t0[0],p00]
            W1=[t1[0],p]
        elif dim==2:
            W=[list(t[0]),list(t[1])]
            W0=[t0[0],t0[1],p00]
            W1=[t1[0],t1[1],p]
        elif dim==3:
            W=[t[0],t[1],t[2]]
            W0=[t0[0],t0[1],t0[2],p00]
            W1=[t1[0],t1[1],t1[2],p]
        WW=bunpu()
        WW.para=W
        WW.mesh=W0
        WW.flatten=W1
        WW.div=divt
        WW.dx=dt
        WW.dim=dim
        WW.xmin=tmin
        WW.xmax=tmax
        WW.pmax=p.max()
        return WW
        
    def bunpu_multiplier(self,multifact=2,divt0=0,neg=0,shw=0):
        xpara=self.para
        xmesh=self.mesh
        xfltn=self.flatten
        dim=self.dim
        xmax=self.xmax
        xmin=self.xmin
        divx=self.div
        dx=self.dx
        if divt0!=0:
            divt=divt0
        else:
            divt=divx
        tmax=[]
        tmin=[]
        t=[]
        divtt=1
        dt=[]
        for i in range(dim):
            if xmax[i]>=0 or multifact%2!=0 or neg==0:
                tmax.append(xmax[i]**multifact)
            else:
                tmax.append(-1*xmax[i]**multifact)
            if xmin[i]>=0 or multifact%2!=0 or neg==0:
                tmin.append(xmin[i]**multifact)
            else:
                tmin.append(-1*xmin[i]**multifact)
            t.append(np.linspace(tmin[i],tmax[i],divt[i]))
            dt.append((tmax[i]-tmin[i])/divt[i])
            divtt*=divt[i]
        menseki=0
        if dim==1:
            x0=[]
            for i in range(divx[0]):
                if xmesh[0][i]>=0 or multifact%2!=0 or neg==0:
                    x0.append(xmesh[0][i]**multifact)
                else:
                    x0.append(-1*xmesh[0][i]**multifact)
            p=np.zeros(divx[0])
            for i in range(divx[0]):
                if i==0:
                    p[i]=xmesh[1][i]*dx[0]/(x0[i]-(xmesh[0][i]-dx)**multifact)
                else:
                    p[i]=xmesh[1][i]*dx[0]/(x0[i]-x0[i-1])
            pintplt = interpolate.Rbf(x0, p, function='linear', smooth=0)
            p00=np.zeros(divt[0])
            for i in range(divt[0]):
                p00[i]=pintplt(t[0][i])
                menseki+=p00[i]*dt[0]
            if menseki>0:
                p0=p00/menseki
            else:
                p0=p00
                
            p1=p0
        elif dim==2:
            x0=[]
            x1=[]
            for i in range(divx[0]):
                if xpara[0][i]>=0 or multifact%2!=0 or neg==0:
                    x0.append(xpara[0][i]**multifact)
                else:
                    x0.append(-1*xpara[0][i]**multifact)
                
            for i in range(divx[1]):
                if xpara[1][i]>=0 or multifact%2!=0 or neg==0:
                    x1.append(xpara[1][i]**multifact)
                else:
                    x1.append(-1*xpara[1][i]**multifact)
                
            xx,yy=np.meshgrid(x0,x1)
            p=np.zeros((divx[1],divx[0]))
            for i in range(divx[1]):
                for j in range(divx[0]):
                    if i==0 and j==0:
                        p[i][j]=xmesh[2][i][j]*dx[0]*dx[1]/((x0[j+1]-x0[j])*(x1[i+1]-x1[i]))
                        if p[i][j]<0:
                            p[i][j]=0
                    elif i==0:

                        p[i][j]=xmesh[2][i][j]*dx[0]*dx[1]/((x1[i+1]-x1[i])*(x0[j]-x0[j-1]))
                        if p[i][j]<0:
                            p[i][j]=0
                    elif j==0:

                        p[i][j]=xmesh[2][i][j]*dx[0]*dx[1]/((x0[j+1]-x0[j])*(x1[i]-x1[i-1]))
                        if p[i][j]<0:
                            p[i][j]=0
                    else:
                        p[i][j]=xmesh[2][i][j]*dx[0]*dx[1]/((x0[j]-x0[j-1])*(x1[i]-x1[i-1]))
                        if p[i][j]<0:
                            p[i][j]=0

            pintplt = interpolate.Rbf(xx, yy, p, function='linear', smooth=0)
            p00=np.zeros((divt[1],divt[0]))
            pp=np.zeros((divt[1],divt[0]))
            t0=[[],[]]
            t0[0],t0[1]=np.meshgrid(t[0],t[1])
            t1=[t0[0].flatten(),t0[1].flatten()]
            for i in range(divt[1]):
                for j in range(divt[0]):

                    pp[i][j]=pintplt(t0[0][i][j],t0[1][i][j])
            haba=3
            hasi=[pp[0][0]*multifact,pp[0][divt[0]-1]**multifact,pp[divt[1]-1][divt[0]-1]**multifact,pp[divt[1]-1][0]*multifact]#
            for i in range(divt[1]):
                for j in range(divt[0]):
                    #
                    k=0
                    ph=0
                    if i<=haba-1:
                        ph=((divt[0]-1-j)*hasi[0]+j*hasi[1])/(divt[0]-1)
                        k=haba-i-1
                    elif j<=haba-1:
                        ph=((divt[1]-1-i)*hasi[0]+i*hasi[3])/(divt[1]-1)
                        k=haba-j-1
                    elif j>=divt[0]-haba:
                        ph=((divt[1]-1-i)*hasi[1]+i*hasi[2])/(divt[1]-1)
                        k=j+haba-divt[0]
                    elif i>=divt[1]-haba:
                        ph=((divt[0]-1-j)*hasi[3]+j*hasi[2])/(divt[0]-1)
                        k=i+haba-divt[1]
                    if k>=0 and ph>=0:
                        p00[i][j]=((haba-k)*pp[i][j]+k*ph)/haba
                    else:
                        p00[i][j]=pp[i][j]
                    #
                    
                    menseki+=p00[i][j]*dt[0]*dt[1]
            if menseki>0:
                p0=p00/menseki
            else:
                p0=p00
            p1=p0.flatten()
                    
            
        elif dim==3:
            x0=[]
            x1=[]
            x2=[]
            for i in range(divx[0]):
                if xpara[0][i]>=0 or multifact%2!=0 or neg==0:
                    x0.append(xpara[0][i]**multifact)
                else:
                    x0.append(-1*xpara[0][i]**multifact)
            for i in range(divx[1]):
                if xpara[1][i]>=0 or multifact%2!=0 or neg==0:
                    x1.append(xpara[1][i]**multifact)
                else:
                    x1.append(-1*xpara[1][i]**multifact)
            for i in range(divx[2]):
                if xpara[2][i]>=0 or multifact%2!=0 or neg==0:
                    x2.append(xpara[2][i]**multifact)
                else:
                    x2.append(-1*xpara[2][i]**multifact)
            xx,yy,zz=np.meshgrid(x0,x1,x2)
            p=np.zeros((divx[0],divx[1],divx[2]))
            for i in range(divx[2]):
                for j in range(divx[1]):
                    for k in range(divx[0]):
                        if i==0 and j==0 and k==0:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/\
                                (((x0[j]-xmesh[0][i][j][k]-dx[0])**multifact)*(x1[i]-(xmesh[1][i][j][k]-dx[1])**multifact)*(x2[i]-(xmesh[2][i][j][k]-dx[2])**multifact))
                                                                            
                            if p[i][j][k]<0:
                                p[i][j][k]=0
                        elif i==0 and j==0:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/\
                                ((x0[k]-x0[k-1])*(x1[i]-(xmesh[1][i][j][k]-dx[1])**multifact)*(x2[i]-(xmesh[2][i][j][k]-dx[2])**multifact))
                            if p[i][j][k]<0:
                                p[i][j][k]=0
                        elif i==0 and k==0:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/\
                                (((x0[j]-xmesh[0][i][j][k]-dx[0])**multifact)*(x1[j]-x1[j-1])*(x2[i]-(xmesh[2][i][j][k]-dx[2])**multifact))
                            if p[i][j][k]<0:
                                p[i][j][k]=0
                        elif j==0 and k==0:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/\
                                ((x0[j]-(xmesh[0][i][j][k]-dx[0])**multifact)*((x1[i]-xmesh[1][i][j][k]-dx[1])**multifact)*(x2[i]-x2[i-1]))
                            if p[i][j][k]<0:
                                p[i][j][k]=0
                        elif i==0:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/(((x2[i]-xmesh[2][i][j][k]-dx[2])**multifact)*(x1[j]-x1[j-1])*(x0[k]-x0[k-1]))
                            if p[i][j][k]<0:
                                p[i][j][k]=0
                        elif j==0:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/(((x1[j]-xmesh[1][i][j][k]-dx[1])**multifact)*(x2[i]-x2[i-1])*(x0[k]-x0[k-1]))
                            if p[i][j][k]<0:
                                p[i][j][k]=0
                        elif k==0:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/((x0[k]-(xmesh[0][i][j][k]-dx[0])**multifact)*(x2[i]-x2[i-1])*(x1[j]-x1[j-1]))
                            if p[i][j][k]<0:
                                p[i][j][k]=0
                        else:
                            p[i][j][k]=xmesh[3][i][j][k]*dx[0]*dx[1]*dx[2]/((x0[k]-x0[k-1])*(x1[j]-x1[j-1])*(x2[i]-x2[i-1]))
            pintplt = interpolate.Rbf(xx, yy, zz, p, function='linear', smooth=0)
            p00=np.zeros((divt[0],divt[1],divt[2]))
            t0=[[],[],[]]
            t0[0],t0[1],t0[2]=np.meshgrid(t[0],t[1],t[2])
            t1=[t0[0].flatten(),t0[1].flatten(),t0[2].flatten()]
            for i in range(divt[2]):
                for j in range(divt[1]):
                    for k in range(divt[0]):
                        p00[i][j][k]=pintplt(t0[0][i][j][k],t0[1][i][j][k],t0[2][i][j][k])
                        menseki+=p00[i][j][k]*dt[0]*dt[1]*dt[2]
            if menseki>0:
                p0=p00/menseki
            else:
                p0=p00
            p1=p0.flatten()
        if shw==0:
            print('multiplier:menseki',menseki)
        if dim==1:
            W=[t[0]]
            W0=[t[0],p0]
            W1=[t[0],p1]
        elif dim==2:
            W=[list(t[0]),list(t[1])]
            W0=[t0[0],t0[1],p0]
            W1=[t1[0],t1[1],p1]
        elif dim==3:
            W=[t[0],t[1],t[2]]
            W0=[t0[0],t0[1],t0[2],p0]
            W1=[t1[0],t1[1],t1[2],p1]
        WW=bunpu()
        WW.para=W
        WW.mesh=W0
        WW.flatten=W1
        WW.div=divt
        WW.dx=dt
        WW.dim=dim
        WW.xmin=tmin
        WW.xmax=tmax
        WW.pmax=p1.max()
        return WW
                


    def bunpu_division(self,other,divt=[],divs=[],corel=[0],shw=0):
        filename='x・y_log'
        logname=filename+'.csv'
        log = open(logname,'w' , encoding='shift_jis' )
        info=['tパラメータ','p値']
        writer = csv.writer(log, lineterminator='\n')
        writer.writerow(info)
        x,x2,x20,xmin,xmax,dim,dx,divx,m,y0min,y0max,divt0,flag0=self.bunpu_filter(other,1,divt)

        if flag0!=0:
            if flag0!=3:
                Y=m
                if divt==[]:
                    divt=divx
                if dim==1:
                    y0=np.linspace(y0min[0],y0max[0],len(x2[0]))
                    t10=[x2[0][i]/y0[i] for i in range(len(x2[0]))]
                    t00=t10
                    tmin=(min(xmin[0]/y0min[0],xmin[0]/y0max[0],xmax[0]/y0min[0],xmax[0]/y0max[0]))
                    tmax=(max(xmax[0]/y0max[0],xmin[0]/y0max[0],xmax[0]/y0min[0],xmax[0]/y0max[0]))
                    p1=x2[1]*(xmax[0]-xmin[0])/(tmax-tmin)
                    p0=x2[1]
                    divtt=divt[0]
                    zfun_smooth_rbf = interpolate.Rbf(t00, x20[1], function='linear', smooth=0)#
                elif dim==2:

                    y0=[np.linspace(y0min[0],y0max[0],len(x2[0])),np.linspace(y0min[0],y0max[0],len(x2[1]))]
                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]/y0[0][i] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]/y0[1][i] for i in range(len(x2[1]))])
                    t00.append([[x20[0][i][j]/y0[0][j] for j in range(len(x20[0,:]))] for i in range(len(x20[:,0]))])
                    t00.append([[x20[1][i][j]/y0[1][i] for j in range(len(x20[0,:]))] for i in range(len(x20[:,0]))])
                    
                    tmin.append(min(xmin[0]/y0min[0],xmin[0]/y0max[0],xmax[0]/y0min[0],xmax[0]/y0max[0]))
                    tmin.append(min(xmin[1]/y0min[0],xmin[1]/y0max[0],xmax[1]/y0min[0],xmax[1]/y0max[0]))
                    tmax.append(max(xmin[0]/y0min[0],xmin[0]/y0max[0],xmax[0]/y0min[0],xmax[0]/y0max[0]))
                    tmax.append(max(xmin[1]/y0min[0],xmin[1]/y0max[0],xmax[1]/y0min[0],xmax[1]/y0max[0]))
                    p1=x2[2]
                    p0=x20[2]
                    divtt=divt[0]*divt[1]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], x20[2], function='linear', smooth=0)#

                elif dim==3:
                    y0=[np.linspace(y0min[0],y0max[0],len(x2[0])),np.linspace(y0min[0],y0max[0],len(x2[1])),np.linspace(y0min[0],y0max[0],len(x2[2]))]
                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]/y0[0][i] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]/y0[1][i] for i in range(len(x2[1]))])
                    t10.append([x2[2][i]/y0[2][i] for i in range(len(x2[2]))])
                    t00.append([[[x20[0][i][j][k]/y0[0][j] for k in range(len(x20[0][0,0,:]))] for j in range(len(x20[0][0,:,0]))] for i in range(len(x20[0][:,0,0]))])
                    t00.append([[[x20[1][i][j][k]/y0[1][i] for k in range(len(x20[1][0,0,:]))] for j in range(len(x20[1][0,:,0]))] for i in range(len(x20[1][:,0,0]))])
                    t00.append([[[x20[2][i][j][k]/y0[2][k] for k in range(len(x20[2][0,0,:]))] for j in range(len(x20[2][0,:,0]))] for i in range(len(x20[2][:,0,0]))])
                    tmin.append(min(xmin[0]/y0min[0],xmin[0]/y0max[0],xmax[0]/y0min[0],xmax[0]/y0max[0]))
                    tmin.append(min(xmin[1]/y0min[0],xmin[1]/y0max[0],xmax[1]/y0min[0],xmax[1]/y0max[0]))
                    tmin.append(min(xmin[2]/y0min[0],xmin[2]/y0max[0],xmax[2]/y0min[0],xmax[2]/y0max[0]))
                    tmax.append(max(xmin[0]/y0min[0],xmin[0]/y0max[0],xmax[0]/y0min[0],xmax[0]/y0max[0]))
                    tmax.append(max(xmin[1]/y0min[0],xmin[1]/y0max[0],xmax[1]/y0min[0],xmax[1]/y0max[0]))
                    tmax.append(max(xmin[2]/y0min[0],xmin[2]/y0max[0],xmax[2]/y0min[0],xmax[2]/y0max[0]))
                    p1=x2[3]
                    p0=x20[3]
                    divtt=divt[0]*divt[1]*divt[2]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], t00[2], x20[3], function='linear', smooth=0)#
                tt=[]
                dt=[]
                for i in range(dim):
                    tt.append(np.linspace(tmin[i],tmax[i],divx[i]))
                    dt.append((tmax[i]-tmin[i])/(divt[i]-1))
                menseki=0
                if dim == 1:
                    t=[t00]
                    t1=t
                    t0=t
                    p0=np.zeros(divt[0])
                    p1=np.zeros(divt[0])
                elif dim == 2:
                    t=[tt[0],tt[1]]
                    ta,tb=np.meshgrid(tt[0],tt[1])
                    t1=[ta.flatten(),tb.flatten()]
                    t0=[ta,tb]
                    p0=np.zeros((divt[1],divt[0]))
                    p1=np.zeros(divtt)
                elif dim == 3:
                    t=[tt[0],tt[1],tt[2]]
                    ta,tb,tc=np.meshgrid(tt[0],tt[1],tt[2])
                    t1=[ta.flatten(),tb.flatten(),tc.flatten()]
                    t0=[ta,tb,tc]
                    p0=np.zeros((divt[2],divt[1],divt[0]))
                    p1=np.zeros(divtt)
                    
                for ti in range(divtt):
                    if dim==1:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti])*(xmax[0]-xmin[0])/(tmax[0]-tmin[0])
                        p0[ti]=p1[ti]
                        menseki+=p1[ti]*dt[0]
                    elif dim==2:
                        
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1]))
                        tq,tmod=divmod(ti,divt[0])

                        p0[tq][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1])
                    elif dim==3:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti],t1[2][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])\
                            /((tmax[0]-tmin[0])*(tmax[1]-tmin[1])*(tmax[2]-tmin[2]))
                        tq,tmod=divmod(ti,divt[0])
                        tq1,tq2=divmod(tq,divt[1])
                        p0[tq1][tq2][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
            #
            else:
                Y=other
                if divt==[]:
                    divt=divx
                if dim==1:
                
                    t10=[x2[0][i]/Y[0] for i in range(len(x2[0]))]
                    t00=[x2[0][i]/Y[0] for i in range(len(x2[0]))]
                    p1=x2[1]
                    p0=x2[1]
                    tmin=[xmin[0]/Y[0]]
                    tmax=[xmax[0]/Y[0]]
                    divtt=divt[0]
                    zfun_smooth_rbf = interpolate.Rbf(t00, x20[1], function='linear', smooth=0)#
                elif dim==2:

                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]/Y[0] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]/Y[0] for i in range(len(x2[1]))])
                    t00.append([[x20[0][i][j]/Y[0] for j in range(len(x20[0][0]))] for i in range(len(x20[0]))])
                    t00.append([[x20[1][i][j]/Y[0] for j in range(len(x20[1][0]))] for i in range(len(x20[1]))])
                    
                    tmin.append(xmin[0]/Y[0])
                    tmin.append(xmin[1]/Y[0])
                    tmax.append(xmax[0]/Y[0])
                    tmax.append(xmax[1]/Y[0])
                    p1=x2[2]
                    p0=x20[2]
                    divtt=divt[0]*divt[1]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], x20[2], function='linear', smooth=0)#
                elif dim==3:
                    t10=[]
                    t00=[]
                    tmin=[]
                    tmax=[]
                    t10.append([x2[0][i]/Y[0] for i in range(len(x2[0]))])
                    t10.append([x2[1][i]/Y[0] for i in range(len(x2[1]))])
                    t10.append([x2[2][i]/Y[0] for i in range(len(x2[2]))])
                    t00.append([[[x20[0][i][j][k]/Y[0] for k in range(len(x20[0][0][0]))] for j in range(len(x20[0][0]))] for i in range(len(x20[0]))])
                    t00.append([[[x20[1][i][j][k]/Y[0] for k in range(len(x20[1][0][0]))] for j in range(len(x20[1][0]))] for i in range(len(x20[1]))])
                    t00.append([[[x20[2][i][j][k]/Y[0] for k in range(len(x20[2][0][0]))] for j in range(len(x20[2][0]))] for i in range(len(x20[2]))])
                    tmin.append(xmin[0]/Y[0])
                    tmin.append(xmin[1]/Y[0])
                    tmin.append(xmin[2]/Y[0])
                    tmax.append(xmax[0]/Y[0])
                    tmax.append(xmax[1]/Y[0])
                    tmax.append(xmax[2]/Y[0])
                    p1=x2[3]
                    p0=x20[3]
                    divtt=divt[0]*divt[1]*divt[2]
                    zfun_smooth_rbf = interpolate.Rbf(t00[0],t00[1], t00[2], x20[3], function='linear', smooth=0)#
                tt=[]
                dt=[]
                for i in range(dim):
                    tt.append(np.linspace(tmin[i],tmax[i],divx[i]))
                    dt.append((tmax[i]-tmin[i])/(divt[i]-1))
                menseki=0
                if dim == 1:
                    t=[t00]
                    t1=t
                    t0=t
                    p0=np.zeros(divt[0])
                    p1=np.zeros(divt[0])
                elif dim == 2:
                    t=[tt[0],tt[1]]
                    ta,tb=np.meshgrid(tt[0],tt[1])
                    t1=[ta.flatten(),tb.flatten()]
                    t0=[ta,tb]
                    p0=np.zeros((divt[1],divt[0]))
                    p1=np.zeros(divtt)
                elif dim == 3:
                    t=[tt[0],tt[1],tt[2]]
                    ta,tb,tc=np.meshgrid(tt[0],tt[1],tt[2])
                    t1=[ta.flatten(),tb.flatten(),tc.flatten]
                    t0=[ta,tb,tc]
                    p0=np.zeros((divt[2],divt[1],divt[0]))
                    p1=np.zeros(divtt)
                    
                for ti in range(divtt):
                    if dim==1:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti])*(xmax[0]-xmin[0])/(tmax[0]-tmin[0])
                        p0[ti]=p1[ti]
                        menseki+=p1[ti]*dt[0]
                    elif dim==2:
                        
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])/((tmax[0]-tmin[0])*(tmax[1]-tmin[1]))
                        tq,tmod=divmod(ti,divt[0])

                        p0[tq][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1])
                    elif dim==3:
                        p1[ti]=zfun_smooth_rbf(t1[0][ti],t1[1][ti],t1[2][ti])*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])\
                            /((tmax[0]-tmin[0])*(tmax[1]-tmin[1])*(tmax[2]-tmin[2]))
                        tq,tmod=divmod(ti,divt[0])
                        tq1,tq2=divmod(tq,divt[1])
                        p0[tq1][tq2][tmod]=p1[ti]
                        menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
                #
        else:
            y2=other.flatten

            if other.xmin==[] or other.xmax==[] or other.dim==0 or other.div==[]:
                yf,py0,y0,dimy,divy,dy,ymin,ymax=other.kaiseki()
                other.xmin=ymin
                other.xmax=ymax
                other.dim=dimy
                other.div=divy
                other.dx=dy
            else:
                ymin=other.xmin
                ymax=other.xmax
                dimy=other.dim
                divy=(other.div)[0]
                dy=other.dx
            divxt=1
            for i in range(dim):
                divxt*=divx[i]

            corner1=np.array(xmin)/ymin[0]
            corner2=np.array(xmin)/ymax[0]
            corner3=np.array(xmax)/ymin[0]
            corner4=np.array(xmax)/ymax[0]
            tmin=[]
            tmax=[]
            dt=[]
            if divt==[]:
                divt=divx
            if divs==[]:
                divs=copy.copy(divt)

            for i in range(dim):
                tmin.append(min(corner1[i],corner2[i],corner3[i],corner4[i]))
                tmax.append(max(corner1[i],corner2[i],corner3[i],corner4[i]))
                dt.append((tmax[i]-tmin[i])/(divt[i]-1))#

            if dim==1:
                divtt=divt[0]
                t=[np.linspace(tmin[0],tmax[0],divt[0])]
                t0=t
                t1=t
                p1=np.zeros(divt[0])
                p0=np.zeros(divt[0])
                lastp=np.zeros(divt[0])
                base_x = interpolate.Rbf(x2[0], x2[1], function='linear', smooth=0)  #
                base_y = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0)  #
                
                signt=[0 if i < divt[0]-1 and np.sign(t1[0][i])!=np.sign(t1[0][i+1]) else np.sign(t1[0][i]) for i in range(divt[0])]
            elif dim ==2:
                divtt=divt[0]*divt[1]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1])]
                tx,ty=np.meshgrid(t[0],t[1])#
                t1=[tx.flatten(),ty.flatten()]
                xx,xy=np.meshgrid(np.linspace(xmin[0],xmax[0],divt[0]),np.linspace(xmin[1],xmax[1],divt[1]))#
                t0=[tx,ty]
                p1=np.zeros(divtt)
                p0=np.zeros((divt[1],divt[0]))
                lastp=np.zeros(divtt)
                base_x = interpolate.Rbf(x20[0], x20[1], x20[2], function='linear', smooth=0)  #
                base_y = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0)  #
            elif dim==3:
                divtt=divt[0]*divt[1]*divt[2]
                t=[np.linspace(tmin[0],tmax[0],divt[0]),np.linspace(tmin[1],tmax[1],divt[1]),np.linspace(tmin[2],tmax[2],divt[2])]
                tx,ty,tz=np.meshgrid(t[0],t[1],t[2])#
                t1=[tx.flatten(),ty.flatten(),tz.flatten()]
                xx,xy,xz=np.meshgrid(np.linspace(xmin[0],xmax[0],divt[0]),np.linspace(xmin[1],xmax[1],divt[1]),np.linspace(xmin[2],xmax[2],divt[2]))#
                t0=[tx,ty,tz]
                p1=[0.]*divtt
                p0=np.zeros((divt[2],divt[1],divt[0]))
                lastp=np.zeros(divtt)
                base_x = interpolate.Rbf(x20[0], x20[1], x20[2], x20[3], function='linear', smooth=0)  #
                base_y = interpolate.Rbf(y2[0], y2[1], function='linear', smooth=0)  #

            menseki=0
            test1x=[[],[]]
            test1y=[[],[]]
            test2x=[[],[]]
            test2y=[[],[]]
            for ti in range(divtt):#
                max_x=0
                max_y=0
                max_p=0
                if dim==1:
                    dts=copy.copy(dt[0])
                    max_s=xmax[0]
                    min_s=xmin[0]
                    min_t=tmin[0]
                    max_t=tmax[0]
                    divst=divtt
                    tti=ti
                if dim==2:
                    tq,tmod=divmod(ti,divt[0])
                    divst, min_t, max_t, min_s, max_s, min_sx, max_sx, min_sy, max_sy, plmin, \
                        plmax=slice2d_bunpu(ti,t1,tmin,tmax,dt,divs,xmin,xmax,ymin,ymax)
                    if plmin!=9 and plmax!=9 and max_t>min_t:
                        if divst<=1:
                            divst=1
                        dts=(max_t-min_t)/divst
                        tspan=(t1[0][ti]**2+t1[1][ti]**2)**0.5#
                        tti=(tspan-min_t)//dts
                        divs[0]=divst
                        divss=int(divs[0])
                        s=np.linspace(min_s,max_s,divss)
                        sx=np.linspace(min_sx,max_sx,divss)
                        sy=np.linspace(min_sy,max_sy,divss)
                        if divss<=1:
                            ds=(max_s-min_s)
                        else:
                            ds=(max_s-min_s)/(divss-1)
                        ps=[]
                        for i in(range(divss)):
                            ps.append(max(0,base_x(sx[i],sy[i])))
                        divs[0]=divss
                        if divs[0]<=1:
                            da=(max_s+ymax[0]-(min_s+ymin[0]))
                        else:
                            da=(max_s+ymax[0]-(min_s+ymin[0]))/(divs[0]-1)#
                        if len(s)>1 and max_s>min_s:
                            base_s = interpolate.Rbf(s, ps, function='linear', smooth=0)
                        else:
                            base_s=interpolate.Rbf([s[0]-ds,s[0],s[0]+ds], [ps[0],ps[0],ps[0]], function='linear', smooth=0)
                        signt=[0 if i < divss-1 and np.sign(s[i])!=np.sign(s[i+1]) else np.sign(s[i]) for i in range(divss)]
                    else:
                        divs[0]=0
                if dim==3:
                    tq,tmod=divmod(ti,divt[0])
                    tq1,tq2=divmod(tq,divt[1])
                    divst, min_t, max_t, min_s, max_s, min_sx, max_sx, min_sy, max_sy, min_sz, max_sz, \
                        plmin, plmax=slice3d_bunpu(ti,t1,tmin,tmax,dt,divs,xmin,xmax,ymin,ymax)
                    if plmin!=9 and plmax!=9  and abs(max_s-min_s)>min((xmax[0]-xmin[0])/divx[0],(xmax[1]-xmin[1])/divx[1],(xmax[2]-xmin[2])/divx[2]):
                        if divst<=1:
                            divst=1
                        dts=(max_t-min_t)/divst
                        tspan=(t1[0][ti]**2+t1[1][ti]**2+t1[2][ti]**2)**0.5#

                        tti=-1*(-1*(tspan-min_t)//dts)
                        divs[0]=divst
                        divss=int(divs[0])
                        s=np.linspace(min_s,max_s,divss)#
                        sx=np.linspace(min_sx,max_sx,divss)
                        sy=np.linspace(min_sy,max_sy,divss)
                        sz=np.linspace(min_sz,max_sz,divss)
                        if divss<=1:
                            ds=(max_s-min_s)
                        else:
                            ds=(max_s-min_s)/(divss-1)
                        ps=[]
                        for i in(range(divss)):
                            ps.append(max(0,base_x(sx[i],sy[i],sz[i])))
                        divs[0]=divss
                        if divs[0]<=1:
                            da=(max_s+ymax[0]-(min_s+ymin[0]))
                        else:
                            da=(max_s+ymax[0]-(min_s+ymin[0]))/(divs[0]-1)#
                        if len(s)>1 and max_s>min_s:
                            base_s = interpolate.Rbf(s, ps, function='linear', smooth=0)
                        else:
                            base_s=interpolate.Rbf([s[0]-ds,s[0],s[0]+ds], [ps[0],ps[0],ps[0]], function='linear', smooth=0)
                        signt=[0 if i < divss-1 and np.sign(s[i])!=np.sign(s[i+1]) else np.sign(s[i]) for i in range(divss)]
                    else:
                        divs[0]=0                        
                #
                if divst<2:
                    divst2=2
                else:
                    divst2=divst
                if divs[0]>0:
                    bi=((divst2-1-tti)*min_s/ymax[0]+tti*max_s/ymin[0])/(divst2-1)
                    bi1=((divst2-1-(tti+1))*min_s/ymax[0]+(tti+1)*max_s/ymin[0])/(divst2-1)
                men=0
                cnt=0
                for js in range(divs[0]):#
                    flag=np.zeros(dim)
                    xh=np.zeros(dim)
                    #
                    if signt[js]>0:
                        if min_s>0 or max_s<0:
                            da=((ymax[0]-ymin[0]+max_s-min_s)/(divst2-1))#
                        else:
                            da=((ymax[0]-ymin[0]+max_s-min_s)/(divst2-1))
                        yij=((ymin[0]+min_s)+js*da)/(bi+1)
                        yij1=((ymin[0]+min_s)+(js+1)*da)/(bi+1)
                        yi1j=((ymin[0]+min_s)+js*da)/(bi1+1)
                        xij=((ymin[0]+min_s)+js*da)/(1/bi+1)
                        yi1j1=((ymin[0]+min_s)+(js+1)*da)/(bi1+1)
                        xij1=((ymin[0]+min_s)+(js+1)*da)/(1/bi+1)
                        xi1j=((ymin[0]+min_s)+js*da)/(1/bi1+1)
                        xi1j1=((ymin[0]+min_s)+(js+1)*da)/(1/bi1+1)
                        spxij=(xi1j1-xij)
                        spyij=(yij1-yi1j)
                    else:
                        if min_s>0 or max_s<0:
                            da=((ymax[0]-ymin[0]+max_s-min_s)/(divst2-1))
                        else:
                            da=((ymax[0]-ymin[0]+max_s-min_s)/(divst2-1))
                        yij=((min_s-ymax[0])+js*da)/(bi-1)
                        yij1=((min_s-ymax[0])+(js+1)*da)/(bi-1)
                        yi1j=((min_s-ymax[0])+js*da)/(bi1-1)
                        xij=((min_s-ymax[0])+js*da)/(1-1/bi)
                        yi1j1=((min_s-ymax[0])+(js+1)*da)/(bi1-1)
                        xij1=((min_s-ymax[0])+(js+1)*da)/(1-1/bi)
                        xi1j=((min_s-ymax[0])+js*da)/(1-1/bi1)
                        xi1j1=((min_s-ymax[0])+(js+1)*da)/(1-1/bi1)
                        spxij=xi1j1-xij
                        spyij=yi1j-yij1
                    darea=(spxij*spyij)/2#
                    
                    flag2=cor2ran(xij,min_s,max_s,yij,ymin[0],ymax[0],corel[0])
                    #
                    if flag2==1.0:
                        if dim==1:

                            px = base_x(xij)
                            py = base_y(yij)
                            pp=px*py*darea/dts
                            p1[ti]+=max(0,pp)
                            p0[ti]=p1[ti]
                            #
                        elif dim==2:
                            px = base_s(xij)
                            py = base_y(yij)
                            pp=(xij/tspan)*px*py*darea/dts
                            p1[ti]+=max(0,pp)
                            p0[tq][tmod]= p1[ti]

                        elif dim==3:
                            px = max(0,base_s(xij))
                            py = base_y(yij)
                            p1[ti]+=((xij/tspan)**2)*px*py*darea/dts
                            p0[tq1][tq1][tmod]= p1[ti]
                        else:
                            print('error')
                if dim==1:
                    menseki+=p1[ti]*dt[0]
                if dim==2:
                    menseki+=p1[ti]*(dt[0]*dt[1])
                elif dim==3:
                    menseki+=p1[ti]*(dt[0]*dt[1]*dt[2])
        #
        if shw==0:
            flaginf=['bunpu/bunpu','bunpu/lean','lean/bunpu','bunpu/vector']
            print(flaginf[flag0],menseki)
        if menseki>0:
            if dim==1:
                p00=p0/menseki
                p=p1/menseki
            else:
                p=p1/menseki
                p00=p0/menseki
        else:
            if dim==1:
                p00=p0
                p=p1
            else:
                p=p1
                p00=p0
        #
        if dim==1:
            W=[t[0]]
            W0=[t0[0],p00]
            W1=[t1[0],p]
        elif dim==2:
            W=[list(t[0]),list(t[1])]
            W0=[t0[0],t0[1],p00]
            W1=[t1[0],t1[1],p]
        elif dim==3:
            W=[t[0],t[1],t[2]]
            W0=[t0[0],t0[1],t0[2],p00]
            W1=[t1[0],t1[1],t1[2],p]
        WW=bunpu()
        WW.para=W
        WW.mesh=W0
        WW.flatten=W1
        WW.div=divt
        WW.dx=dt
        WW.dim=dim
        WW.xmin=tmin
        WW.xmax=tmax
        WW.pmax=max(p)
        return WW
                
        
    def bunpu_shrink(self,limitp=0.001,divlow=[]):

        X=self.para
        X1=self.flatten
        X0=self.mesh
        #
        xf,px0,x0,dim,divx,dx,xmin,xmax=self.kaiseki()#
        if self.xmin==[] or self.xmax==[] or self.dim==0 or self.div==[] or self.dx==[]:
            self.xmin=xmin
            self.xmax=xmax
            self.dim=dim
            self.div=divx
            self.dx=dx
        else:
            xmin=self.xmin
            xmax=self.xmax
            dim=self.dim
            #divx=self.div
            dx=self.dx
        if divlow==[]:
            divlow=[10*dim]
        
        #
        if dim==1:
            menseki_low=0
            long_low=0
            cnt_low=0
            menseki_hi=0
            long_hi=0
            cnt_hi=divx[0]
            for i in range(divx[0]):
                menseki_low += X0[1][i]*dx[0]
                
                if menseki_low < limitp[0] and cnt_low+1<cnt_hi:
                    cnt_low=long_low//dx[0]
                long_low += dx[0]
                menseki_hi += X0[1][divx[0]-i-1]*dx[0]

                if menseki_hi < limitp[0] and cnt_low+1<cnt_hi:
                    cnt_hi=divx[0]-long_hi//dx[0]
                long_hi += dx[0]
            t0=[X0[0][i] for i in range(divx[0]) if i >= cnt_low-1 and i <= cnt_hi+1]
            p0=[X0[1][i] for i in range(divx[0]) if i >= cnt_low-1 and i <= cnt_hi+1]
            #t0=[X0[0][i] for i in range(divx[0]) if i <= cnt_hi]
            #p0=[X0[1][i] for i in range(divx[0]) if i <= cnt_hi]
            if cnt_hi-cnt_low < divlow[0] and len(t0)>=2:
                t=np.linspace(min(t0),max(t0),divlow[0])
                base = interpolate.Rbf(t0, p0, function='linear', smooth=0)
                p=[base(t[i]) for i in range(divlow[0])]
            else:
                t=t0
                p=p0
            W=[t,p]
            W0=[t,p]
            W1=[t,p]
            tmax=[max(t)]
            tmin=[min(t)]
            divt=[len(t)]
            dt=[(tmax[0]-tmin[0])/(divt[0]-1)]
        elif dim>=2:
            print('under constract')
            
        WW=bunpu()
        WW.para=W
        WW.mesh=W0
        WW.flatten=W1
        WW.div=divt
        WW.dx=dt
        WW.dim=dim
        WW.xmin=tmin
        WW.xmax=tmax

        return WW        
                   
    
    def gene2(self,xmin,xmax,xmean,xdev,xskew,div=100,filename=''):#
        dim=len(xmean)
        if dim==1:
            ZZ=[]
            dx=(xmax-xmin)/div
            x=np.linspace(xmin,xmax,div)
            y=norm.pdf(x.xmean,xdev)

    def divide(self,dx):
        x=self.para
        x0=self.mesh
        dim=self.dim
        div1=self.div
        xmax=[]
        xmin=[]
        for i in range(dim):
            xmin.append(min(x[i]))
            xmax.append(max(x[i]))
        if div1 == 0:
            for i in range(dim):
                div1=len(x[i])

        if dim==1:
            zfun_smooth_rbf = interpolate.Rbf(x0[0], x0[1], function='linear', smooth=0)
            x2=np.arange(xmin[0],xmax[0],dx[0])
            y2=[]
            if np.max(x2) < xmax[0]:
                x2=np.append(x2,xmax[0])
            divx=len(x2)
            for i in range(divx):
                y2.append(zfun_smooth_rbf(x2[i]))#
            self.para2=[x2,y2]
            self.div2=[divx]
            self.dx2=[dx]
            return [x2,y2],[divx],[dx]
        elif dim>=2:
            
            x2=[]
            y2=[]
            divx=[]
            divm=1

            for i in range(dim):
                x2.append(list(np.arange(xmin[i],xmax[i],dx[i])))
                if np.max(x2[i]) < xmax[i]:
                    x2[i]=list(np.append(x2[i],xmax[i]))
                divx.append(len(x2[i]))
                divm*=divx[i]
                
            ytmp=[]
            if dim==2:

                p0=np.zeros((divx[1],divx[0]))
                xx2,xy2=np.meshgrid(x2[0],x2[1])
                xxf=xx2.flatten()
                xyf=xy2.flatten()
                zfun_smooth_rbf = interpolate.Rbf(x0[0], x0[1], x0[2], function='linear', smooth=0)
                for j in range(divm):
                    tq,tmod=divmod(j,divx[0])

                    ytmp=zfun_smooth_rbf(xxf[j],xyf[j])
                    if ytmp > 0:
                        p0[tq][tmod]=ytmp
                    else:
                        p0[tq][tmod]=0
                self.para2=[xx2,xy2,p0]
                self.div2=divx
                self.dx2=dx
                return [xx2,xy2,p0],divx,dx
            elif dim==3:
                p0=np.zeros((divx[2],divx[1],divx[0]))
                xx2,xy2,xz2=np.meshgrid(x2[0],x2[1],x2[2])
                xxf=xx2.flatten()
                xyf=xy2.flatten()
                xzf=xz2.flatten()
                zfun_smooth_rbf = interpolate.Rbf(x0[0], x0[1], x0[2], x0[3],function='linear', smooth=0)  #
                for j in range(divm):
                    tq,tmod=divmod(j,divx[0])
                    tq1,tq2=divmod(tq,divx[1])
                    ytmp=zfun_smooth_rbf(xxf[j],xyf[j],xzf[j])
                    if ytmp > 0:
                        p0[tq1][tq2][tmod]=ytmp
                    else:
                        p0[tq1][tq2][tmod]=0

                self.para2=[xx2,xy2,xz2,p0]
                self.div2=divx
                self.dx2=dx
                return [xx2,xy2,xz2,p0],divx,dx
            

    
    def bunpu_graph(self,filename='' ,gr=0,levels=5,linewidth=4,fsize=(14,7),yscale=0,xscale=0):
        fp = FontProperties(fname=r'C:\Windows\Fonts\msgothic.ttc',size=24)
        x=self.flatten
        if self.dim!=[]:
            dim=self.dim
        else:
            dim=len(x)-1
        xp=self.para
        xm=self.mesh
        dx=self.dx
        divn=self.div
        if dim==1:
            fig=plt.figure(figsize=(14,7))
            ax = fig.add_subplot(111)
            if filename!='':
                outgraphname=filename+'.png'
            if yscale != 0:
                ax.ylim(yscale)
            if xscale != 0:
                ax.xlim(xscale)
            ax.plot(x[0], x[1], label="origin",linewidth = 4)
            ax.set_title(u''+'分布', fontproperties=fp,fontsize=24)
            ax.set_xlabel(u''+'X軸', fontproperties=fp,fontsize=24)
        
        if dim==2:
            if gr==0:

                if filename!='':
                    outgraphname=filename+'.png'
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.plot_surface(xm[0], xm[1], xm[2],cmap='viridis',linewidth=0)
                ax.set_xlabel('X axis')
                ax.set_ylabel('Y axis')
                ax.set_zlabel('Z axis')
            elif gr==1:

                if filename!='':
                    outgraphname=filename+'_h.png'
                fig=plt.figure(figsize=(14,7))
                cont=plt.contour(xm[0],xm[1],xm[2],levels,colors=['black'])#
                cont.clabel(fmt='%1.6f', fontsize=14)
                
                plt.xlabel('X axis', fontsize=24)
                plt.ylabel('Y axis', fontsize=24)

        elif dim==3:
            if filename!='':
                outgraphname=filename+'.png'
            p_mins = min(xm[3].flatten())
            p_maxs = max(xm[3].flatten())
            if p_mins==p_maxs:
                p_maxs=p_mins+1
            fig = plt.figure()
            xmin0=xp[0].min()
            xmax0=xp[0].max()
            xmin1=xp[1].min()
            xmax1=xp[1].max()
            xmin2=xp[2].min()
            xmax2=xp[2].max()
            gr=3
            if gr==0 or gr==2:
                cmap = cm.binary
                cmap_data = cmap(np.arange(cmap.N))
                lnmap=len(cmap_data)
                for i in range(lnmap):
                    trn=i/lnmap
                    cmap_data[i, 3] = trn 
                customized_binary = colors.ListedColormap(cmap_data)
                scam = plt.cm.ScalarMappable(
                    norm=cm.colors.Normalize(p_mins, p_maxs),
                    cmap=customized_binary
                )
                ax  = fig.gca(projection='3d')
                    
                div=len(xm[0])
                for i in range(div):
                    q,mod=divmod(i,div)
                    q1,q2=divmod(q,div)
                    scam.set_array([])
                    ax.plot_surface(
                        xm[0][i],xm[1][i],xm[2][i],#
                        facecolors  = scam.to_rgba(xm[3][i]),
                        antialiased = False,
                        rstride=1, cstride=1,#
                        edgecolor='None'
                    )
                ax.view_init(elev=30, azim=10)
            elif gr==1 or gr==3:
                ax = fig.add_subplot(111, projection='3d')
                plt.figaspect(1)
                cmap = cm.binary
                cmap_data = cmap(np.arange(cmap.N))
                lnmap=len(cmap_data)
                for i in range(lnmap):
                    trn=i/lnmap
                    cmap_data[i, 3] = trn #
                customized_binary = colors.ListedColormap(cmap_data)
                xm0,xm1,xm2=np.mgrid[xmin0:xmax0:(xmax0-xmin0)/divn[0],xmin1:xmax1:(xmax1-xmin1)/divn[1],xmin2:xmax2:(xmax2-xmin2)/divn[2]]
                sc = ax.scatter(xm0, xm1, xm2, c=xm[3].ravel(), cmap=customized_binary)
                fig.colorbar(sc)
            if gr==2 or gr==3:
                ax.set_xlim((x[0].min(),x[0].max()))
                ax.set_ylim((x[1].min(),x[1].max()))
                ax.set_zlim((x[2].min(),x[2].max()))
                
                plotdat = np.sum(xm[3], axis=2)
                plotdat = plotdat / np.max(plotdat)
                plotx, ploty=np.mgrid[xmin0:xmax0:(xmax0-xmin0)/divn[0],xmin1:xmax1:(xmax1-xmin1)/divn[1]]
                
                ax.contour(plotx, ploty, plotdat, offset=xp[2].min(), zdir='z')
                plotdat = np.sum(xm[3], axis=1) 
                plotdat = plotdat / np.max(plotdat)
                plotx, plotz=np.mgrid[xmin0:xmax0:(xmax0-xmin0)/divn[0],xmin2:xmax2:(xmax2-xmin2)/divn[2]]
                ax.contour(plotx, plotdat, plotz, offset=xp[1].max(), zdir='y')
                plotdat = np.sum(xm[3], axis=0) 
                plotdat = plotdat / np.max(plotdat)
                ploty, plotz=np.mgrid[xmin1:xmax1:(xmax1-xmin1)/divn[1],xmin2:xmax2:(xmax2-xmin2)/divn[2]]
                ax.contour(plotdat, ploty, plotz, offset=xp[0].min(), zdir='x')
                
        if filename!='':
            plt.savefig(outgraphname)
            plt.close()
        else:
            plt.show()

    def skew(self,filename = ''):
        s,m,xf,p,x,dim,divt,dx,xmin,xmax=self.sdev()

        sk0=np.zeros(dim)
        sk=[]
        for i in range(dim):
            for j in range(divt):
                sk0[i] += ((xf[i][j]-m[i])**3)*(np.clip(p[j],0,None))*dx[i]/(sdev**3)

            s.append(s0[i]**0.5)
            
        if filename != '':
            g = open(filename, 'w')
            writer = csv.writer(g, lineterminator='\n') 
            writer.writerow(['平均',m])
            writer.writerow(['分散',s])

        return sk,s,m,xf,p,x,dim,divt,dx,xmin,xmax
    
    def sdev(self,filename = ''):
        m,xf,p,x,dim,divt,dx,xmin,xmax=self.mean()
        
        s0=np.zeros(dim)
        s=[]
        for i in range(dim):
            for j in range(divt):
                s0[i] += ((xf[i][j]-m[i])**2)*np.clip(p[j],0,None)*dx[i]

            s.append(s0[i]**0.5)
        if filename != '':
            g = open(filename, 'w')
            writer = csv.writer(g, lineterminator='\n')
            writer.writerow(['平均',m])
            writer.writerow(['分散',s])
        return s,m,xf,p,x,dim,divt,dx,xmin,xmax
    
    def mean(self,filename = ''):

        dim=self.dim
        xf=self.flatten
        divt=len(xf[0])
        dx=self.dx
        m=np.zeros(dim)
        if dim==1:
            dxx=dx[0]
        elif dim==2:
            dxx=dx[0]*dx[1]
        elif dim==3:
            dxx=dx[0]*dx[1]*dx[2]
        for i in range(dim):
            for j in range(divt):

                m[i]+=xf[i][j]*(np.clip(xf[dim][j],0,None))*dxx
        if filename != '':
            g = open(filename, 'w')
            writer = csv.writer(g, lineterminator='\n') 
            writer.writerow(['平均',m])
        return m

    def bunpu_file(self,outfilename):
        x0=self.para
        x=self.mesh
        dim=len(x)-1
        outfile=outfilename
        if outfilename != '':
            if dim==1:
                df1=pd.DataFrame(x0)
                df2=pd.DataFrame(x[1])
            elif dim==2:
                df1=pd.DataFrame(x0)
                df2=pd.DataFrame(x[2])
            elif dim==3:
                df1=pd.DataFrame(x0)
                xx=[]
                for i in range(len(x[3])):
                    for j in range(len(x[3][0])):
                        xx.append(x[3][i][j])
                    xx.append(['*']*len(x[3][i][j]))
                df2=pd.DataFrame(xx)

            df1.to_csv(outfilename+'t.csv')

            df2.to_csv(outfilename+'p.csv')

                
        else:
            return df
    
    def kaiseki(self):
        #X=np.array(self.mesh)
        XX=self.mesh
        if type(XX)==list:
            X=np.array(XX)
        else:
            X=XX
        X1=self.para
        dim=len(X)-1
        xmin=[]
        xmax=[]
        x=[]
        div=[]
        divt=[]
        dx=[]
        xf=[]
        if dim == 1:
            p=X[dim]
        else:
            p=X[dim].flatten()
        dmen=1
        for i in range(dim):
            if dim==1:
                x.append(X[i])
                div.append(len(X[i]))
                divt.append(div[i])
                xmin.append(min(X[i]))
                xmax.append(max(X[i]))
                #xf.append(X[i].tolist())
                xf.append(X[i])
            elif dim==2:
                if i==0:
                    x.append(X[0][0,:])
                if i==1:
                    x.append(X[1][:,0])
                    
                div.append(len(x[i]))
                divt.append(len(x[i]))
                xmin.append(min(x[i]))
                xmax.append(max(x[i]))
                xf.append(X[i].flatten())

            elif dim==3:
                if i==0:
                    x.append(X[0][0,:,0])
                if i==1:
                    x.append(X[1][:,0,0])
                if i==2:
                    x.append(X[2][0,0,:])
                div.append(len(x[i]))
                divt.append(div[i])
                xmin.append(min(x[i]))
                xmax.append(max(x[i]))
                xf.append(X[i].flatten())

            div.append(len(x[i]))
            dx.append((xmax[i]-xmin[i])/(div[i]-1))
            dmen*=dx[i]
        return xf[0],p,x,dim,divt,dx,xmin,xmax

    def bunpu_split(self):
        xf,px,x,dim,divt,dx,xmin,xmax=self.kaiseki()
        m=np.zeros((divt,divt))
        tmin=np.array(xmin)/2
        tmax=np.array(xmax)/2
        dt=np.array(dx)/2
        
        ca=self.case
        W=bunpu()
        W=self.bunpu_muls(0.5)
        W.case=ca
        V=bunpu()
        V=W.bunpu_muls(2)
        vf,pv,v,dimv,divv,dxv,vmin,vmax=V.kaiseki()
        ratiod=0.01
        diffp=ratiod*(np.array(px)-np.array(pv))

        df=pd.DataFrame(diffp)
        df.to_csv('diff.csv')

        for i in range(divt):
            if dim == 1:
                if (tmax[0]+tmin[0])/2 > xf[0][i]/2:
                    tkmin=tmin[0]
                    tkmax=xf[0][i]-tmin[0]#
                else:
                    tkmin=xf[0][i]-tmax[0]
                    tkmax=tmax[0]

                for j in range(divt):
                    if xf[0][j]/2>=tkmin and xf[0][j]/2<=tkmax:
                        m[j][i]=1
        df=pd.DataFrame(m)
        df.to_csv('m.csv')
        if np.linalg.matrix_rank(m) == divt:
            print('lank ok')
            inv_m=np.linalg.inv(m)
        else:
            print('lank down')#
            mc=m+np.eye(divt)*2
            inv_m=np.linalg.inv(mc)
            print(np.linalg.matrix_rank(m))

        df=pd.DataFrame(inv_m)
        df.to_csv('inv_m.csv')

        test=np.dot(m,inv_m)

        df=pd.DataFrame(test)
        df.to_csv('test.csv')
        
        dpr=np.reshape(diffp,(1,divt))
        dprt=dpr.T
        hosei=np.dot(inv_m,dprt)
        wp=W.flatten[dim]
        W.flatten[dim]=np.array(wp)+hosei.flatten()
        return W




    def file2cor(self,infilename,skipline,colum1,colum2):
    
        if type(infilename) is list:
            if len(infilename)==2:
                dirname=[]
                infile=infilename[0]
                sht=infilename[1]#
                if 'csv' in sht:
                    dirname=infilename[0]
                    infile=infilename[1]
            elif len(infilename)==3:
                dirname=infilename[0]
                infile=infilename[1]
                sht=infilename[2]
                
        else:
            infile=infilename
            dirname=[]
        filename=os.path.splitext(infile)
        outfilename=os.path.splitext(outdataname)
        if filename[1]=='.csv':
            if dirname==[]:
                df = pd.read_csv(infile,usecols=[colum1,colum2],skiprows=skipline,na_values=999)
                dt = df.values
            else:
                dt=[]
                allFiles = glob.glob(dirname+infile)
                for file_ in allFiles:
                    df = pd.read_csv(file_,usecols=[colum1,colum2],skiprows=skipline,na_values=999)
                    dt.extend(df.values)
            
        elif filename[1]=='.xlsx':
            if dirname==[]:
                df = pd.read_excel(infile,usecols=[colum1,colum2],skiprows=skipline,sheet_name=sht,na_values=999)
                dt = df.values
            else:
                dt=[]
                allFiles = glob.glob(dirname+infile)
                for file_ in allFiles:
                    df = pd.read_excel(file_,usecols=[colum1,colum2],skiprows=skipline,sheet_name=sht,na_values=999)
                    dt.extend(df.values)
        data=np.array(dt)
        x=data[:,0]
        y=data[:,1]
        meanx=sum(x)/len(x)
        meany=sum(y)/len(y)
        disco=0
        disx=0
        disy=0
        for i in range(len(data)):
            dx=x[i]-meanx
            dy=y[i]-meany
            disco+=dx*dy
            disx+=dx*dx
            disy+=dy*dy
        corel=disco/(disx**0.5*disy**0.5)
        return corel

    def bunpu_contour(self,vdiv):
        dim=self.dim
        div=self.div

        xs=(self.flatten[dim]).tolist()
        max_value = max(xs)
        max_index = xs.index(max_value)
        xe=(self.mesh[dim])
        if dim==1:
            xsmin=max(xs[0],xs[div[0]-1])
            div_index=[[],[]]
            div_value=[[],[]]

        elif dim==2:
            xsmin=max(np.concatenate([xe[0],xe[div[0]-1],xe[:,0],xe[:,div[1]-1]],0))#
            div_index=[[],[],[],[]]
            div_value=[[],[],[],[]]
            xst=[(self.flatten[0]).tolist(),(self.flatten[1]).tolist()]
        elif dim==3:
            xsmin=max(np.concatenate([xe[0].flatten(),xe[div[0]-1].flatten(),xe[:,0].flatten(),\
                                      xe[:,div[1]-1].flatten(),xe[:,:,0].flatten(),xe[:,:,div[2]-1].flatten()],0))
            div_index=[[],[],[],[],[],[],[],[]]
            div_value=[[],[],[],[],[],[],[],[]]
            xst=[(self.flatten[0]).tolist(),(self.flatten[1]).tolist(),(self.flatten[2]).tolist()]
        dp=round((max_value-xsmin)/(vdiv-1),2)
        indexes=[]#
        layer=[]
        for j in range(0,vdiv+1):
            if j==0:
                xdmin=-0.1
                xdmax=xsmin
                
            elif j==vdiv:
                xdmin=xsmin+dp*(j-1)
                xdmax=xsmin+dp*j+0.01
            else:
                xdmin=xsmin+dp*(j-1)
                xdmax=xsmin+dp*j
            layer.append(xdmax)
            indexes.append([i for i, e in enumerate(xs) if e > xdmin and e<=xdmax])#

            if dim==1:
                div_index[0].append([i for i in indexes[j] if i<=max_index])
                div_index[1].append([i for i in indexes[j] if i>max_index])
                div_value[0].append([xs[i] for i in div_index[0][j]])
                div_value[1].append([xs[i] for i in div_index[1][j]])
            if dim==2:
                div_index[0].append([i for i in indexes[j] if xst[0][i]<=xst[0][max_index] and xst[1][i]<=xst[1][max_index]])
                div_index[1].append([i for i in indexes[j] if xst[0][i]>xst[0][max_index] and xst[1][i]<=xst[1][max_index]])
                div_index[2].append([i for i in indexes[j] if xst[0][i]<=xst[0][max_index] and xst[1][i]>xst[1][max_index]])
                div_index[3].append([i for i in indexes[j] if xst[0][i]>xst[0][max_index] and xst[1][i]>xst[1][max_index]])
                div_value[0].append([xs[i] for i in div_index[0][j]])
                div_value[1].append([xs[i] for i in div_index[1][j]])
                div_value[2].append([xs[i] for i in div_index[2][j]])
                div_value[3].append([xs[i] for i in div_index[3][j]])

            if dim==3:
                div_index[0].append([i for i in indexes[j] if xst[0][i]<=xst[0][max_index] and xst[1][i]<=xst[1][max_index] and xst[2][i]<=xst[2][max_index]])
                div_index[1].append([i for i in indexes[j] if xst[0][i]>xst[0][max_index] and xst[1][i]<=xst[1][max_index] and xst[2][i]<=xst[2][max_index]])
                div_index[2].append([i for i in indexes[j] if xst[0][i]<=xst[0][max_index] and xst[1][i]>xst[1][max_index] and xst[2][i]<=xst[2][max_index]])
                div_index[3].append([i for i in indexes[j] if xst[0][i]>=xst[0][max_index] and xst[1][i]>=xst[1][max_index] and xst[2][i]<xst[2][max_index]])
                div_index[4].append([i for i in indexes[j] if xst[0][i]<xst[0][max_index] and xst[1][i]<xst[1][max_index] and xst[2][i]>=xst[2][max_index]])
                div_index[5].append([i for i in indexes[j] if xst[0][i]<=xst[0][max_index] and xst[1][i]>xst[1][max_index] and xst[2][i]>xst[2][max_index]])
                div_index[6].append([i for i in indexes[j] if xst[0][i]>xst[0][max_index] and xst[1][i]<=xst[1][max_index] and xst[2][i]>xst[2][max_index]])
                div_index[7].append([i for i in indexes[j] if xst[0][i]>xst[0][max_index] and xst[1][i]>xst[1][max_index] and xst[2][i]>xst[2][max_index]])
                div_value[0].append([xs[i] for i in div_index[0][j]])
                div_value[1].append([xs[i] for i in div_index[1][j]])
                div_value[2].append([xs[i] for i in div_index[2][j]])
                div_value[3].append([xs[i] for i in div_index[3][j]])
                div_value[4].append([xs[i] for i in div_index[4][j]])
                div_value[5].append([xs[i] for i in div_index[5][j]])
                div_value[6].append([xs[i] for i in div_index[6][j]])
                div_value[7].append([xs[i] for i in div_index[7][j]])
        return indexes, max_value, max_index,div_value,div_index,layer
                                    

    def bunpu_contour2(self,vdiv):
        dim=self.dim
        div=self.div
        
        xs=(self.flatten[dim]).tolist()
        max_value = max(xs)
        max_index = xs.index(max_value)
        xe=(self.mesh[dim])#
        xsmin=[]
        xsmax=[]
        xst=[]
        paras=[]
        parae=[]
        div_paras=[]
        div_paras2=[]
        div_try=[]
        if dim == 1:
            hk=0.005
        elif dim==2:
            hk=0.018
        elif dim==3:
            hk=0.1
            
        hinum=(np.geomspace(1,1+hk,vdiv)-1)/hk
        lonum=1-np.array(list(reversed(hinum)))
        for i in range(dim):

            paras.append(self.flatten[i])
            parae.append(self.mesh[i])
            xsmin.append(min(self.flatten[i]))
            xsmax.append(max(self.flatten[i]))
            xst.append((self.flatten[i]))
            div_paras.append((np.linspace(xsmin[i],paras[i][max_index],vdiv)).tolist()+(np.delete(np.linspace(paras[i][max_index],xsmax[i],vdiv),0)).tolist())
            div_paras2.append((lonum*(paras[i][max_index]-xsmin[i])+xsmin[i]).tolist()+(np.delete(hinum*(xsmax[i]-paras[i][max_index])+paras[i][max_index],0)).tolist())
        if dim==1:
            div_value=[]
            div_parae=div_paras2
            base_dis = interpolate.Rbf(parae[0],xe,function='linear', smooth=0)  #
            for j in range(len(div_paras[0])):
                div_value.append(max(0,base_dis(div_paras[0][j])))
            div_valus=div_value
            div_paras=div_paras2
        elif dim==2:
            div_parae=[[],[]]
            div_parae[0],div_parae[1]=np.meshgrid(div_paras[0],div_paras[1])
            base_dis = interpolate.Rbf(parae[0],parae[1],xe,function='linear', smooth=0)  #
            div_value=[]
            for j in range(len(div_parae[0])):
                value_tmp=[]
                for k in range(len(div_parae[0][j])):
                    value_tmp.append(max(0,base_dis(div_parae[0][j][k],div_parae[1][j][k])))
                div_value.append(value_tmp)
            div_parae[0],div_parae[1]=np.meshgrid(div_paras2[0],div_paras2[1])
            div_paras[0]=div_parae[0].flatten()
            div_paras[1]=div_parae[1].flatten()
            div_valus=(np.array(div_value)).flatten()
        elif dim==3:
            div_parae=[[],[],[]]
            div_parae[0],div_parae[1],div_parae[2]=np.meshgrid(div_paras[0],div_paras[1],div_paras[2])
            base_dis = interpolate.Rbf(parae[0],parae[1],parae[2],xe,function='linear', smooth=0)  #
            div_value=[]
            for j in range(len(div_parae[0])):
                value_tmp1=[]
                for k in range(len(div_parae[0][j])):
                    value_tmp2=[]
                    for h in range(len(div_parae[0][j][k])):
                        value_tmp2.append(max(0,base_dis(div_parae[0][j][k][h],div_parae[1][j][k][h],div_parae[2][j][k][h])))
                    value_tmp1.append(value_tmp2)
                div_value.append(value_tmp1)
            div_parae[0],div_parae[1],div_parae[2]=np.meshgrid(div_paras2[0],div_paras2[1],div_paras2[2])
            div_paras[0]=div_parae[0].flatten()
            div_paras[1]=div_parae[1].flatten()
            div_paras[2]=div_parae[2].flatten()
            div_valus=(np.array(div_value)).flatten()
        return max_value, max_index,div_parae,div_value,div_paras,div_valus
                                    
    def bunpu_contour3(self,vdiv):
        dim=self.dim
        div=self.div

        
        xs=(self.flatten[dim]).tolist()
        max_value = max(xs)
        max_index = xs.index(max_value)
        xe=(self.mesh[dim])
        xsmin=[]
        xsmax=[]
        xst=[]
        plistlo=[]
        plisthi=[]
        paras=[]
        parae=[]
        div_paras=[]
        minlo=[]
        minhi=[]
        conlistlo=[]
        conlisthi=[]
        if dim==1:
            tlistlo=[]
            tlisthi=[]
            paralistlo=[[]]
            paralisthi=[[]]
            plistlo.append([xe[i] for i in range(len(xe)) if i <= max_index])
            plisthi.append([xe[i] for i in range(len(xe)) if i >= max_index])
            parae.append(self.mesh[0])
            tlistlo.append([parae[0][i] for i in range(len(parae[0])) if i <= max_index])
            tlisthi.append([parae[0][i] for i in range(len(parae[0])) if i >= max_index])
            base_dislo00 = interpolate.Rbf(plistlo[0],tlistlo[0],function='linear', smooth=0) 
            base_dishi00 = interpolate.Rbf(plisthi[0],tlisthi[0],function='linear', smooth=0) 
            minlo.append(min(plistlo[0]))
            minhi.append(min(plisthi[0]))
            conlistlo.append((np.linspace(minlo[0],max_value,vdiv)).tolist())
            conlisthi.append((np.linspace(max_value,minhi[0],vdiv)).tolist())
            for j in range(vdiv):
                paralistlo[0].append(base_dislo00(conlistlo[0][j]))
                paralisthi[0].append(base_dishi00(conlisthi[0][j]))
        elif dim ==2:
            q,mod=divmod(max_index,div[0])
            plistlo.append([xe[q][i] for i in range(len(xe[q])) if i <= mod and (i>=1 and xe[q][i]>xe[q][i-1])])
            plisthi.append([xe[q][i] for i in range(len(xe[q])) if i >= mod and (i<=len(xe[q])-2 and xe[q][i]>xe[q][i+1])])
            plistlo.append([xe[i][mod] for i in range(len(xe)) if i <= q and (i>=1 and xe[i][mod]>xe[i-1][mod])])
            plisthi.append([xe[i][mod] for i in range(len(xe)) if i >= q and (i<=len(xe)-2 and xe[i][mod]>xe[i+1][mod])])
            parae.append(self.mesh[0])
            parae.append(self.mesh[1])
            tlistlo=[]
            tlisthi=[]
            
            for j in range(dim):
                minlo.append(min(plistlo[j]))
                minhi.append(min(plisthi[j]))
                conlistlo.append((np.linspace(minlo[j],max_value,vdiv)).tolist())
                conlisthi.append((np.linspace(max_value,minhi[j],vdiv)).tolist())
            tlistlo.append([parae[0][q][i] for i in range(len(parae[0][q])) if i <= mod and (i>=1 and xe[q][i]>xe[q][i-1])])
            tlisthi.append([parae[0][q][i] for i in range(len(parae[0][q])) if i >= mod and (i<=len(xe[q])-2 and xe[q][i]>xe[q][i+1])])
            tlistlo.append([parae[1][i][mod] for i in range(len(parae[1])) if i <= q and (i>=1 and xe[i][mod]>xe[i-1][mod])])
            tlisthi.append([parae[1][i][mod] for i in range(len(parae[1])) if i >= q and (i<=len(xe)-2 and xe[i][mod]>xe[i+1][mod])])

            base_dislo00 = interpolate.Rbf(plistlo[0],tlistlo[0],function='linear', smooth=0) 
            base_dishi00 = interpolate.Rbf(plisthi[0],tlisthi[0],function='linear', smooth=0) 
            base_dislo11 = interpolate.Rbf(plistlo[1],tlistlo[1],function='linear', smooth=0) 
            base_dishi11 = interpolate.Rbf(plisthi[1],tlisthi[1],function='linear', smooth=0) 
            paralistlo=[[],[]]
            paralisthi=[[],[]]
            for j in range(vdiv):
                paralistlo[0].append(base_dislo00(conlistlo[0][j]))
                paralisthi[0].append(base_dishi00(conlisthi[0][j]))
                paralistlo[1].append(base_dislo11(conlistlo[1][j]))
                paralisthi[1].append(base_dishi11(conlisthi[1][j]))

        elif dim ==3:
            q,mod=divmod(max_index,div[0])
            q1,q2=divmod(q,div[1])
            plistlo.append([xe[q1][q2][i] for i in range(len(xe[q1][q2])) if i <= mod])
            plisthi.append([xe[q1][q2][i] for i in range(len(xe[q1][q2])) if i >= mod])
            plistlo.append([xe[q1][i][mod] for i in range(len(xe[q1])) if i <= q2])
            plisthi.append([xe[q1][i][mod] for i in range(len(xe[q1])) if i >= q2])
            plistlo.append([xe[i][q2][mod] for i in range(len(xe)) if i <= q1])
            plisthi.append([xe[i][q2][mod] for i in range(len(xe)) if i >= q1])
            parae.append(self.mesh[0])
            parae.append(self.mesh[1])
            parae.append(self.mesh[2])
            tlistlo=[[],[],[]]
            tlisthi=[[],[],[]]
            for j in range(dim):
                minlo.append(min(plistlo[j]))
                minhi.append(min(plisthi[j]))
                conlistlo.append((np.linspace(minlo[j],max_value,vdiv)).tolist())
                conlisthi.append((np.linspace(max_value,minhi[j],vdiv)).tolist())
            tlistlo[0].append([parae[0][q1][q2][i] for i in range(len(parae[0][q1][q2])) if i <= mod])
            tlisthi[0].append([parae[0][q1][q2][i] for i in range(len(parae[0][q1][q2])) if i >= mod])
            tlistlo[1].append([parae[1][q1][i][mod] for i in range(len(parae[1][q1])) if i <= q2])
            tlisthi[1].append([parae[1][q1][i][mod] for i in range(len(parae[1][q1])) if i >= q2])
            tlistlo[2].append([parae[2][i][q2][mod] for i in range(len(parae[2])) if i <= q1])
            tlisthi[2].append([parae[2][i][q2][mod] for i in range(len(parae[2])) if i >= q1])
            base_dislo00 = interpolate.Rbf(plistlo[0],tlistlo[0],function='linear', smooth=0) 
            base_dishi00 = interpolate.Rbf(plisthi[0],tlisthi[0],function='linear', smooth=0) 
            base_dislo11 = interpolate.Rbf(plistlo[1],tlistlo[1],function='linear', smooth=0) 
            base_dishi11 = interpolate.Rbf(plisthi[1],tlisthi[1],function='linear', smooth=0) 
            base_dislo22 = interpolate.Rbf(plistlo[2],tlistlo[2],function='linear', smooth=0) 
            base_dishi22 = interpolate.Rbf(plisthi[2],tlisthi[2],function='linear', smooth=0) 
            paralistlo=[[],[],[]]
            paralisthi=[[],[],[]]
            for j in range(vdiv):
                paralistlo[0].append(base_dislo00(conlistlo[0][j]))
                paralisthi[0].append(base_dishi00(conlisthi[0][j]))
                paralistlo[1].append(base_dislo11(conlistlo[1][j]))
                paralisthi[1].append(base_dishi11(conlisthi[1][j]))
                paralistlo[2].append(base_dislo22(conlistlo[2][j]))
                paralisthi[2].append(base_dishi22(conlisthi[2][j]))

        for i in range(dim):
            paras.append(self.para[i].tolist())
            xst.append((self.para[i]).tolist())
            paralisthi[i].pop(0)
            div_paras.append(paralistlo[i]+paralisthi[i])
            
        if dim==1:
            div_value=[]
            div_parae=div_paras[0]
            base_dis = interpolate.Rbf(parae[0],xe,function='linear', smooth=0)  #
            for j in range(len(div_paras[0])):
                div_value.append(max(0,base_dis(div_paras[0][j])))
            div_valus=div_value
        elif dim==2:
            div_parae=[[],[]]
            div_parae[0],div_parae[1]=np.meshgrid(div_paras[0],div_paras[1])
            base_dis = interpolate.Rbf(parae[0],parae[1],xe,function='linear', smooth=0)  #
            div_value=[]
            for j in range(len(div_parae[0])):
                value_tmp=[]
                for k in range(len(div_parae[0][j])):
                    value_tmp.append(max(0,base_dis(div_parae[0][j][k],div_parae[1][j][k])))
                div_value.append(value_tmp)
            div_paras[0]=div_parae[0].flatten()
            div_paras[1]=div_parae[1].flatten()
            div_valus=(np.array(div_value)).flatten()
        elif dim==3:
            div_parae=[[],[],[]]
            div_parae[0],div_parae[1],div_parae[2]=np.meshgrid(div_paras[0],div_paras[1],div_paras[2])
            base_dis = interpolate.Rbf(parae[0],parae[1],parae[2],xe,function='linear', smooth=0)  #
            div_value=[]
            for j in range(len(div_paras[0])):
                value_tmp1=[]
                for k in range(len(div_paras[0][j])):
                    value_tmp2=[]
                    for h in range(len(div_paras[0][j][k])):
                        div_tmp2.append(max(0,base_dis(div_parae[0][j][k][h],div_parae[1][j][k][h],div_parae[2][j][k][h])))
                    value_tmp1.append(value_tmp2)
                div_value.append(value_tmp1)
            div_paras[0]=div_parae[0].flatten()
            div_paras[1]=div_parae[1].flatten()
            div_paras[2]=div_parae[2].flatten()
            div_valus=(np.array(div_value)).flatten()

        return max_value, max_index,div_parae,div_value,div_paras,div_valus
                                    

    def bunpu_filter(self,other,calc,divt0=[]):
        flag=0
        flag0=[]
        if type(other)==bunpu:
            m=[]
            xran=[]
            yran=[]
            ratio=12
            xdim=self.dim
            ydim=other.dim
            for i in range(xdim):
                if calc==0:
                    xran.append(self.xmax[i]-self.xmin[i])
                    yran.append(other.xmax[i]-other.xmin[i])
                    if yran[i]==0 or xran[i]/yran[i]>ratio:
                        flag0.append(1)
                    elif xran[i]==0 or yran[i]/xran[i]>ratio:
                        flag0.append(2)
                    else:
                        flag0.append(0)
                elif calc==1:
                    if self.xmax[i]+self.xmin[i] != 0:
                        xran.append((self.xmax[i]-self.xmin[i])/(self.xmax[i]+self.xmin[i]))
                    else:
                        xran.append((self.xmax[i]-self.xmin[i])/self.xmax[i])
                    if i < ydim:
                        if other.xmax[i]+other.xmin[i] != 0:
                            yran.append((other.xmax[0]-other.xmin[0])/(other.xmax[0]+other.xmin[0]))
                        else:
                            yran.append((other.xmax[0]-other.xmin[0])/other.xmax[0])
                    if yran[0]==0 or xran[i]/yran[0]>ratio:
                        flag0.append(1)
                    elif xran[i]==0 or yran[0]/xran[i]>ratio:
                        flag0.append(2)
                    else:
                        flag0.append(0)
            if flag0==[]:
                flag=0
            elif (xdim==1 and flag0[0]==1):
                flag=1
                if len(self.para[0])==1:
                    m.append(self.para[0])
                else:
                    m.append(self.mean())
            elif (xdim==1 and flag0[0]==2):
                flag=2
                if len(other.para[0])==1:
                    m.append(other.para[0])
                else:
                    m.append(other.mean())
            elif (xdim==2 and (flag0[0]==1 or flag0[1]==1)):
                flag=1
                for i in range(xdim):
                    if len(self.para[0])==1:
                        m.append(self.para[0])
                    else:
                        m.append(self.mean())
            elif (xdim==2 and (flag0[0]==2 or flag0[1]==2)):
                flag=2
                for i in range(xdim):
                    if len(other.para[0])==1:
                        m.append(other.para[0])
                    else:
                        m.append(other.mean())
            elif (xdim==3 and (flag0[0]==1 or flag0[1]==1 or flag0[2]==1)):
                flag=1
                for i in range(xdim):
                    if len(self.para[0])==1:
                        m.append(self.para[0])
                    else:
                        m.append(self.mean())
            elif (xdim==3 and (flag0[0]==2 or flag0[1]==2 or flag0[2]==2)):
                flag=2
                for i in range(xdim):
                    if len(other.para[0])==1:
                        m.append(other.para[0])
                    else:
                        m.append(other.mean())
            a0=bunpu()
            if flag==1:

                X=self.para
                X1=self.flatten
                X0=self.mesh
                #
                xf,px0,x0,dim,divx,dx,x0min,x0max=self.kaiseki()#
                if self.xmin==[] or self.xmax==[] or self.dim==0 or self.div==[] or self.dx==[]:
                    
                    self.xmin=x0min
                    self.xmax=x0max
                    self.dim=dim
                    self.div=divx
                    self.dx=dx
                    ymin=other.xmin
                    ymax=other.xmax
                else:
                    
                    x0min=self.xmin
                    x0max=self.xmax
                    dim=self.dim
                    dx=self.dx
                    ymin=other.xmin
                    ymax=other.xmax
                if divt0==[]:
                    for i in range(dim):
                        divt0.append(divx[i])

                
            elif flag==2:#

                X=other.para
                X1=other.flatten
                X0=other.mesh
                #
                xf,px0,x0,dim,divx,dx,x0min,x0max=other.kaiseki()#
                yf,py0,y0,dimy,divy,dy,ymin,ymax=self.kaiseki()#
                if other.xmin==[] or other.xmax==[] or other.dim==0 or other.div==[] or other.dx==[]:
                    
                    other.xmin=x0min
                    other.xmax=x0max
                    other.dim=dim
                    other.div=divx
                    other.dx=dx
                    ymin=self.xmin
                    ymax=self.xmax
                else:
                    
                    x0min=other.xmin
                    x0max=other.xmax
                    dim=other.dim
                    divx=other.div
                    dx=other.dx
                    ymin=self.xmin
                    ymax=self.xmax
                if divt0==[]:
                    for i in range(dim):
                        divt0.append(divx[i])
                
            elif flag==0:
                #x
                X=self.para
                X1=self.flatten
                X0=self.mesh
                #
                m=0
                xf,px0,x0,dim,divx,dx,x0min,x0max=self.kaiseki()#
                yf,py0,y0,dimy,divy,dy,ymin,ymax=other.kaiseki()#
                if self.xmin==[] or self.xmax==[] or other.xmin==[] or other.xmax==[] or self.dim==0 or self.div==[] or self.dx==[]:
                    
                    self.xmin=x0min
                    self.xmax=x0max
                    self.dim=dim
                    self.div=divx
                    self.dx=dx
                    other.xmin=ymin
                    other.xmax=ymax
                    other.div=divy
                    other.dx=dy
                else:
                    
                    x0min=self.xmin
                    x0max=self.xmax
                    ymin=other.xmin
                    ymax=other.xmax
                    dim=self.dim
                    dx=self.dx
                if divt0==[]:
                    for i in range(dim):
                        divt0.append(divx[i])
        else:
            X=self.para
            X1=self.flatten
            X0=self.mesh
            flag=3
            x0min=self.xmin
            x0max=self.xmax
            dim=self.dim
            divx=self.div
            dx=self.dx
            m=other
            ymin=m
            ymax=m

        return X,X1,X0,x0min,x0max,dim,dx,divx,m,ymin,ymax,divt0,flag
                
    def bunpu_percentile(self, percent, d):

        dim=self.dim
        x=self.mesh
        menseki=0
        lastmen=0
        dx=self.dx
        if dim==1:
            limitpara=0
            flag=0
            if d[0]>0:
                per=percent
            elif d[0]<0:
                per=1-percent
            for i in range(self.div[0]):
                lastmen=menseki
                menseki += x[1][i]*dx[0]
                if menseki>=per and flag==0:
                    flag=1
                    #limitpara=(x[0][i]*(menseki-per)+x[0][i-1]*(per-lastmen))/dx[0]
                    limitpara=x[0][i]
        elif dim==2:
            a=1
                    
        else:
            print('under constract')
            limitpara=0
        return limitpara


    def bunpu_percent(self, filename=[],pos1=[],pos2=[], dirc=[]):
        dim=self.dim
        x1=self.flatten
        x=self.mesh
        xmax=self.xmax
        xmin=self.xmin
        menseki=0
        lastmen=0
        dx=self.dx
        div=self.div
        if dim==1:
            flag=0
            flag1=0
            if pos1[0]==pos2[0] and dirc[0]>=0:
                flag=1
                flag1=1
                pos2[0]=xmax[0]+dx[0]
            elif pos1[0]==pos2[0] and dirc[0]<0:
                flag=3
                flag1=1
                pos2[0]=xmin[0]-dx[0]
            elif pos1[0]<pos2[0] and dirc[0]>=0:
                flag=1
            elif pos1[0]<pos2[0] and dirc[0]<0:
                flag=2
            elif pos1[0]>pos2[0] and dirc[0]>=0:
                flag=3
            elif pos1[0]>pos2[0] and dirc[0]<0:
                flag=4
            for i in range(div[0]):
                lastmen=menseki
                if (flag==1 and x[0][i]>=pos1[0] and x[0][i]<=pos2[0]) or (flag==2 and (x[0][i]<=pos1[0] or x[0][i]>=pos2[0])) or \
                  (flag==3 and x[0][i]<=pos1[0] and x[0][i]>=pos2[0]) or (flag==4 and (x[0][i]>=pos1[0] or x[0][i]<=pos2[0])):
                    menseki += x[1][i]*dx[0]
                elif i>=1:#posがi~i-1の間にある場合
                    if (flag==1 and x[0][i-1]>=pos1[0] and x[0][i-1]<=pos2[0]) or (flag==4 and (x[0][i-1]>=pos1[0] or x[0][i-1]<=pos2[0])):#x[i-1]<pos1<x[i]
                        menseki+=(x[1][i]*(pos1[0]-x[0][i-1])+x[1][i-1]*(x[0][i]-pos1[0]))
                    elif (flag==2 and (x[0][i-1]<=pos1[0] or x[0][i-1]>=pos2[0])) or (flag==3 and x[0][i-1]<=pos1[0] and x[0][i-1]>=pos2[0]):#x[i-1]<pos2[i]<x[i]
                        menseki+=(x[1][i]*(pos2[0]-x[0][i-1])+x[1][i-1]*(x[0][i]-pos2[0]))
                elif i<=div[0]-2:#posがi~i+1の間にある場合
                    if (flag==2 and (x[0][i+1]<=pos1[0] or x[0][i+1]>=pos2[0])) or (flag==3 and x[0][i+1]<=pos1[0] and x[0][i+1]>=pos2[0]):#x[i]<pos1<x[i+1]
                        menseki+=(x[1][i]*(x[0][i+1]-pos1[0])+x[1][i+1]*(pos1[0]-x[0][i]))
                    elif (flag==1 and x[0][i+1]>=pos1[0] and x[0][i+1]<=pos2[0]) or (flag==4 and (x[0][i+1]>=pos1[0] or x[0][i+1]<=pos2[0])):#x[i]<pos2<x[i+1]
                        menseki+=(x[1][i]*(x[0][i+1]-pos2[0])+x[1][i+1]*(pos2[0]-x[0][i]))
            if flag1==1:
                outbunpu=[[pos1[0],1],[pos1[0],1]]
            else:
                outbunpu=[[pos1[0],1],[pos2[0],1]]

            #グラフ
            #分布×ベクトル
            if filename==[]:
                filename=['output','unit']
                
            fg=plt.figure(figsize=(14,7))
            ax1 = fg.add_subplot(111)
            ln1 = ax1.plot(x[0],x[1],'b-',label=filename[0])
            
            ax2 = ax1.twinx()
            if flag1==0:
                ln2 =ax2.plot([xmin[0],pos1[0],pos1[0],pos2[0],pos2[0],xmax[0]],[0,0,1,1,0,0],'r-', label="border")
            elif flag1==1:
                ln2 =ax2.plot([xmin[0],pos1[0],pos1[0],pos1[0],pos1[0],xmax[0]],[0,0,1,1,0,0],'r-', label="border")
            h1, l1 = ax1.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax1.legend(h1+h2, l1+l2, loc='upper right')
            
            fp = FontProperties(fname=r'C:\Windows\Fonts\msgothic.ttc',size=24)
            ax1.set_xlabel(u''+filename[1], fontproperties=fp)
            ax1.set_ylabel(r'probability')
            ax1.grid(True)
            ax2.set_ylabel(r'border')
            
            plt.title(u''+filename[0], fontproperties=fp) 
            plt.xlabel(u''+filename[1], fontproperties=fp)
            posx=pos1[0]
            posy1=0.5

            plt.text(posx,posy1,'menseki='+str(menseki),fontproperties=fp)

            outfilename=filename[0]
            outgraphname=outfilename+'.png'
            plt.savefig(outgraphname)
                
        elif dim==2:
            #
            #
            inbunpu=[[],[]]
            outbunpu=[[],[]]
            poslg1=pos1[0]*dirc[0]+pos1[1]*dirc[1]
            poslg2=pos2[0]*dirc[0]+pos2[1]*dirc[1]
            for i in range(len(x1[0])):
                xlg=x1[0][i]*dirc[0]+x1[1][i]*dirc[1]
                if (poslg1==poslg2 and xlg>=poslg1) or (poslg1>poslg2 and xlg>=poslg2 and xlg<=poslg1) or (poslg1<poslg2 and xlg<=poslg2 and xlg>=poslg1):
                    menseki+=x1[2][i]*dx[0]*dx[1]
                    inbunpu[0].append(((x1[0][i]-dirc[0]*xlg/(dirc[0]**2+dirc[1]**2))**2+(x1[1][i]-dirc[1]*xlg/(dirc[0]**2+dirc[1]**2))**2)**0.5)
                    inbunpu[1].append(x1[2][i])

            if menseki==0:
                outbunpu[0].append(0)
                outbunpu[1].append(0)
            else:
                maxb=max(inbunpu[0])
                minb=min(inbunpu[0])
                nin=len(inbunpu[0])
                divout=int((maxb-minb)//min(dx[0],dx[1]))
                if divout>1:
                    dxb=(maxb-minb)/(divout-1)
                    outbunpu[0]=(np.linspace(maxb,minb,divout)).tolist()
                    #
                    j=0
                    lastp=0
                    fg=0
                    for j in range(divout):
                        num=0
                        pin=0
                        for i in range(len(inbunpu[0])):
                            if inbunpu[0][i]>=outbunpu[0][j]-dxb/2 and inbunpu[0][i]<outbunpu[0][j]+dxb/2:
                                num+=1
                                pin+=inbunpu[1][i]
                        if num==0 and j>0 and j<divout-1:
                            outbunpu[1].append(lastp)
                            fg=1
                        elif num!=0:
                            outbunpu[1].append(pin/num)        
                            lastp=pin/num
                            if fg==1:
                                lastoutbunpu=outbunpu[1][j-1]
                                outbunpu[1][j-1]=(lastoutbunpu+lastp)/2
                                fg=0
                        else:
                            outbunpu[1].append(0)        
                            
                else:
                    pin=0
                    xin=0
                    for i in range(nin):
                        xin+=inbunpu[0][i]
                        pin+=inbunpu[1][i]
                    outbunpu[0].append(xin/nin)
                    outbunpu[1].append(pin/nin)
        else:
            print('under constract')
            limitpara=0
            menseki=0
            outbunpu=0
        return menseki,outbunpu

    
        
    def __exit__(self, exception_type, exception_value, traceback):
        print('END')
