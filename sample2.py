# -*- coding: utf-8 -*-


"""
ストレスストレングス
作動頻度と生涯作動時間から生涯作動回数を求めて、耐久条件を求める
"""
from bunpu import *
import matplotlib.pyplot as plt
from scipy.stats import weibull_min
import numpy as np

div=[200]
#目標故障率
tgtf=0.002
#ワイブル分布の何%タイルを耐久条件とするか
taikyutile=0.6
#ワイブル形状パラメータ
c=6
#ワイブル分布の初期位置（故障率が目標より高ければ、これより右にシフト）
scale=1.0
#作動頻度、ここではダミー分布を使うが、実際には個体の作動頻度（作動回数/時間など）データから母集団の作動頻度分布を作成
hindo=bunpu()
hindo.bunpu_gene([0.0],[4.1],[0.4],[0.5],div,'hindo')

#生涯作動時間、ここではダミー分布を使うが、実際には個体の生涯作動時間データから母集団の生涯作動時間分布を作成
lifetime=bunpu()
lifetime.bunpu_gene([0.0],[12000],[4200],[2000],div,'lifetime')
#生涯作動回数
lifecnt=bunpu()
lifecnt=lifetime*hindo

lifetime.bunpu_graph('lifetime')
hindo.bunpu_graph('hindo')
lifecnt.bunpu_graph('lifecnt')


x1=np.linspace(lifecnt.xmin[0],lifecnt.xmax[0],div[0])
x2=np.linspace(weibull_min.ppf(0.0001, c),weibull_min.ppf(0.9999, c),div[0])#パーセントポイント
y2=weibull_min.pdf(x2, c)#確率分布
y3=weibull_min.cdf(x2, c)#累積分布
fy = interpolate.Rbf(x2, y3, function='linear', smooth=0)
dx=(x1[div[0]-1]-x1[0])/div[0]
menseki=1
while menseki>tgtf:
    x3=[]
    #ワイブル累積分布（ストレングス）
    w=[]
    #生涯作動回数&ワイブル累積
    z=[]
    menseki=0
    for i in range(div[0]):
        x3=x2[0]+(x2[i]-x2[0])/scale
        w.append(fy(x3))
        z.append(w[i]*lifecnt.mesh[1][i])
        menseki+=z[i]*dx
    scale *= 1.2
scale=scale/1.2
taikyucnt=x1[div[0]-1]*scale*weibull_min.ppf(taikyutile, c)/weibull_min.ppf(0.9999, c)

#グラフ
fg=plt.figure(figsize=(14,7))
ax1 = fg.add_subplot(111)

ln1 = ax1.plot(lifecnt.mesh[0],lifecnt.mesh[1],'b-',label="stress")
ln1 = ax1.plot(lifecnt.mesh[0],z,'r-', label="failure")

ax2 = ax1.twinx()
ln2 =ax2.plot(x1, w, 'g-',label="wibull")

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc='upper right')

fp = FontProperties(fname=r'C:\Windows\Fonts\msgothic.ttc',size=24)
ax1.set_xlabel(u''+'history', fontproperties=fp)
ax1.set_ylabel(r'stress failure')
ax1.grid(True)
ax2.set_ylabel(r'wibull')

plt.title(u''+'stress strength', fontproperties=fp) 
plt.xlabel(u''+'history', fontproperties=fp)
posx=dx*0.3*div[0]
posy1=0.9*max(w)
posy2=0.8*max(w)
posy3=0.7*max(w)
plt.text(posx,posy1,'menseki='+str(menseki),fontproperties=fp)
plt.text(posx,posy2,'percentile='+str(taikyutile),fontproperties=fp)
plt.text(posx,posy3,'taikyu kaisu='+str(taikyucnt),fontproperties=fp)
outfilename=['stressstrength']
outgraphname=outfilename[0]+'.png'
plt.savefig(outgraphname)

#対数表示
ax1.semilogy()

outfilename=['stressstrength_log']
outgraphname=outfilename[0]+'.png'
plt.savefig(outgraphname)
            
