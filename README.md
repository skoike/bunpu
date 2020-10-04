
#"分布演算の環境を提供する"

#Providing an environment for distribution calculation


　ビックデータを解析する中で従来の数値演算では大きなバラツキが存在する現実世界をモデル化するには精度限界があることに気づき、それを対策する方法の具体例を様々な論文にて発表してきた。
その方法は25年以上前から様々な応用例を提案してきたが、いまだにその必要性が理解されていないので、今回、会社の許可をもらって誰にでも使えるツールとして広く使ってもらえることでその価値の理解を促進しようと考えた。

　従来の数値演算は、バラツキのある複数データの演算結果から分布を求めても再現性のある分布が得られない、それはデータ数を増やしても改善しない。それはモンテカルロシミュレーションが正しい分布を提供できないことで
知られている。その理由は、数値演算からバラツキを把握しようとしても、多少でも独立性があるパラメータ間の演算では情報欠損が大きいので、非現実的な数の冗長なデータが無いとまともな精度の分布が得られないからである。
対策として、少数のデータでも、それぞれのパラメータのヒストグラムから分布を作成して、分布演算を行うことで、数値演算の情報欠損が補間されて精度の高い演算が可能となる。
　今回そのツールの一部を公開するが、現在はベクトルの四則演算や単純な積分ができる程度である。今後様々な分野の方に協力して頂き、開発を進めることができれば、様々な微分方程式を含む数学全体に拡張できると考えている。
そのことはナビエストークスをはじめとする、様々な基本的な物理解析技術に改善をもたらし、理解を深めるスコープとして活用できる可能性がある。

このツールの目的は、バラツキがあるデータからヒストグラムで分布を抽出して、その分布を演算要素とする演算体系を構築することである。現時点でこのツールは、以下のことが可能である。
ライセンスは、そのまま利用するだけであればフリー、変更する場合は公開前提で提供頂ければフリー、変更や参考としたものの配布や商用利用する場合は知財権利において有償とします。

###ファイルの特定列から抽出したデータのヒストグラムと分布を生成
###範囲や平均値、標準偏差を指定して分布を生成
###四則演算（相関係数に応じた補正が可能）
###分布のグラフ表示
###分布要素をファイル出力
###以上の1次元から3次元のベクトル分布処理

　While analyzing big data, I noticed that there is a limit to the accuracy of modeling the real world with conventional numerical operations, and as a real example of how to deal with it, I have been published in various papers in Japan.I have been proposing various application examples of this method for more than 20 years, but Few fellow understand its importance, so I promote understanding of its value by making it widely used as a tool that anyone can use this time with the permisson of my company. I thought about it.

In the conventional numerical calculation, a reproducible distribution cannot be obtained, if the distribution is obtained from the calculation results of a plurality of various data, which does not improve even if the number of data is increased. That is because the Monte Carlo simulation cannot provide the correct distribution, Are known.
The reason is that even if we try to grasp the variation from the numerical calculation, the information loss is large in the calculation between the parameters that have some independence, so that a highly accurate distribution cannot be obtained without an unrealistic number of data.
As a countermeasure, even for a small amount of data, by creating a distribution from the histogram of each parameter and performing the distribution calculation, the information loss of the numerical calculation is interpolated and the calculation with high accuracy becomes possible.

 Currently, the method is limited to four arithmetic operations and simple integration of vectors, but I believe that it can be extended to the entire mathematics including various differential equations by advancing development.This may bring improvements to various basic physical analysis techniques such as Navier-Stokes and can be used as a scope to deepen understanding.

The method is to construct a calculation system in which the distributed data is distributed in a histogram and the distribution is used as a calculation element.
At the moment this tool can:

### Generate a distribution from the histogram of the measured data
### Generate distribution by specifying range, mean, and standard deviation
### Four arithmetic operations
### Graph display of distribution
### Output distribution elements to a file
### 1D to 3D vector distribution processing above

##


##論文
##literature

ビックデータを活用した制御リスク設計,自動車技術会20年秋季大会学術講演会
Risk design of control system utilizing big data,Society of Automotive Engineers of Japan 

##プログラム言語とインポートするモジュール
python
glob
re
csv
os
pandas
copy
scipy
numpy
matplotlib
mpl_toolkits
sklearn

##ツールの利用方法
##How to use

from bunpu import *   
a=bunpu()
a.bunpu_data('input','output',1,[1,2],0,[20,10],[0],1.5)
b=bunpu()
b.bunpu_gene([-1.0,2.0],[0.5,4.0],[-0.1,2.8],[0.2,0.3],[20,20],'b')
c=bunpu()
c.bunpu_gene([2.0],[4.0],[2.8],[0.3],[20],'c')
d=bunpu()
d=b*c+a
d.bunpu_graph()

    各行の機能
    1 ツールのインポート
    2,4,6,8 分布メソッドの定義
    3 ファイルデータから分布生成、引数:'データファイル名','分布出力名',無視する行,取り込む列,分割数,取り込む範囲の制限,カーネル分布のバンド幅
    5,7 範囲や平均を指定して分布を生成、引数:最小値,最大値,平均値,標準偏差,分割数,分布出力名
    9 分布演算（積商は多次元*1次元、和差は多次元*多次元）
    10 グラフ表示
    


##ライセンス
##License

このソフトウェアをそのままの複製として利用する場合、本ソフトウェアの利用を無償で許可します。

このソフトウェアは未完成で、改善の提案や機能拡張の協力を求めています、このソフトの改善や協力の為にに、変更、追加、結合、移植を含む派生を、
利用可能な情報とともに、公開を前提として、このソフトの作者または著作権者にその情報提供を行うことで、その利用を許可します。

このソフトを利用・参考にする場合は、このソフトの著作権と特許出願（PCT/JP2020/034566とその分割、関連出願）における権利を尊重ください。
このソフトウェアの一部分を利用または参考にして、変更、追加、結合、継承や移植を含む派生を、非公開で配布、または商用利用する場合は作者に相談してください。

ソフトウェアは、未完成で、何らの保証もなく提供されます。
ここでいう保証とは、商品性、特定の目的への適合性、および権利非侵害についての保証も含みますが、それに限定されるものではありません。 
このソフト作者または著作権者は、契約行為、不法行為、またはそれ以外であろうと、ソフトウェアに起因または関連し、あるいはソフトウェアの使用または
その他の扱いによって生じる一切の請求、損害、その他の義務について何らの責任も負わないものとします。

以上の表示および本許諾表示を、ソフトウェアのすべての複製または部分の利用または分布処理を参考とする場合に、作成される著作物に記載するものとします。


Permission is hereby granted, free of charge, to any person obtaining a exact copy of this software and associated documentation files (the "Software"), 
to deal in the Software with restriction.

This software is incomplete and we are seeking suggestions for improvement and cooperation in enhancements.
Anyone who creates and uses, including changes, additions, mergers, and translation, to improve or cooperate with this software, 
by providing it with the information available to the creator or copyright holder of this software on the premise of publication. Is possible.

When using or referring to this software, please correspond the copyright of this software and the rights in patent applications(PCT/JP2020/034566 and divisional other).
Please contact with the author if you want to use or refer to a part of this software and distribute it privately or use it for commercial purposes.

The software is incomplete and is provided without warranty.Warranties here include, but are not limited to, warranties of merchantability, 
fitness for a particular purpose, and non-infringement.
The author or copyright holder of this software, whether contractual, tort, or otherwise, is due to or related to the software, or uses or uses the software.
We shall not be liable for any claims, damages or other obligations arising from any other dealings.

The above copyright notice and this permission notice shall be included in all copies, portions or reference of the software related to dstribution .


