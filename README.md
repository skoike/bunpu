# 分布演算（確率分布ベクトル解析）
# Distribution calculation(Probability distribution vector analysis)


　ビックデータを解析する中で従来の数値演算では大きなバラツキが存在する現実世界をモデル化するには精度限界があることに気づき、それを対策する方法として分布演算を提案してきた。
分布演算は、機能を拡張することでベクトル解析と確率統計を包含する演算体系になると考えている。従来の方法ではできなかった両者のはざまにある確率演算が可能となる。
実際に今まで様々な応用例を提案してきたが、説明するだけではその効果が十分に理解されていないので、誰にでも使えるツールとして広く使ってもらえることでその価値の理解を促進しようと考えた。

　ここで言う分布演算とは、数値に代えて確率分布を演算要素とした、演算方法である。従来から連続した関数表現の確率分布を畳込み積分によって演算する方法が提案されているが、
今回提案する方法によって、計測された実データのヒストグラムを演算対象として容易に四則演算ができるので、その応用範囲は広がる。
現在、物理現象や社会現象を推定する多くの方程式はパラメータを数値演算されているが、多くのパラメータはバラツキを持っており、その範囲が狭ければよいが、ある程度の範囲がある場合は
その中でやみくもに答えを求めていることになる。従って分布演算に拡張することで精度向上につながる。
従来から使われているモンテカルロシミュレーションはバラツキを考慮することができるが、同様に演算結果は無視できないレベルで誤差を持っており、正しい分布を出力することができない。
微分方程式を分布演算に拡張することで、初めて正確な演算結果の分布を求める事ができる。

　従来の数値演算では、相互に独立なバラツキのあるデータ間の演算結果から分布を求めても再現性のある分布が得られない、それはデータ数を増やしても改善しない。
その理由は、数値演算からバラツキを把握しようとしても、多少でも独立性があるパラメータ間の演算では、分布を形成するために必要な組合せ情報に対して欠損している情報量の割合が大きいので、
偶然によって分布の形状が大きく影響を受ける可能性が高い、データ量や演算量が増えることで必要な組合せ情報量が指数的に増加するので、データを増やしても正しい分布は得られない。
対策として、少数のデータでも、それぞれのパラメータのヒストグラムから分布を作成して、分布演算を行うことで、数値演算の情報欠損が補間されて精度の高い演算が可能となる。

　今回そのツールの一部を公開するが、現在はベクトルの四則演算や単純な積分ができる程度である。今後様々な分野の方に協力して頂き、開発を進めることができれば、
様々な微分方程式を含む数学全体に拡張できると考えている。このことは様々なバラツキを扱う物理解析技術や社会現象解析に改善をもたらし、理解を深めるスコープとして活用できる可能性がある。

このツールの目的は、バラツキがあるデータからヒストグラムで分布を抽出して、その分布を演算要素とする演算体系を構築することである。現時点でこのツールは、以下のことが可能である。

- データファイルの特定列から抽出したデータのヒストグラムと分布を生成
- 範囲や平均値、標準偏差を指定して分布を生成
- 前記生成された分布の四則演算（相関係数に応じた補正が可能）
- 分布のグラフ表示
- 分布要素をファイル出力
- 以上の1次元から3次元のベクトル分布処理

ライセンスは、このソフトをそのまま利用するだけであればフリー、ソフトの変更や参考にして作成したものの配布や商用利用する場合は知財権利と著作権にご配慮ください。


　While analyzing big data, I noticed that there is a limit to the accuracy of modeling the real world with conventional numerical operations,
and I have proposed distribution calculation to deal with it.
I believe that distribution arithmetic will become an arithmetic system that includes vector analysis and probability statistics by expanding its functions.
It is possible to perform probability calculations between the two, which was not possible with conventional methods.
I have been proposing various application examples of this method, but Few fellow understand its advantageous effect,
so I promote understanding of its value by making it widely used as a tool that anyone can use. I thought about it. 

 Distribution calculation is a method of solving solutions such as four arithmetic operations
and differential equations using a probability distribution as an arithmetic element instead of a numerical value.
Conventionally, it has been proposed to calculate the probability distribution of continuous function representation by convolution integral.
But,by using the method proposed this time it is possible to perform calculation using the histogram of the measured actual data as the calculation element,
the possibilities are endless.
Currently, many equations for estimating physical and social phenomena have their parameters calculated numerically,
but many parameters have variations, If there is a certain margin of error, in that process, anyone should be blindly seeking an answer.
So to extend to distribution calculation leads to improved accuracy.
The Monte Carlo simulation that has been used conventionally can take into consideration the variation,
but for same reason, the calculation result has an error at a non-negligible level, and the correct distribution cannot be output.
It is possible to obtain an accurate distribution of calculation results for the first time by performing a distribution calculation.

In the conventional numerical calculation, a reproducible distribution cannot be obtained, if the distribution is obtained
from the calculation results of a plurality of various data, which does not improve even if the number of data is increased.
The reason is that even if we try to grasp the variation from the numerical calculation, in operations between independent parameters,
the ratio of the amount of missing information to the combinational information required to form the distribution is large,
so there is a high possibility that the shape of the distribution will be greatly affected by chance.
As the amount of data and the amount of calculation increase, the amount of combination information required increases exponentially,
so even if the amount of data is increased, the correct distribution cannot be obtained.

As a answer, even for a small amount of data, by creating a distribution from the histogram of each parameter and performing the distribution calculation,
the information loss of the numerical calculation is interpolated and the calculation with high accuracy becomes possible.
 Currently, the method is limited to four arithmetic operations and simple integration of vectors,
 but I believe that it can be extended to the entire mathematics including various differential equations by advancing development.
 This may bring improvements to various physical analysis technology and social phenomenon analysis that handle various variations
 and can be used as a scope to deepen understanding.

The method is to construct a calculation system in which the distributed data is distributed in a histogram and the distribution is used as a calculation element.
At the moment this tool can:

- Generate a distribution from the histogram of the measured data
- Generate distribution by specifying range, mean, and standard deviation
- Four arithmetic operations between generated distributions
- Graph display of distribution
- Output distribution elements to a file
- 1D to 3D vector distribution processing above all

The license is free if you just use it as it is, If you change or use this soft as reference for a distribution or commercial purposes,
you should respond to your obligations for the intellectual property rights and the copyright of this software.


## 論文
## literature

- ビックデータを活用した制御リスク設計,自動車技術会20年秋季大会学術講演会
- Risk design of control system utilizing big data,Society of Automotive Engineers of Japan 

## URL

- http://www.na.rim.or.jp/~syn/kakuritsu.html


##プログラム言語とインポートするモジュール


- python3.6
- pandas0.25.2
- scipy1.2.0
- numpy1.16.5
- matplotlib3.0.2
- scikit-learn0.21.3


## ツールの利用方法
## How to use


1. from bunpu import *   
2. a=bunpu() 
3. a.bunpu_data('input','output',1,[1,2],0,[20,10],[0],1.5) 
4. b=bunpu() 
5. b.bunpu_gene([-1.0,2.0],[0.5,4.0],[-0.1,2.8],[0.2,0.3],[20,20],'b') 
6. c=bunpu() 
7. c.bunpu_gene([2.0],[4.0],[2.8],[0.3],[20],'c') 
8. d=bunpu() 
9. d=b*c+a 
10. d.bunpu_graph() 




各行の機能 

1. ツールのインポート 
2. 分布メソッドの定義 
3. ファイルデータから2次元分布生成する。引数:'データファイル名','分布出力名',無視する行,取り込む列,分割数,取り込む範囲の制限,カーネル分布のバンド幅
4. 分布メソッドの定義 
5. 範囲や平均を指定して2次元分布を生成、引数:最小値,最大値,平均値,標準偏差,分割数,分布出力名
6. 分布メソッドの定義
7. 範囲や平均を指定して1次元分布を生成、引数:最小値,最大値,平均値,標準偏差,分割数,分布出力名
8. 分布メソッドの定義 
9. 分布演算（積商は多次元×1次元＝多次元、和差は多次元＋多次元＝多次元） 
10. グラフ表示 
    
    
9の演算後に表示されるものの意味
- bunpu+bunpu:通常の分布演算
- bunpu+lean:分布範囲の差が大きい時に分布として演算されなかった場合
- bunpu+vector:分布とベクトルの演算
- 上記表示に続く数字:演算結果の素の分布面積（体積、超体積）、相関係数0での分布演算結果の面積は必ず1に近い数字になる、1からのズレ2割以上大きい場合は結果が信頼できない

課題
- python3.8.5だと3D表示がおかしい（python3.6.9だと正常）→matplotlib2.2.2だと正常、新しいとsavefigがおかしい、サンプルのsavefigを表示だけにする


## ライセンス
## License
© 2020 Shin Koike  bunpu@a1.rim.or.jp

このソフトウェアをそのままの複製として利用する場合、本ソフトウェアおよび今後作成されるものを含めたそのブランチの利用を無償で許可します。

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
its branches and associated documentation files (the "Software"), to deal in the Software with restriction.

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




