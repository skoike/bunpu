# 分布演算（確率分布ベクトル解析）
# Calculation by distribution (Probability distribution vector analysis)


　既存の様々な数学の数値パラメータを分布に置換えて確率演算を包含する分布演算を実現する。
既存のモンテカルロシミュレーションや確率過程では全てのパスを演算できなければ正しい分布が求まらないが、それはランダムウォークなどの特殊な分布を除いて困難。
分布演算はそもそもパスを求める必要がなく、ヒストグラムなどから求めた任意の分布から最も確からしい結果の分布を求めることができる。

　分布演算は、モンテカルロシミュレーションや確率過程の様に、個別のパスを演算しないが、時系列の分布間に相関関係を持たせることが可能で、速度vと位置xの
時系列の相関関係で様々な運動の特徴を表すことができる。例えば、ランダムウォークは、xが正相関の加算、vが独立加算。マルコフ連鎖はx,v正相関の加算、
ホワイトノイズはxが独立加算、フィードバック制御は逆相関の加算で模擬できると考えている。

　今回そのツールの一部を公開するが、現在はベクトルの四則演算や単純な時系列積分などができる程度である。今後様々な分野の方に協力して頂き、開発を進めることができれば、
様々な微分方程式を含む数学全体に拡張できると考えている。このことは様々なバラツキを扱う物理解析技術や社会現象解析に改善をもたらし、理解を深めるスコープとして活用できる可能性がある。

現時点でこのツールは、以下のことが可能である。

- データファイルの特定列から抽出したデータのヒストグラムと分布を生成
- 範囲や平均値、標準偏差を指定して分布を生成
- 前記生成された分布の四則演算や時系列積分（相関係数に応じた補正が可能）
- 分布のグラフ表示
- 分布要素をファイル出力
- 以上の1次元から3次元のベクトル分布処理

ライセンスは、このソフトをそのまま利用するだけであればフリー、ソフトの変更や参考にして作成したものの配布や商用利用する場合は知財権利と著作権にご配慮ください。

 Realize the calculation by distribution including probability operations by replacing various existing mathematical parameters with distributions.
In the existing Monte Carlo simulation and stochastic process, the correct distribution cannot be obtained unless all paths can be calculated,
but it is difficult except for special distributions such as random walks.
In the calculation by distribution, it is not necessary to obtain the path in the first place,
and the distribution of the most probable result can be obtained from an arbitrary distribution obtained from a histogram or the like.

 the calculation by distribution does not calculate individual paths like Monte Carlo simulation and stochastic process,
 but it is possible to have a correlation between time series distributions,
 and the velocity v and position x Various movement characteristics can be represented by time-series correlation.
 For example, I think that random walk can be simulated by adding x positive correlation and v independent correlation,
 Markov chains adding x, v positive correlations, white noise adding independent x and feedback control adding inverse correlation.

  Currently, the method is limited to four arithmetic operations and simple integration of vectors(time series analysis),
 but I believe that it can be extended to the entire mathematics including various differential equations by advancing development.
 This may bring improvements to various physical analysis technology and social phenomenon analysis that handle various variations
 and can be used as a scope to deepen understanding.

The method is to construct a calculation system in which the distributed data is distributed in a histogram and the distribution is used as a calculation element.
At the moment this tool can:

- Generate a distribution from the histogram of the measured data
- Generate distribution by specifying range, mean, and standard deviation
- Four arithmetic operations between generated distributions and time series integral
- Graph display of distribution
- Output distribution elements to a file
- 1D to 3D vector distribution processing above all

The license is free if you just use it as it is, If you change or use this soft as reference for a distribution or commercial purposes,
you should respond to your obligations for the intellectual property rights and the copyright of this software.


## 論文
## literature

- 市場走行データを活用した設計方法,Toyota Technical Review 2018/5 Vol64.p95
  - Design Method Using Real-World Vehicle Data,Toyota Technical Review 2018/9 Vol64.p96
  - https://shop.ohmsha.co.jp/shopdetail/000000005659/
  - アクチュエータの走行距離あたりの作動頻度分布と生涯走行距離分布の積から生涯作動回数のストレスストレングス分布を作成して、耐久条件を求める.
  - レーダクルーズのブレーキ制御頻度分布と日当たり走行分布とお客様の入庫日数分布から、入庫時に記録が上書きされないメモリーと制御記録要件を設計する.

- ビックデータを活用した制御リスク設計,自動車技術会20年秋季大会学術講演会
  - Risk design of control system utilizing big data,Society of Automotive Engineers of Japan
  - https://www.jstage.jst.go.jp/article/jsaeronbun/52/1/52_20214041/_article/-char/ja
  - レーダー認識距離分布と車線変更制御完了時間から安全な車線変更制御設計を行う.

- 確率分布ベクトル解析について、情報処理学会第83回全国大会（2021 3/18）
  - https://www.ipsj.or.jp/event/taikai/83/ipsj_web2021/data/pdf/1B-03.html
  - モンテカルロシミュレーションと分布演算、それぞれを使って放物線運動を行う飛距離の分布を求めて比較、正確な分布を求めるために分布演算が必要であることを示す.

## URL
- http://www.na.rim.or.jp/~syn/kakuritsu.html


##プログラム言語とインポートするモジュール
##import module

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
11. d.bunpu_file('output')
12. print(d.mesh)



各行の機能 

1. 分布クラスツールのインポート 
2. 分布クラスからインスタンスの定義
3. ファイルデータから2次元分布生成する。引数:'データファイル名','分布出力名',無視する行,取り込む列,分割数,取り込む範囲の制限,カーネル分布のバンド幅
4. 分布クラスからインスタンスの定義
5. 範囲や平均を指定して2次元分布を生成、引数:最小値,最大値,平均値,標準偏差,分割数,分布出力名
6. 分布クラスからインスタンスの定義
7. 範囲や平均を指定して1次元分布を生成、引数:最小値,最大値,平均値,標準偏差,分割数,分布出力名
8. 分布クラスからインスタンスの定義
9. 分布演算（積商は多次元×1次元＝多次元、和差は多次元＋多次元＝多次元） 
10. グラフ表示 
11. 分布のパラメータ、確率値をファイル出力
12. 分布のパラメータと確率値（メソッド）をコンソールに表示
    
    
9の演算後に表示されるものの意味
- bunpu+bunpu:通常の分布演算
- bunpu+lean:分布範囲の差が大きい時に分布として演算されなかった場合
- bunpu+vector:分布とベクトルの演算
- 上記表示に続く数字:演算結果の素の分布面積（体積、超体積）、相関係数0での分布演算結果の面積は必ず1に近い数字になる、1からのズレ2割以上大きい場合は結果が信頼できない、相関係数が0以外の数字を指定した場合１以下になる。


サンプル

- sample0.py　最小値、最大値、平均値、標準偏差を指定して分布作成、その分布を四則演算、1次元～3次元
- sample1.py　ファイルからデータを読み込んで分布作成、その分布を四則演算、1次元～3次元
- sample2.py　作動頻度分布×生涯寿命分布＝生涯作動回数分布をもとめ、ストレスストレングスで目標故障率を満足する耐久回数を求める
- sample3.py　拡散方程式の解、
- sample4.py　運動方程式の解、3次元放物線運動




課題
- python3.8.5だと3D表示がおかしい（python3.6.9だと正常）→matplotlib2.2.2だと正常、新しいとsavefigがおかしい、サンプルのsavefigを表示だけにする


## ライセンス

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

## License

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




