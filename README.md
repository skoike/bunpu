# ＜バラツキの対処法＞
# ＜How to deal with bariance＞

#  分布演算によりバラツキを厳密に扱う設計、解析、演算結果の判断方法、捉え方の提案

# Mathematical method for handling dispersion
#  Proposal of new design, new analysis, and understanding using Calculation by distribution

2023年2月18日技術評論社から発売の“バラツキの対処法”で説明した演算を行うためのソフトウェアです。
書籍ではこのソフトウェアを電卓と呼びましたが、これは一般的な電卓で扱う数字を分布に置換えて、
分布どうしをを演算するための電卓です。

　このソフトウェアはそのままの複製を学習や研究を目的として利用する場合に限り、フリーに使ってもらえます。
それ以外の以下のケースなどは、ライセンス記述にあるアドレスに相談ください。
個別のニーズに応じた技術サポートを、可能な範囲で2023年5月以降から有償で対応できるようにします。

・本技術を利用したモノやサービスを産業活動（商用）として行う場合
・本技術の関数スクリプト(python)を入手して自由度の高い活用を行う場合。
・本技術を参考にして類似のソフトウェアを開発・配布する場合。
・本技術を活用する為のコンサルティング、説明、講演などが必要な場合。
・ユーザーニーズにあわせて本技術を応用するためのツールカスタマイズ。


Windows10で開発していましたが、Windows11ではうまく動かなかったので、
それぞれの環境でコンパイルしたものを用意しました。

Windows10用→bunpu_win10.exe
Windows11用→bunpu_win11.exe

その他csvファイルが多数ありますが、これは上記ツールで練習用に使うダミーデータです。
使い方は“バラツキの対処法”を参照ください。



現時点でこのツールは、以下のことが可能である。

- データファイルの特定列から抽出したデータのヒストグラムと分布を生成
- 範囲や平均値、標準偏差を指定して分布を生成
- 前記生成された分布の四則演算や時系列積分（シミュレーション）
- 分布のグラフ表示
- 分布要素をファイル出力
- 以上の1次元から3次元のベクトル分布処理

ライセンスは、このソフトをそのまま利用するだけであればフリー、ソフトの変更や参考にして作成したものの配布や商用利用する場合は知財権利と著作権にご配慮ください。


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
  - https://shop.ohmsha.co.jp/shopdetail/000000005824/Z-EY/page1/order/
  - ブレーキアクチュエータの走行距離あたりの作動頻度分布(回/km)と生涯走行距離分布(km)の積から生涯作動回数の分布(回)を求め、ストレスストレングス分布(回)を作成して、耐久条件を求める.
  - レーダクルーズのブレーキ制御頻度分布(回/km)と日当たり走行分布(km/day)とお客様の入庫日数分布(day)から、入庫時に記録が上書きされないメモリーと制御記録要件を設計する.

- ビックデータを活用した制御リスク設計,自動車技術会20年秋季大会学術講演会
  - Risk design of control system utilizing big data,Society of Automotive Engineers of Japan
  - https://www.jstage.jst.go.jp/article/jsaeronbun/52/1/52_20214041/_article/-char/ja
  - レーダー認識位置分布(m)と相対速度分布(m/sec)と車線変更制御完了時間(sec)から安全な車線変更制御設計を行う.

- 確率分布ベクトル解析について、情報処理学会第83回全国大会（2021 3/18）
  - https://www.ipsj.or.jp/event/taikai/83/ipsj_web2021/data/pdf/1B-03.html
  - モンテカルロシミュレーションと分布演算、それぞれを使って放物線運動を行う飛距離の分布を求めて比較、正確な分布を求めるために分布演算が必要であることを示す.
  - 6/2にこの論文発表にて優秀賞を受賞しました。以下の83回大会優秀賞に全文が公開されています。
  - https://www.ipsj.or.jp/award/taikaiyusyu.html

## URL
- http://www.na.rim.or.jp/~syn/kakuritsu.html

##ツール

ツールはpython モジュールに加えて、windows10の実行形式を追加した。
windows10の実行形式はGUIによって分布の生成と四則演算など基本的な機能が可能。
pythonモジュールは当面引き上げて、bunpu.exeだけとする。


## windows10実行形式
- bunpu.exe
## How to use bunpu.exe

- インプット分布の選択
  - parameter:パラメータを指定して分布を生成[最小値],[最大値],[平均値],[標準偏差],[分割数]　例　[-1.0,2.0],[0.5,4.0],[-0.1,2.8],[0.2,0.3],[20,20]
  - data_file:ファイルデータから分布を生成、'データファイル名','分布出力名(選択)',無視する行,[取り込む列],[分割数],0,カーネル分布のバンド幅　例　toyota_sony_nikkei.txt,2,[5,11],[20,20],0,20
  - distribution_file:アウトプットした分布の再利用、.npzファイルを選択
  - 上記2番目と3番目はselect_file1,2で選択しても良い
- 演算の選択
  - add:加算
  - sub:減算
  - product:乗算
  - division:除算
  - percent:分布のパーセンタイルを求める
  - percent2:分布と分布のバランスを求める（演習動画２参照）
- execute
  - アウトプットデータ名を指定して演算
  


## ツールの利用方法
## How to use python module


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


課題
- python3.8.5だと3D表示がおかしい（python3.6.9だと正常）→matplotlib2.2.2だと正常、新しいとsavefigがおかしい、サンプルのsavefigを表示だけにする


## ライセンス

© 2020 Shin Koike  bunpu@a1.rim.or.jp

このソフトウェアをそのままの複製を学習や研究を目的として利用する場合、本ソフトウェアおよび今後作成されるものを含めたそのブランチの利用を無償で許可します。

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



=========================================================
Python license


Copyright © 2001-2020 Python Software Foundation; All Rights Reserved



Copyright (c) 2005-2022, NumPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    * Neither the name of the NumPy Developers nor the names of any
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



License agreement for matplotlib versions 1.3.0 and later
=========================================================

1. This LICENSE AGREEMENT is between the Matplotlib Development Team
("MDT"), and the Individual or Organization ("Licensee") accessing and
otherwise using matplotlib software in source or binary form and its
associated documentation.

2. Subject to the terms and conditions of this License Agreement, MDT
hereby grants Licensee a nonexclusive, royalty-free, world-wide license
to reproduce, analyze, test, perform and/or display publicly, prepare
derivative works, distribute, and otherwise use matplotlib
alone or in any derivative version, provided, however, that MDT's
License Agreement and MDT's notice of copyright, i.e., "Copyright (c)
2012- Matplotlib Development Team; All Rights Reserved" are retained in
matplotlib  alone or in any derivative version prepared by
Licensee.

3. In the event Licensee prepares a derivative work that is based on or
incorporates matplotlib or any part thereof, and wants to
make the derivative work available to others as provided herein, then
Licensee hereby agrees to include in any such work a brief summary of
the changes made to matplotlib .

4. MDT is making matplotlib available to Licensee on an "AS
IS" basis.  MDT MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR
IMPLIED.  BY WAY OF EXAMPLE, BUT NOT LIMITATION, MDT MAKES NO AND
DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS
FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF MATPLOTLIB
WILL NOT INFRINGE ANY THIRD PARTY RIGHTS.

5. MDT SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF MATPLOTLIB
 FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR
LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING
MATPLOTLIB , OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF
THE POSSIBILITY THEREOF.

6. This License Agreement will automatically terminate upon a material
breach of its terms and conditions.

7. Nothing in this License Agreement shall be deemed to create any
relationship of agency, partnership, or joint venture between MDT and
Licensee.  This License Agreement does not grant permission to use MDT
trademarks or trade name in a trademark sense to endorse or promote
products or services of Licensee, or any third party.

8. By copying, installing or otherwise using matplotlib ,
Licensee agrees to be bound by the terms and conditions of this License
Agreement.




Pyinstaller License

Can I use PyInstaller for my commercial, closed-source, Python application?
Yes.

If I use PyInstaller for my commercial Python application, will I have to distribute my source code as well?
Absolutely not. You can ship the executables created with PyInstaller with whatever license you want.






