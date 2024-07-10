# ＜バラツキの対処法＞
# ＜The method for bariance＞

#  バラツキを厳密に扱い、確率を正しく演算するための数学、分布演算を実現するツール


# Mathematical method to handle dispersion strictly and calculate probability correctly, 
#       that is, a tool that realizes distribution calculations.

“バラツキの対処法 ～品質を最大限に引き出す数学～”の中で説明した演算を行うためのソフトウェアです。
このソフトは出版前に公開したものですが、出版後に開発した最新のソフトは以下のアドレスにあります。

This software is used to perform the calculations explained in "Baratsuki no taisyohou (How to deal with variance - Mathematics to maximize quality)"
This software was released before publication, but the latest software developed after publication is at the address below.

https://github.com/skoike/bunpu2
https://github.com/skoike/bunpu3

ソフト比較と使い方/Software comparison and how to use
https://www.youtube.com/watch?v=pnTQOt4HYwg

バラツキの対処法、概要と応用例/How to deal with variance, overview and application examples
https://www.youtube.com/watch?v=YHU92qjDdDA&t=0s  

バラツキの対処法、数学における位置づけ/How to deal with variance, Position in mathematics
https://www.youtube.com/watch?v=bGJhkxFb-2M&t=0s

ソフトの置き場/Software storage
https://github.com/skoike/bunpu
https://github.com/skoike/bunpu2
https://github.com/skoike/bunpu3


　出版前に公開したソフトは、プルダウンメニューで関数を指定するものですが、最新のソフトはユーザーが作成した
プログラムを実行するもので、関数を組合わせることで自由なシミュレーションソフトを作成できます。
初めて使う人は、ここにある出版前のソフトを使ってみるのが良いでしょう。
更に高度な演算やシミュレーションを行う場合は、上記のアドレスにある最新のソフトを使ってください。

The software released before publication allows you to specify functions using a pull-down menu, 
but the latest software executes programs created by users, allowing you to create your own simulation software by combining functions.
If you are a first-time user, you may want to try using the pre-publication software available here.
If you wish to perform more advanced calculations or simulations, please use the latest software available at the address above.

　既存の数値演算ではバラツキを正しく扱うことができないので、見えない大きな誤差が残ります。
標準偏差は、正規分布以外では正しい確率を保証しないし、一般的なバラツキが正規分布に従うことなどありません。
モンテカルロシミュレーションや確率過程も同様に一般的なバラツキに適用すると誤差があります。特にバラツキを持つ要素が複数ある場合には誤差が拡大します。
バラツキの対処法は、そういった既存の数値演算に無い、バラツキの分布形状を演算結果に厳密に反映させる分布演算を扱う環境を提供します。

Existing numerical calculations cannot handle variance correctly, so invisible errors remain.
The standard deviation does not guarantee the correct probability for anything other than a normal distribution, and the general variation does not follow a normal distribution.
Monte Carlo simulations and stochastic processes similarly have errors when applied to general variance. 
Especially when there are multiple elements with variance, the error increases.
This variation handling method provides an environment for handling distribution calculations that strictly reflects the distribution shape of the variation in the calculation results, 
which is not available in existing numerical calculations.

　このソフトウェアはそのままの複製を学習や研究を目的として利用する場合に限り、フリーに使ってもらえます。
それ以外の以下のケースなどは、ライセンス記述にあるアドレスに相談ください。
個別のニーズへの対応は、主に法人を対象として行います。
メールに対する回答は、その要否や期限についてこちらで判断させていただきます。
(3)~(6)はライセンス契約が前提となります。

- 本技術に対するご意見、間違いや改善の提案。(1)
- 本技術を活用する為のコンサルティング、説明、講演などが必要な場合。(2)
- それぞれのニーズにあわせて本技術を応用するためのツールカスタマイズ。(3)
- 本技術を利用したモノやサービスを産業活動（商用）として行う場合(4)
 （検討や試行での利用は自由です）。
- 本技術の関数スクリプト(python)を入手・参考にして自由度の高い活用を行う。(5)
- 本技術を参考にして類似のソフトウェアを開発・配布する場合。(6)

This software may be used free of charge only if you use an exact copy for educational or research purposes.
For other cases such as those listed below, please contact the address listed in the license description.
We respond to individual needs primarily for corporations.
We will decide whether or not a response to the email is necessary and the deadline.
For (3) to (6) a license agreement must be concluded.


- Opinions, mistakes and suggestions for improvement regarding this technology. (1)
- If consulting, explanations, lectures, etc. are required to utilize this technology. (2)
- Tool customization to apply this technology to suit each individual's needs. (3)
- When products and services using this technology are carried out as industrial activities (commercial) (4)
  (You are free to use it for reviewing or trial use.)
- Get or refer to the function script (python) of this technology and use it with a high degree of freedom.(5)
- When developing and distributing similar software using this technology as a reference. (6)

Windows10、Windows11それぞれの環境でコンパイルしたものを用意しました。
まだ、開発途上なので、全ての演算を精度良くカバーできているわけではありません。
現状では起動や演算に時間がかかるので、お待ちください。

I have prepared versions compiled for Windows10 and Windows11 environments.
Since it is still under development, it does not cover all calculations with high accuracy.
Currently, it takes time to start up and calculate, so please wait.

- Windows10用→bunpu_win10.exe
- Windows11用→bunpu_win11.exe

その他csvファイルが多数ありますが、これは上記ツールで練習用に使うダミーデータです。
使い方は“バラツキの対処法”を参照ください。

There are many csv files, this is dummy data used for practice with the above tools.
Please refer to "Baratsuki no taisyohou (How to deal with variance - Mathematics to maximize quality)" for how to use it.


現時点でこのツールは、以下のことが可能である。

- データファイルの特定列から抽出したデータのヒストグラムとカーネル分布を生成
- 範囲や平均値、標準偏差を指定して分布を生成
- 前記生成された分布間の四則演算や時系列積分（シミュレーション）
- 前記生成された分布間の関係を確率的に比較
- 分布のグラフ表示
- 分布要素をファイル出力
- 以上の1次元から3次元のベクトル分布処理（シミュレーションを除く）

ライセンスは、このソフトをそのまま利用するだけであればフリー、ソフトの変更や参考にして作成したものの配布や商用利用する場合は知財権利と著作権にご配慮ください。


At the moment this tool can:

- Generate a distribution from the histogram of the measured data
- Generate distribution by specifying range, mean, and standard deviation
- Four arithmetic operations between generated distributions and time series integral
- Probabilistically compare the relationship between the generated distributions
- Graph display of distribution
- Output distribution elements to a file
- 1D to 3D vector distribution processing above all

The license is free if you just use it as it is, If you change or use this soft as reference for a distribution or commercial purposes,
you should respond to your obligations for the intellectual property rights and the copyright of this software.



## ライセンス

© 2020-2024 Shin Koike  bunpu@a1.rim.or.jp

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


本技術のベースコンセプトは作者がトヨタ自動車在籍中のものですが、この汎用ツール実現の為の全ての発明行為は、トヨタ自動車の知的財産部と自動運転先進安全統括部の間で職務外発明であることの確認エビデンスを取り交わした後に職務外で行ったものです。出版やソフト公開などの行為はトヨタ自動車人事部の確認、指導、了承に基づき行ったものです。


## License

© 2020-2024 Shin Koike  bunpu@a1.rim.or.jp

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

The concept of this technology was developed while I was employed at Toyota Motor Corporation, but all the inventions to realize this general-purpose tool were performed outside of my job. In this regard, Toyota Motor Corporation's Intellectual Property Department and Autonomous Driving Advanced Safety Management Department are exchanging evidence to confirm that this is an off-the-job invention. Actions such as publication and software release were conducted based on the confirmation, guidance, and approval of Toyota Motor Corporation's Human Resources Department.

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






