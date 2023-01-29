# ＜バラツキの対処法＞
# ＜How to deal with bariance＞

#  分布演算によりバラツキを厳密に扱う設計、解析、演算結果の判断方法、捉え方の提案

# Mathematical method for handling dispersion
#  Proposal of new design, new analysis, and understanding using Calculation by distribution

　2023年2月18日技術評論社から発売の“バラツキの対処法”で説明した演算を行うためのソフトウェアです。
書籍ではこのソフトウェアを電卓と呼びましたが、これは一般的な電卓で扱う数字を分布に置換えて、
分布どうしをを演算したり、演算結果を比較するための電卓です。

　既存の数値演算ではバラツキを正しく扱うことができないので、見えない大きな誤差が残ります。
標準偏差は、正規分布以外では正しい領域を保証しないし、そもそもσの位置以外では何も保証されません。
モンテカルロシミュレーションや確率過程も同様です。バラツキの対処法は、そういった既存の数値演算に
無い、バラツキを厳密に処理する数学を提供します。

　このソフトウェアはそのままの複製を学習や研究を目的として利用する場合に限り、フリーに使ってもらえます。
それ以外の以下のケースなどは、ライセンス記述にあるアドレスに相談ください、
個別のニーズへの対応は、5月以降から限られた時間の中で可能な範囲で、主に法人を対象として行います。
メールに対する回答は、その要否や期限についてこちらで判断させていただきます。
(3)~(6)はライセンス契約が前提となります。

- 本技術の間違いや改善の提案。(1)
- 本技術を活用する為のコンサルティング、説明、講演などが必要な場合。(2)
- それぞれのニーズにあわせて本技術を応用するためのツールカスタマイズ。(3)
- 本技術を利用したモノやサービスを産業活動（商用）として行う場合(4)
 （検討や試行での利用は自由です）。
- 本技術の関数スクリプト(python)を入手・参考にして自由度の高い活用を行う。(5)
- 本技術を参考にして類似のソフトウェアを開発・配布する場合。(6)

Windows10で開発していましたが、Windows11ではうまく動かなかったので、
それぞれの環境でコンパイルしたものを用意しました。
まだ、開発途上なので、全ての演算を精度良くカバーできているわけではありません。
現状では起動や演算に時間がかかるので、お待ちください。

- Windows10用→bunpu_win10.exe
- Windows11用→bunpu_win11.exe

その他csvファイルが多数ありますが、これは上記ツールで練習用に使うダミーデータです。
使い方は“バラツキの対処法”を参照ください。



現時点でこのツールは、以下のことが可能である。

- データファイルの特定列から抽出したデータのヒストグラムとカーネル分布を生成
- 範囲や平均値、標準偏差を指定して分布を生成
- 前記生成された分布間の四則演算や時系列積分（シミュレーション）
- 前記生成された分布間の関係を確率的に比較
- 分布のグラフ表示
- 分布要素をファイル出力
- 以上の1次元から3次元のベクトル分布処理（シミュレーションを除く）

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






