# ʬ�۱黻�δĶ����󶡤���
# Providing an environment for distribution calculation


���ӥå��ǡ�������Ϥ�����ǽ���ο��ͱ黻�Ǥ��礭�ʥХ�ĥ���¸�ߤ��븽���������ǥ벽����ˤ����ٸ³������뤳�Ȥ˵��Ť���������к�������ˡ�Ȥ���ʬ�۱黻����Ƥ��Ƥ�����
����Ϻ��ޤ��͡��ʱ��������Ƥ��Ƥ�������������������ǤϤ���ɬ��������ʬ�����򤵤�Ƥ��ʤ��Τǡ�ï�ˤǤ�Ȥ���ġ���Ȥ��ƹ����ȤäƤ�館�뤳�ȤǤ��β��ͤ������¥�ʤ��褦�ȹͤ�����

������ο��ͱ黻�Ǥϡ��Х�ĥ��Τ���ǡ����֤α黻��̤���ʬ�ۤ���Ƥ�Ƹ����Τ���ʬ�ۤ������ʤ�������ϥǡ����������䤷�Ƥ�������ʤ�������ϥ��ƥ�����ߥ�졼�����������ʬ�ۤ��󶡤Ǥ��ʤ����Ȥ�
�Τ��Ƥ��롣������ͳ�ϡ����ͱ黻����Х�ĥ����İ����褦�Ȥ��Ƥ⡢¿���Ǥ���Ω��������ѥ�᡼���֤α黻�Ǥϡ�ʬ�ۤ�������뤿���ɬ�פ��ȹ礻������Ф��Ʒ�»���Ƥ�������̤γ�礬�礭���Τǡ�
�����ˤ�ä�ʬ�ۤη������礭���ƶ���������ǽ�����⤤���ǡ����̤�黻�̤������뤳�Ȥ�ɬ�פ��ȹ礻�����̤��ؿ�Ū�����ä���Τǡ��ǡ��������䤷�Ƥ�������ʬ�ۤ������ʤ���
�к��Ȥ��ơ������Υǡ����Ǥ⡢���줾��Υѥ�᡼���Υҥ��ȥ���फ��ʬ�ۤ�������ơ�ʬ�۱黻��Ԥ����Ȥǡ����ͱ黻�ξ����»����֤�������٤ι⤤�黻����ǽ�Ȥʤ롣

�����󤽤Υġ���ΰ�����������뤬�����ߤϥ٥��ȥ�λ�§�黻��ñ�����ʬ���Ǥ������٤Ǥ��롣�����͡���ʬ������˶��Ϥ���ĺ������ȯ��ʤ�뤳�Ȥ��Ǥ���С�
�͡�����ʬ��������ޤ�������Τ˳�ĥ�Ǥ���ȹͤ��Ƥ��롣���Τ��Ȥ��͡��ʥХ�ĥ��򰷤�ʪ�����ϵ��Ѥ�Ҳ񸽾ݲ��Ϥ˲�����⤿�餷������򿼤�륹�����פȤ��Ƴ��ѤǤ����ǽ�������롣

���Υġ������Ū�ϡ��Х�ĥ�������ǡ�������ҥ��ȥ�����ʬ�ۤ���Ф��ơ�����ʬ�ۤ�黻���ǤȤ���黻�ηϤ��ۤ��뤳�ȤǤ��롣�������Ǥ��Υġ���ϡ��ʲ��Τ��Ȥ���ǽ�Ǥ��롣

### �ǡ����ե�����������󤫤���Ф����ǡ����Υҥ��ȥ�����ʬ�ۤ�����
### �ϰϤ�ʿ���͡�ɸ���к�����ꤷ��ʬ�ۤ�����
### �����������줿ʬ�ۤλ�§�黻����ط����˱�������������ǽ��
### ʬ�ۤΥ����ɽ��
### ʬ�����Ǥ�ե��������
### �ʾ��1��������3�����Υ٥��ȥ�ʬ�۽���

�饤���󥹤ϡ����Υ��եȤ򤽤Τޤ����Ѥ�������Ǥ���Хե꡼�����եȤ��ѹ��仲�ͤˤ��ƺ���������Τ����ۤ侦�����Ѥ�������κ⸢���ˤ�����ͭ���Ȥ��ޤ���


��While analyzing big data, I noticed that there is a limit to the accuracy of modeling the real world with conventional numerical operations,
and as a real example of how to deal with it, I have been published in various papers in Japan.I have been proposing various application examples of this method,
but Few fellow understand its importance, so I promote understanding of its value by making it widely used as a tool that anyone can use. I thought about it.

In the conventional numerical calculation, a reproducible distribution cannot be obtained, if the distribution is obtained
from the calculation results of a plurality of various data, which does not improve even if the number of data is increased.
That is because the Monte Carlo simulation cannot provide the correct distribution, Are known.
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

### Generate a distribution from the histogram of the measured data
### Generate distribution by specifying range, mean, and standard deviation
### Four arithmetic operations between generated distributions
### Graph display of distribution
### Output distribution elements to a file
### 1D to 3D vector distribution processing above all

The license is free if you just use it as it is, If you change or use this soft as reference for a distribution or commercial purposes,
you will be charged for the intellectual property rights.


## ��ʸ
## literature

- �ӥå��ǡ�������Ѥ�������ꥹ���߷�,��ư�ֵ��Ѳ�20ǯ�������ؽѹֱ��
- Risk design of control system utilizing big data,Society of Automotive Engineers of Japan 

## URL

- http://www.na.rim.or.jp/~syn/kakuritsu.html


##�ץ�������ȥ���ݡ��Ȥ���⥸�塼��


- python
- glob
- re
- csv
- os
- pandas
- copy
- scipy
- numpy
- matplotlib
- mpl_toolkits
- sklearn


## �ġ����������ˡ
## How to use


1 from bunpu import *   
2 a=bunpu() 
3 a.bunpu_data('input','output',1,[1,2],0,[20,10],[0],1.5) 
4 b=bunpu() 
5 b.bunpu_gene([-1.0,2.0],[0.5,4.0],[-0.1,2.8],[0.2,0.3],[20,20],'b') 
6 c=bunpu() 
7 c.bunpu_gene([2.0],[4.0],[2.8],[0.3],[20],'c') 
8 d=bunpu() 
9 d=b*c+a 
10 d.bunpu_graph() 


�ƹԤε�ǽ 
1 �ġ���Υ���ݡ��� 
2,4,6,8 ʬ�ۥ᥽�åɤ���� 
3 �ե�����ǡ�������ʬ������������:'�ǡ����ե�����̾','ʬ�۽���̾',̵�뤹���,��������,ʬ���,�������ϰϤ�����,�����ͥ�ʬ�ۤΥХ���� 
5,7 �ϰϤ�ʿ�Ѥ���ꤷ��ʬ�ۤ�����������:�Ǿ���,������,ʿ����,ɸ���к�,ʬ���,ʬ�۽���̾ 
9 ʬ�۱黻���Ѿ���¿����*1�������º���¿����*¿������ 
10 �����ɽ�� 
    


## �饤����
## License
��� 2020 Shin Koike  bunpu@a1.rim.or.jp

���Υ��եȥ������򤽤Τޤޤ�ʣ���Ȥ������Ѥ����硢�ܥ��եȥ���������Ӻ������������Τ�ޤ᤿���Υ֥��������Ѥ�̵���ǵ��Ĥ��ޤ���

���Υ��եȥ�������̤�����ǡ���������Ƥ䵡ǽ��ĥ�ζ��Ϥ���Ƥ��ޤ������Υ��եȤβ����䶨�Ϥΰ٤ˤˡ��ѹ����ɲá���硢�ܿ���ޤ�������
���Ѳ�ǽ�ʾ���ȤȤ�ˡ�����������Ȥ��ơ����Υ��եȤκ�Ԥޤ�������Ԥˤ��ξ����󶡤򤪴ꤤ���ޤ���
�������Ƥϸ������˴�Ť����ܥ��եȤޤ��Ϥ��Υ֥�����ȿ�Ǥ����Ƥ����ޤ���

���Υ��եȤ����ѡ����ͤˤ�����ϡ����Υ��եȤ�������õ��д��PCT/JP2020/034566�Ȥ���ʬ�䡢��Ϣ�д�ˤ���Ӥ��ζ��ϼԤˤ����븢����º�Ť���������
���Υ��եȥ������ΰ���ʬ�����Ѥޤ��ϻ��ͤˤ��ơ��ѹ����ɲá���硢�Ѿ���ܿ���ޤ����������ۤޤ��Ͼ������Ѥ�����Ϻ�Ԥ����̤��Ƥ���������

���եȥ������ϡ�̤�����ǡ�������ݾڤ�ʤ��󶡤���ޤ���
�����Ǥ����ݾڤȤϡ����������������Ū�ؤ�Ŭ����������Ӹ����󿯳��ˤĤ��Ƥ��ݾڤ�ޤߤޤ���������˸��ꤵ����ΤǤϤ���ޤ��� 
���Υ��եȺ�Ԥޤ�������Ԥϡ�����԰١���ˡ�԰١��ޤ��Ϥ���ʳ��Ǥ����ȡ����եȥ������˵����ޤ��ϴ�Ϣ�������뤤�ϥ��եȥ������λ��Ѥޤ���
����¾�ΰ����ˤ�ä���������ڤ����ᡢ»��������¾�ε�̳�ˤĤ��Ʋ������Ǥ�����ʤ���ΤȤ��ޤ���

�ʾ��ɽ��������ܵ���ɽ���򡢥��եȥ������Τ��٤Ƥ�ʣ���ޤ�����ʬ�����Ѥޤ���ʬ�۽����򻲹ͤȤ�����ˡ��������������ʪ�˵��ܤ����ΤȤ��ޤ���

��� 2020 Shin Koike  bunpu@a1.rim.or.jp

Permission is hereby granted, free of charge, to any person obtaining a exact copy of this software,
its branches and associated documentation files (the "Software"), to deal in the Software with restriction.

This software is incomplete and we are seeking suggestions for improvement and cooperation in enhancements.
For the improvement and cooperation of this software, please provide the derivation
including modification, addition, mergers,combination, translation with available information to the author
or copyright holder of this software on the assumption that it will be published.
The contents will be reflected in this software and its branches based on public nature and my leeway.

When using or referring to this software, please correspond the copyright of this software
and the rights in patent applications(PCT/JP2020/034566 and divisional other).
Please contact with the author if you want to use or refer to a part of this software and distribute it privately or use it for commercial purposes.

The software is incomplete and is provided without warranty.Warranties here include, but are not limited to, warranties of merchantability, 
fitness for a particular purpose, and non-infringement.
The author or copyright holder of this software, whether contractual, tort, or otherwise, is due to or related to the software, or uses or uses the software.
We shall not be liable for any claims, damages or other obligations arising from any other dealings.

The above copyright notice and this permission notice shall be included in all copies, portions or reference of the software related to dstribution .


