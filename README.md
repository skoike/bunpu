
#"���z���Z�̊���񋟂���"

#Providing an environment for distribution calculation


�@�r�b�N�f�[�^����͂��钆�ŏ]���̐��l���Z�ł͑傫�ȃo���c�L�����݂��錻�����E�����f��������ɂ͐��x���E�����邱�ƂɋC�Â��A�����΍􂷂���@�̋�̗��l�X�Ș_���ɂĔ��\���Ă����B
���̕��@��25�N�ȏ�O����l�X�ȉ��p����Ă��Ă������A���܂��ɂ��̕K�v������������Ă��Ȃ��̂ŁA����A��Ђ̋���������ĒN�ɂł��g����c�[���Ƃ��čL���g���Ă��炦�邱�Ƃł��̉��l�̗����𑣐i���悤�ƍl�����B

�@�]���̐��l���Z�́A�o���c�L�̂��镡���f�[�^�̉��Z���ʂ��番�z�����߂Ă��Č����̂��镪�z�������Ȃ��A����̓f�[�^���𑝂₵�Ă����P���Ȃ��B����̓����e�J�����V�~�����[�V���������������z��񋟂ł��Ȃ����Ƃ�
�m���Ă���B���̗��R�́A���l���Z����o���c�L��c�����悤�Ƃ��Ă��A�����ł��Ɨ���������p�����[�^�Ԃ̉��Z�ł͏�񌇑����傫���̂ŁA�񌻎��I�Ȑ��̏璷�ȃf�[�^�������Ƃ܂Ƃ��Ȑ��x�̕��z�������Ȃ�����ł���B
�΍�Ƃ��āA�����̃f�[�^�ł��A���ꂼ��̃p�����[�^�̃q�X�g�O�������番�z���쐬���āA���z���Z���s�����ƂŁA���l���Z�̏�񌇑�����Ԃ���Đ��x�̍������Z���\�ƂȂ�B
�@���񂻂̃c�[���̈ꕔ�����J���邪�A���݂̓x�N�g���̎l�����Z��P���Ȑϕ����ł�����x�ł���B����l�X�ȕ���̕��ɋ��͂��Ē����A�J����i�߂邱�Ƃ��ł���΁A�l�X�Ȕ������������܂ސ��w�S�̂Ɋg���ł���ƍl���Ă���B
���̂��Ƃ̓i�r�G�X�g�[�N�X���͂��߂Ƃ���A�l�X�Ȋ�{�I�ȕ�����͋Z�p�ɉ��P�������炵�A������[�߂�X�R�[�v�Ƃ��Ċ��p�ł���\��������B

���̃c�[���̖ړI�́A�o���c�L������f�[�^����q�X�g�O�����ŕ��z�𒊏o���āA���̕��z�����Z�v�f�Ƃ��鉉�Z�̌n���\�z���邱�Ƃł���B�����_�ł��̃c�[���́A�ȉ��̂��Ƃ��\�ł���B
���C�Z���X�́A���̂܂ܗ��p���邾���ł���΃t���[�A�ύX����ꍇ�͌��J�O��Œ񋟒�����΃t���[�A�ύX��Q�l�Ƃ������̂̔z�z�⏤�p���p����ꍇ�͒m�������ɂ����ėL���Ƃ��܂��B

###�t�@�C���̓���񂩂璊�o�����f�[�^�̃q�X�g�O�����ƕ��z�𐶐�
###�͈͂╽�ϒl�A�W���΍����w�肵�ĕ��z�𐶐�
###�l�����Z�i���֌W���ɉ������␳���\�j
###���z�̃O���t�\��
###���z�v�f���t�@�C���o��
###�ȏ��1��������3�����̃x�N�g�����z����

�@While analyzing big data, I noticed that there is a limit to the accuracy of modeling the real world with conventional numerical operations, and as a real example of how to deal with it, I have been published in various papers in Japan.I have been proposing various application examples of this method for more than 20 years, but Few fellow understand its importance, so I promote understanding of its value by making it widely used as a tool that anyone can use this time with the permisson of my company. I thought about it.

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


##�_��
##literature

�r�b�N�f�[�^�����p�������䃊�X�N�݌v,�����ԋZ�p��20�N�H�G���w�p�u����
Risk design of control system utilizing big data,Society of Automotive Engineers of Japan 

##�v���O��������ƃC���|�[�g���郂�W���[��
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

##�c�[���̗��p���@
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

    �e�s�̋@�\
    1 �c�[���̃C���|�[�g
    2,4,6,8 ���z���\�b�h�̒�`
    3 �t�@�C���f�[�^���番�z�����A����:'�f�[�^�t�@�C����','���z�o�͖�',��������s,��荞�ޗ�,������,��荞�ޔ͈͂̐���,�J�[�l�����z�̃o���h��
    5,7 �͈͂╽�ς��w�肵�ĕ��z�𐶐��A����:�ŏ��l,�ő�l,���ϒl,�W���΍�,������,���z�o�͖�
    9 ���z���Z�i�Ϗ��͑�����*1�����A�a���͑�����*�������j
    10 �O���t�\��
    


##���C�Z���X
##License

���̃\�t�g�E�F�A�����̂܂܂̕����Ƃ��ė��p����ꍇ�A�{�\�t�g�E�F�A�̗��p�𖳏��ŋ����܂��B

���̃\�t�g�E�F�A�͖������ŁA���P�̒�Ă�@�\�g���̋��͂����߂Ă��܂��A���̃\�t�g�̉��P�⋦�ׂ͂̈ɂɁA�ύX�A�ǉ��A�����A�ڐA���܂ޔh�����A
���p�\�ȏ��ƂƂ��ɁA���J��O��Ƃ��āA���̃\�t�g�̍�҂܂��͒��쌠�҂ɂ��̏��񋟂��s�����ƂŁA���̗��p�������܂��B

���̃\�t�g�𗘗p�E�Q�l�ɂ���ꍇ�́A���̃\�t�g�̒��쌠�Ɠ����o��iPCT/JP2020/034566�Ƃ��̕����A�֘A�o��j�ɂ����錠���𑸏d���������B
���̃\�t�g�E�F�A�̈ꕔ���𗘗p�܂��͎Q�l�ɂ��āA�ύX�A�ǉ��A�����A�p����ڐA���܂ޔh�����A����J�Ŕz�z�A�܂��͏��p���p����ꍇ�͍�҂ɑ��k���Ă��������B

�\�t�g�E�F�A�́A�������ŁA����̕ۏ؂��Ȃ��񋟂���܂��B
�����ł����ۏ؂Ƃ́A���i���A����̖ړI�ւ̓K�����A����ь�����N�Q�ɂ��Ă̕ۏ؂��܂݂܂����A����Ɍ��肳�����̂ł͂���܂���B 
���̃\�t�g��҂܂��͒��쌠�҂́A�_��s�ׁA�s�@�s�ׁA�܂��͂���ȊO�ł��낤�ƁA�\�t�g�E�F�A�ɋN���܂��͊֘A���A���邢�̓\�t�g�E�F�A�̎g�p�܂���
���̑��̈����ɂ���Đ������؂̐����A���Q�A���̑��̋`���ɂ��ĉ���̐ӔC������Ȃ����̂Ƃ��܂��B

�ȏ�̕\������і{�����\�����A�\�t�g�E�F�A�̂��ׂĂ̕����܂��͕����̗��p�܂��͕��z�������Q�l�Ƃ���ꍇ�ɁA�쐬����钘�앨�ɋL�ڂ�����̂Ƃ��܂��B


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


