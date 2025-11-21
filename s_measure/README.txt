0.Unpack the archive into a folder on your computer

If you want to check the accuracy of the program, do Step 1.1.

If you want to work with real DNA sequences, do Step 1.2.

1.1 Estimate of selection coefficient from simulated sequence data

1) Open file main.m
2) Set generate=true (by default, generate=false)
3) Set the values of parameters (r,M,L,N,run2,NUbs2,tf,s0,f0,ac,bc,tsec1(2,3),C,appr,apprR)
4) Launch main.m. For example, write main in command line and press "ENTER"
Program main.m will call function recomb_2022.m to generate sequences and function s_measure.m to calculate selection coefficients from them.

1.2 Estimate of selection coefficient from DNA sequence data

1) Download sequences and group them (e.g., using MEGA11) to make a file .fas contain individual DNA from all indepdendent populations at one time period. 
2) Put these three files into folders period1, period2, period3.
3) Open file ACTGtranslate.m and set the parameters defined in the beginning of the file in comments
4) Launch ACTGtranslate.m. The program will generate data file BinData.mat.
5) Open file main.m
6) Set generate = false (it is default)
7) Set parameters (C,aprr,apprR); the other parameters can be ignored for Step 1.2.
8) Launch main.m. The program will call function s_measure.m and make figure similar to Figures  2 and 4 in the text.

After working with a data set, we recommend deleting work space.*ACTGtranslate.m code is upgraded compared to that in the main text for better filtering of deletions and unknown sites.
