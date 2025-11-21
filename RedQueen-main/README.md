# RedQueen
MATLAB codes for the paper on vaccination effect on virus evolution (Rouzine 2025):

vaccine_Myrthe.m     % the last vaccination only (Fig. 2A,B,C,F). 

vaccine.m            % frequent vaccination limit with ovelapping effects (Fig 3).  

vaccine_w_Dz.m       % last vaccination only, peak width comparable to the vaccine lag (Fig 2D,E and S6 Text). 

vaccine_exact_cz.m   % last vaccination only, vaccine coordinate right at peak for a finite width of peak w (S5 Text, S7 Figure). 

stabilityanalysis.m  % eigenvalues of red queen problem (S4 Text, S3 Figure). 

scaling_genealogy.m  % testing Eq. S24 in S3 Text with Monte Carlo simulation (S3 Text, S1 Figure).  

recomb_train.m       % Called by scaling_genealogy.m. Monte Carlo evolution program to generate binary genome sequences. 

recomb_train2.m       % Called by testS3andS5.m. A version of recomb_train.m.

main_fitting.m       % Fitting b_y and U_b to data on influenza, no vaccine (S3 Text, S2 Figure). 

call_fitting.m       % Called by main_fitting. Calculates the mean square diff between V and TMRCA and the experimental values. 

testS3andS5.m        % Testing Eq. S26 for the substitution rate in case of uniform distribution of s (S3 Text, S6 figure). 
