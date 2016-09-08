
Predict-MoRF directory Contains:
1. Demo.m - used to run test.m to predict MoRF scores 
2. test1.m  - used to predict MoRF scores for query sequence residues
3. Model_1to4_Ratio_1_to_2  - trained model 1 to 4 for sample ratio 1:2
4. Model_6to8_Ratio_1_to_2  - trained model 6 to 8 for sample ratio 1:2
5. Model_5and9_Ratio_1_to_1 - trained model 5 and 9 for sample ratio 1:1
6. query_hmm_for_demo - sample hmm profile for demo

Train and test Data and HMM profiles are available at the following web-link:
1. 421 train sequences and 419 test sequences collected by (Disfani et al 2012) (web-link:https://github.com/roneshsharma/Train-and-Test-sequence-data)
2. Test sequences hmm profiles (output from hhblits) (web-link:https://github.com/roneshsharma/Test-hmm)
3. Test sequences hmm profiles filtered and processed, each profile is of size L(sequence length) by 30) (web-link:https://github.com/roneshsharma/Test_hmm_filtered_and-_processed)
4. Train sequences hmm profiles (output from hhblits) (web-link:https://github.com/roneshsharma/Train-hmm-)
5. Train sequences hmm profiles filtered and processed, each profile is of size L(sequence length) by 30) (web-link:https://github.com/roneshsharma/Train-hmm-filtered-and-procesed)


To use test1.m and predict MoRF scores :

1. Install matlab 2015 and above.
2. Download LibSVM (Chang and lin 2011) http://www.csie.ntu.edu.tw/~cjlin/libsvm/
   and build the binaries.
 
Demo example is given to run test1.m:
To use demo.m, LibSVM matlab path must be added in lines 10 and 11



