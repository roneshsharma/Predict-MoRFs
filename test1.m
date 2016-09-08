function MoRFscores = test1(hmmprofile)
% MoRFscores = test1(hmmprofile)
%
% MoRF prediction algorithm 
%
% Input
% hmmprofile: HMM profile of size L x 20(where L = length of query protein seq).
%
% Output
% MoRFscore: scores for each residue of the query sequence
%
% Ronesh Sharma, USP, Fiji. 
% Email: sharmaronesh@yahoo.com
% Ref. Sharma et al., Predicting MoRFs in Protein Sequences using HMM Profiles,
% BMC Bioinformatics InCoB2016, 2016

load Model_1to4_Ratio_1_to_2; %load trained model 1 to 4 for sample ratio 1:2
load Model_6to8_Ratio_1_to_2; %load trained model 6 to 8 for sample ratio 1:2
load Model_5and9_Ratio_1_to_1; %load trained model 5 and 9 for sample ratio 1:1
%--------------------------------------------------------------------------
kk=1;Pop_seq(1,1)=1;hl=1; win_flank_siz= 12;mat4=hmmprofile ;T_len=size(mat4,1);
 for e=1:T_len
    mf_1=zeros(((win_flank_siz*2)+1),size(mat4,2)); 
    mf_2=zeros(((win_flank_siz*2)+1),size(mat4,2)); 
    if e<win_flank_siz+1 % seq at starting
      if e>1  
      mf_1((win_flank_siz+2)-e:win_flank_siz,:)=  mat4(1:e-1,:); 
      mf_1((win_flank_siz+1):((win_flank_siz*2)+1),:)=  mat4(e:e+win_flank_siz,:);
      end
      if e==1
      mf_1((win_flank_siz+1):((win_flank_siz*2)+1),:)=  mat4(1:e+win_flank_siz,:) ;
      end 
    sample_d{kk,:}= mf_1; 
    elseif e> T_len-win_flank_siz %seq at ending
        mf_2(1:(win_flank_siz+1),:)= mat4(e-win_flank_siz:e,:);
        if e~=T_len
        mf_2((win_flank_siz+2):(win_flank_siz+1)+(T_len-e),:)= mat4(e+1:end,:);
        end        
    sample_d{kk,:}= mf_2;
    else % seq in middle
    sample_d{kk,:}= mat4(e-win_flank_siz:e+win_flank_siz,:);
    end
    kk=kk+1;
 end
Pop_seq(hl,2) =kk-1;
hl=hl+1;
Pop_seq(hl,1) =kk;
mat4=[];
%--------------------------------------------------------------------------
for ip=1:size(sample_d,1)%feature vectors based on windows size
 flanksreq=1; ALL_hmm3= roundn(sample_d{ip,1}(win_flank_siz-flanksreq+1:win_flank_siz+flanksreq+1,:),-3); %win 3
 flanksreq=2; ALL_hmm5= roundn(sample_d{ip,1}(win_flank_siz-flanksreq+1:win_flank_siz+flanksreq+1,:),-3); %win 5
 flanksreq=3; ALL_hmm7= roundn(sample_d{ip,1}(win_flank_siz-flanksreq+1:win_flank_siz+flanksreq+1,:),-3); %win 7
 flanksreq=4; ALL_hmm9= roundn(sample_d{ip,1}(win_flank_siz-flanksreq+1:win_flank_siz+flanksreq+1,:),-3); %win 9
 flanksreq=5; ALL_hmm11= roundn(sample_d{ip,1}(win_flank_siz-flanksreq+1:win_flank_siz+flanksreq+1,:),-3); %win 11
 flanksreq=6; ALL_hmm13= roundn(sample_d{ip,1}(win_flank_siz-flanksreq+1:win_flank_siz+flanksreq+1,:),-3); %win 13
 feature_pp3(ip,:) =  ALL_hmm3(:)';
 feature_pp5(ip,:) =  ALL_hmm5(:)';
 feature_pp7(ip,:) =  ALL_hmm7(:)';
 feature_pp9(ip,:) =  ALL_hmm9(:)';
 feature_pp11(ip,:)=  ALL_hmm11(:)';
 feature_pp13(ip,:)=  ALL_hmm13(:)';
end
%--------------------------------------------------------------------------
% Predict and combine scores 
[predict_label_L, accuracy_L1, dec_values_1] = svmpredict( ones(size(feature_pp11,1),1), feature_pp11,model_1to4_Ratio_1_to_2.model1,['-q -b 1']);
[predict_label_L, accuracy_L2, dec_values_2] = svmpredict( ones(size(feature_pp7,1),1), feature_pp7,model_1to4_Ratio_1_to_2.model2,['-q -b 1']);
[predict_label_L, accuracy_L3, dec_values_3] = svmpredict( ones(size(feature_pp3,1),1), feature_pp3,model_1to4_Ratio_1_to_2.model3,['-q -b 1']);
[predict_label_L, accuracy_L4, dec_values_4] = svmpredict( ones(size(feature_pp13,1),1), feature_pp13,model_1to4_Ratio_1_to_2.model4,['-q -b 1']);
[predict_label_L, accuracy_L5, dec_values_5] = svmpredict( ones(size(feature_pp9,1),1), feature_pp9,model_5and9_Ratio_1_to_1.model5,['-q -b 1']);
[predict_label_L, accuracy_L6, dec_values_6] = svmpredict( ones(size(feature_pp5,1),1), feature_pp5,model_6to8_Ratio_1_to_2.model6,['-q -b 1']);
[predict_label_L, accuracy_L7, dec_values_7] = svmpredict( ones(size(feature_pp7,1),1), feature_pp7,model_6to8_Ratio_1_to_2.model7,['-q -b 1']);
[predict_label_L, accuracy_L8, dec_values_8] = svmpredict( ones(size(feature_pp13,1),1), feature_pp13,model_6to8_Ratio_1_to_2.model8,['-q -b 1']);
[predict_label_L, accuracy_L9, dec_values_9] = svmpredict( ones(size(feature_pp7,1),1), feature_pp7,model_5and9_Ratio_1_to_1.model9,['-q -b 1']);
MoRFscores = ((dec_values_1(:,1) + dec_values_2(:,1)+ dec_values_3(:,1)+ dec_values_4(:,1)+ dec_values_5(:,1)+ dec_values_6(:,1)+ dec_values_7(:,1)+ dec_values_8(:,1)+ dec_values_9(:,1)) * 1/9)';
clear sample_d;clear feature_pp3;clear feature_pp5; clear feature_pp7;
clear feature_pp9; clear feature_pp11;clear feature_pp13;
end
%##############################