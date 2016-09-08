clear all
close all

%Demo script to run test1.m
% loads query hmm profile query_hmm_for_demo.txt
% uses test1.m to predict MoRF scores
% saves residues and scores in a text file scores.txt
% to run demo.m path to LibSVM must be added in lines 10 and 11

%addpath('......./libsvm-3.21/matlab'); %add matlab path for Libsvm 
%addpath('......./libsvm-3.21');

A1= importdata('query_hmm_for_demo.txt');T= struct2cell(A1); % load query hmm profile
seq= T{2,1}(2:end,2); %residues of query sequence
hmmprofile= T{1,1}(:,1:20); %select first 20 columns of hmm profile

MoRFscores = test1(hmmprofile); %run test1

% save residues and scores in text file
fileID = fopen('scores.txt','w');
fprintf(fileID,'%3s  %3s  %6s \n','No:', 'residues','Predict_MoRF_scores' );
for rry=1:size(MoRFscores,2)
fprintf(fileID,'%0.1f  %3s  %f \n',rry,seq{rry,1}(1,1),MoRFscores(rry) );
end
fclose(fileID);
clear rry;clear seq; clear A1;
%##########################################