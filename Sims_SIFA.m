%% 100 Simulations, data simulated from R.
n=200; %sample size
p=400; % dim of Y, V
q=40;  % dim of X, B
Sim_dataY=cell(1,100);
Sim_dataX=cell(1,100);
for i=1:100
Sim_dataY{i}=DataYAJIVESSn((i-1)*n+1:i*n,:);  
Sim_dataX{i}=DataXAJIVESSn((i-1)*n+1:i*n,:);
end

%% 100 replications, output the estimated matrix B_hat and V_hat into cell (large matrix)
SIFA_V_list=cell(1,100);
SIFA_B_list=cell(1,100);

replication=1;
while replication<=100
   rep_n=replication;
   rep_dataY=Sim_dataY{rep_n};
   rep_dataX=Sim_dataX{rep_n};
   Y_SIFA=cell(1,4);
   for k=1:4
       Y_SIFA{k}=rep_dataY(:,(k-1)*100+1:k*100);
   end
  %SCARF setting:
%   paramstruct.sparsity = 1; 
%   [B0,B,V_joint, V_ind,se2,Sf0,Sf,EU]=SIPCA_A(rep_dataX,Y_SIFA,2,[1,1,0,0],paramstruct);
%   SIFA_B=[B{1},B{2},B0];
%   SIFA_V=[blkdiag(V_ind{1},V_ind{2},V_ind{3},V_ind{4}),V_joint];
%    
  %SIFA / AJIVE setting:
  paramstruct.sparsity = 1; 
  [B0,B,V_joint, V_ind,se2,Sf0,Sf,EU]=SIPCA_A(rep_dataX,Y_SIFA,1,[1,1,1,0],paramstruct);
  SIFA_B=[B0,B{1},B{2},B{3}];
  SIFA_V=[V_joint,blkdiag(V_ind{1},V_ind{2},V_ind{3},V_ind{4})];
    
   SIFA_V_list{rep_n}=SIFA_V;
   SIFA_B_list{rep_n}=SIFA_B;
   
   disp([num2str(replication),' simulations done.']);
   replication=replication+1;
end

%% Save data to csv file, for evaluation in r
for i=1:100
V_SIFA_list((i-1)*p+1:i*p,:)=SIFA_V_list{i};  
B_SIFA_list((i-1)*q+1:i*q,:)=SIFA_B_list{i};
end

csvwrite('SIFA_V_list_AJIVE_S_Sn.csv',V_SIFA_list);
csvwrite('SIFA_B_list_AJIVE_S_Sn.csv',B_SIFA_list);

%% Apply Melanoma Data
paramstruct.sparsity = 0;
SifaY=cell(1,4);
SifaY{1}=MelaSifaY(:,1:18);
SifaY{2}=MelaSifaY(:,19:36);
SifaY{3}=MelaSifaY(:,37:49);
SifaY{4}=MelaSifaY(:,50:62);
[B0,B,V_joint, V_ind,se2,Sf0,Sf,EU]=SIPCA_A(MelaSifaX,SifaY,2,[2,4,2,3],paramstruct);
SIFA_B=[B0,B{1},B{2},B{3}];
SIFA_V=[V_joint,blkdiag(V_ind{1},V_ind{2},V_ind{3},V_ind{4})];
%%
csvwrite('Mela_SIFA_B.csv',SIFA_B);
csvwrite('Mela_SIFA_V.csv',SIFA_V);

%% Apply BRCA Data
paramstruct.sparsity=1;

SifaY_train=cell(1,3);
SifaY_train{1}=DataYtrain(:,1:200);
SifaY_train{2}=DataYtrain(:,201:400);
SifaY_train{3}=DataYtrain(:,401:600);
[B0,B,V_joint, V_ind,se2,Sf0,Sf,EU]=SIPCA_A(DataXtrain,SifaY_train,4,[10,9,3],paramstruct);
SIFA_B=[B0,B{1},B{2},B{3}];
SIFA_V=[V_joint,blkdiag(V_ind{1},V_ind{2},V_ind{3})];

%%
csvwrite('BRCA_SIFA_B_train.csv',SIFA_B);
csvwrite('BRCA_SIFA_V_train.csv',SIFA_V);




