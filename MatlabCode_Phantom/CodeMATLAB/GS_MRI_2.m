clear;
%% Gold standard T1 and T2 measurements
 
folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/GS/';
 
% % Gold standard T2 measurements 
MESE(:,:,1)= double(dicomread([folder 'IM_0107']));
MESE(:,:,2)= double(dicomread([folder 'IM_0108']));
MESE(:,:,3)= double(dicomread([folder 'IM_0109']));
MESE(:,:,4)= double(dicomread([folder 'IM_0110']));
MESE(:,:,5)= double(dicomread([folder 'IM_0111']));
MESE(:,:,6)= double(dicomread([folder 'IM_0112']));
MESE(:,:,7)= double(dicomread([folder 'IM_0113']));
MESE(:,:,8)= double(dicomread([folder 'IM_0114']));
%%
% Floating Point DICOM 
% Gold standard T1 measurements 
T1_1= double(dicomread([folder 'IM_0116']));
infoT1_1 =  dicominfo([folder 'IM_0116']);
T1_1_FP = T1_1/infoT1_1.Private_2005_100e  -  infoT1_1.Private_2005_100d/infoT1_1.Private_2005_100e;

T1_2= double(dicomread([folder 'IM_0119']));
infoT1_2 =  dicominfo([folder 'IM_0119']);
T1_2_FP = T1_2/infoT1_2.Private_2005_100e  -  infoT1_2.Private_2005_100d/infoT1_2.Private_2005_100e;

T1_3= double(dicomread([folder 'IM_0122']));
infoT1_3 =  dicominfo([folder 'IM_0122']);
T1_3_FP = T1_3/infoT1_3.Private_2005_100e  -  infoT1_3.Private_2005_100d/infoT1_3.Private_2005_100e;

T1_4= double(dicomread([folder 'IM_0125']));
infoT1_4 =  dicominfo([folder 'IM_0125']);
T1_4_FP = T1_4/infoT1_4.Private_2005_100e  -  infoT1_4.Private_2005_100d/infoT1_4.Private_2005_100e;

T1_5= double(dicomread([folder 'IM_0128']));
infoT1_5 =  dicominfo([folder 'IM_0128']);
T1_5_FP = T1_5/infoT1_5.Private_2005_100e  -  infoT1_5.Private_2005_100d/infoT1_5.Private_2005_100e;

% T1_6= double(dicomread([folder 'IM_0160']));
% infoT1_6 =  dicominfo([folder 'IM_0160']);
% T1_6_FP = T1_6/infoT1_6.Private_2005_100e  -  infoT1_6.Private_2005_100d/infoT1_6.Private_2005_100e;
% % T1_6_FP = T1_6_FP./max(T1_6_FP(:));
% 
% T1_7= double(dicomread([folder 'IM_0163']));
% infoT1_7 =  dicominfo([folder 'IM_0163']);
% T1_7_FP = T1_7/infoT1_7.Private_2005_100e  -  infoT1_7.Private_2005_100d/infoT1_7.Private_2005_100e;
% % T1_7_FP = T1_7_FP./max(T1_7_FP(:));
% 
% T1_8= double(dicomread([folder 'IM_0166']));
% infoT1_8 =  dicominfo([folder 'IM_0166']);
% T1_8_FP = T1_8/infoT1_8.Private_2005_100e  -  infoT1_8.Private_2005_100d/infoT1_8.Private_2005_100e;
% % T1_8_FP = T1_8_FP./max(T1_8_FP(:));
% 
% T1_9= double(dicomread([folder 'IM_0169']));
% infoT1_9 =  dicominfo([folder 'IM_0169']);
% T1_9_FP = T1_9/infoT1_9.Private_2005_100e  -  infoT1_9.Private_2005_100d/infoT1_9.Private_2005_100e;
% % T1_9_FP = T1_9_FP./max(T1_9_FP(:));
% 
% T1_10= double(dicomread([folder 'IM_0172']));
% infoT1_10 =  dicominfo([folder 'IM_0172']);
% T1_10_FP = T1_10/infoT1_10.Private_2005_100e  -  infoT1_10.Private_2005_100d/infoT1_10.Private_2005_100e;
% % T1_10_FP = T1_10_FP./max(T1_10_FP(:));

VTI(:,:,1)= T1_1_FP;
VTI(:,:,2)= T1_2_FP;
VTI(:,:,3)= T1_3_FP;
VTI(:,:,4)= T1_4_FP;
VTI(:,:,5)= T1_5_FP; 
% VTI(:,:,6)= T1_6_FP;
% VTI(:,:,7)= T1_7_FP;
% VTI(:,:,8)= T1_8_FP;
% VTI(:,:,9)= T1_9_FP;
% VTI(:,:,10)= T1_10_FP;

TI = [100, 400, 800, 1500, 3000]';

% TI = [35, 75, 100, 125,150, 250, 1000, 1500, 2000, 3000]';
TE = 22:22:22*8;
TE = TE';
N = size(VTI,1);

VTI = double(permute(reshape(VTI,N^2,size(VTI,3)),[2,1]));
MESE = permute(reshape(MESE,N^2,size(MESE,3)),[2,1]);

% Define the exponential model for T1
% model = @(params, TI)abs(params(1)*((1-(1 + 1)*exp(-TI/params(2))) + exp(-4500/params(2))));
model = @(params, TI)abs(params(1)*((1 - (1 + params(2))*exp(-TI/params(3))) + params(2)*exp(-4500/params(3))));

% Define the parameter bounds and starting points
lb = [1, 0, 1];  
ub = [Inf, 1, 4000];
params0 = [3000, 1, 600];

T1 = zeros(N^2,1);
invEff = zeros(N^2,1);
ms = MultiStart('UseParallel',true);

%Set up the problem for MultiStart.
for m = 1:N^2
    problem = createOptimProblem('lsqcurvefit','x0',params0,'objective',model,...
        'lb',lb,'ub',ub,'xdata',TI,'ydata',VTI(:,m)); %,'options',opts
    [paramMulti,~] = run(ms,problem,50);
    T1(m) = paramMulti(3);
   % invEff(m) = paramMulti(2); 
end 


% Define the exponential model for T2
modelT2 = @(paramsT2, TE)(paramsT2(1)*exp(-TE/paramsT2(2)));

% Define the parameter bounds and starting points
lbT2 = [1, 1];  
ubT2 = [Inf, 4000];
params0T2 = [1000, 500];

T2= zeros(N^2,1);
%Set up the problem for MultiStart.
for m = 1:N^2
    problem = createOptimProblem('lsqcurvefit','x0',params0T2,'objective',modelT2,...
        'lb',lbT2,'ub',ubT2,'xdata',TE,'ydata',MESE(:,m)); %,'options',opts
    [paramMultiT2,~] = run(ms,problem,50);
    T2(m) = paramMultiT2(2); 
end 

%%
addpath('C:\Users\lucya\MSC_PROJECT\MatlabCode_Phantom\CodeMATLAB\utils');


T1 = reshape(T1, N, N);
% invEff = reshape(invEff, N, N);
T2 = reshape(T2, N,N);
mask = double(dicomread([folder 'IM_0107']));
mask(mask < 60) = 0;
mask(mask ~= 0) = 1;

T1map= bsxfun(@times,T1,mask);
T1map(isnan(T1map))=0;
figure; imshow3(T1map, [0 3000]);
figure; imshow3(T1map, [0 3000]);
colormap(inferno);

% 
% invEffmap= bsxfun(@times,invEff,mask);
% invEffmap(isnan(invEffmap))=0;
% figure; imshow3(invEffmap, [0 1]);

T2map= bsxfun(@times,T2,mask);
T2map(isnan(T2map))=0;
figure; imshow3(T2map, [0 1500]);
figure; imshow3(T2map, [0 1500]);
colormap(turbo);




% Masked maps already computed
save('T1T2_results.mat', 'T1', 'T2', 'T1map', 'T2map', 'mask');


