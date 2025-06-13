%% Magnetic Resonance Fingerprinting reconstruction code
% SUMMARY: Code shows one example of an in-vivo abdominal MRF
% reconstruction on a 1.5T MR-linac. Reconstruction includes, dictionary
% generation and mrf reconstruction.
%
% Notes:
%  - Make sure you have write permission in the working directory
%  - You must have BART version 0.4.04 or higher, otherwhise no support for
%       pics temporal base.
%  - The code for the GIRF corrections is not included in the script yet.
%  - The dictionary generation code is slow, on my own machine I use Julia based code, which I am not able to provide. 
% 
% Dependencies: 
%  - I use some code from Jakob's Asslaender's repository on low-rank ADMM
%  MRF: https://bitbucket.org/asslaender/nyu_mrf_recon/src/master/
%  - The BART toolbox from berkeley: https://bitbucket.org/asslaender/nyu_mrf_recon/src/master/
%  - Extended phase graph code from M. Weigel: http://epg.matthias-weigel.net/
%
% If you find this code useful please cite:
%   Bruijnen, T., van der Heide, O., Intven, M. P. W., Mook, S., Lagendijk, J. J. W., van den Berg, C. A. T., & Tijssen, R. H. N. (2020).
%   Technical feasibility of Magnetic Resonance Fingerprinting on a 1.5T MRI-Linac.
%   Phys Med Biol. 2020 Sep 25. doi: 10.1088/1361-6560/abbb9d. Epub ahead of print. PMID: 32977318.
%
% Contact: T.Bruijnen@umcutrecht.nl | University Medical Center Utrecht
% Department of Radiotherapy, University Medical Center Utrecht, Utrecht, the Netherlands
% Computational Imaging Group for MRI diagnostics and therapy, Centre for
% Image Sciences, University Medical Center Utrecht, Utrecht, The Netherlands
%
% Init: T.Bruijnen - 20200101
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Initialize paths, check for bart toolbox and load data

% Required: BART (installed in WSL)
%bart_loc = '/home/mnuixe/BART/bart-0.8.00/'; % Make sure this ends with a forward slash
%bart_loc = '/home/awetscherek/bart-0.8.00/'; % Make sure this ends with a forward slash
bart_loc = '/mnt/c/Users/lucya/MSC_PROJECT/bart/';


% matlab folder of BART contains a script to call BART in WSL
%addpath(genpath('D:\Codes\MATLAB\bart-0.8.00\matlab\'))
%addpath(genpath('G:/My Drive/MATLAB/bart-master/matlab/'))
addpath(genpath('C:/Users/lucya/MSC_PROJECT/bart/matlab/'));
setenv('TOOLBOX_PATH',bart_loc);


% check if bart toolbox is installed
if ~isempty(which('bart')); else;disp('>> BART toolbox not found in path. Please check BART path.');end

% add the "utils" subfolder to the PATH:
%addpath(genpath('C:/Codes/MATLAB/MRF/mrf-mrl-master/'))
%addpath(genpath('C:/data/MRFData_Lucy/CodeMATLAB/utils'))
addpath(genpath('C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/CodeMATLAB/utils'));

% folder that contains the raw data:
%folder = 'D:\1_Data\2025\Article1\Phantom\Data\Session1\';
%folder = 'C:/data/MRFData_Lucy/';
folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/';

%% Load data  
data = loadRawKspacePhilips([folder 'raw_001.list']);

% raw data and image dimensions (according to .list file):
Nkx = abs(data.kspace_properties.kx_range(2) - data.kspace_properties.kx_range(1)) + 1;
Nky = abs(data.kspace_properties.ky_range(2) - data.kspace_properties.ky_range(1)) + 1;
% Nkz = abs(data.kspace_properties.kz_range(2) - data.kspace_properties.kz_range(1)) + 1;
Nx = abs(data.kspace_properties.X_range(2) - data.kspace_properties.X_range(1)) + 1;
Ny = abs(data.kspace_properties.Y_range(2) - data.kspace_properties.Y_range(1)) + 1;
% Nz = abs(data.kspace_properties.Z_range(2) - data.kspace_properties.Z_range(1)) + 1;
Nc  = numel(unique(data.chan));
% dyn=numel(unique(data.dyn));
dyn=1;
neco = data.kspace_properties.number_of_echoes;
slices=numel(unique(data.loca));

%prompt = "What is the distance of the slice from the isocentre (in mm)? ";
%distanceSlice = input(prompt); % I think for this data set it should be 0 here ...
%distanceSlice = distanceSlice/1000;
distanceSlice = 0;

% extract noise adjustment data:
noise = complex(zeros(numel(data.complexdata{find(cell2mat(data.typ) == 'NOI', 1)}), Nc));
for ii = 1:numel(data.complexdata)
    if strcmp(data.typ{ii}, 'NOI')
        noise(:, data.chan(ii) + 1) = data.complexdata{ii};
    end
end

% put raw data in more useful structure. ky index is problematic, we use the knowledge about the
% order in which k-space data is acquired:
raw = complex(zeros(Nkx,Nc,Nky,dyn));
for ii = 1:(numel(data.complexdata))
    if strcmp(data.typ{ii}, 'STD')
               raw(:,  data.chan(ii)+1, data.ky(ii) +1 ,data.dyn(ii)+1) = data.complexdata{ii};
    end
end

%% Data preprocessing
U.RawKspaceData=squeeze(raw(:,:,:,1:dyn));
U.RawKspaceData= permute(U.RawKspaceData,[1 3 4 2]);
U.RawKspaceData = noise_prewhitening(U.RawKspaceData,noise);
kdim            = c12d(size(U.RawKspaceData));

%prompt = "What is the slice orientation (1- coronal, 2- transverse, 3- sagittal)? ";
%sliceOrientation = input(prompt); % for the example data here, this should be coronal (based on the DICOM attribute Image Position Patient 1/0/0/0/0/-1)
sliceOrientation = 1;

%if sliceOrientation ~= 1 && sliceOrientation ~= 2 && sliceOrientation ~= 3 
   %disp("Wrong entry");
%end 

kx_GIRF = fread(fopen([folder 'kx_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 
ky_GIRF = fread(fopen([folder 'ky_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 
kz_GIRF = fread(fopen([folder 'kz_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 
b0_GIRF = fread(fopen([folder 'b0_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 

b0_GIRF= b0_GIRF.*(2* pi* 42.57746778); %scale the error by multiplying by 2* pi* gyromagnetic ratio 42.57746778 MHz/T same as Philips 
trajx_GIRF = reshape(kx_GIRF,[Nkx,Nky*dyn]);
trajy_GIRF = reshape(ky_GIRF,[Nkx,Nky*dyn]);
trajz_GIRF = reshape(kz_GIRF,[Nkx,Nky*dyn]);
b0_GIRF = reshape(b0_GIRF,[Nkx,Nky*dyn]);

if sliceOrientation == 1 
    trajx_GIRF = trajx_GIRF * 2 * pi * distanceSlice; % put the output of the GIRF in rad (for coronal slice)
    trajGIRF(1,:,:)= trajz_GIRF; 
    trajGIRF(2,:,:) = trajy_GIRF;
    trajGIRF(3,:,:) = 0;
elseif sliceOrientation == 2
    trajz_GIRF = trajz_GIRF * 2 * pi * distanceSlice; % put the output of the GIRF in rad (for transverse slice)
    trajGIRF(1,:,:)= trajx_GIRF;  
    trajGIRF(2,:,:) = trajy_GIRF;
    trajGIRF(3,:,:) = 0;
elseif sliceOrientation == 3
    trajy_GIRF = trajy_GIRF * 2 * pi * distanceSlice; % put the output of the GIRF in rad (for sagittal slice)
    trajGIRF(1,:,:)= trajx_GIRF; 
    trajGIRF(2,:,:) = trajz_GIRF;
    trajGIRF(3,:,:) = 0;
else 
    disp('An error occured')
end

trajGIRF        = repmat(trajGIRF,[1 1 1 1 kdim(6)]);
dcf             = permute(radial_density(trajGIRF(:,:,:,:,1)),[1 2 3 5 4]);
% N               = ceil(max(abs(trajGIRF(:)))/2);
N = 124;
comb_raw=zeros(kdim(1),Nc,kdim(2)*dyn);
location=0;
step=1000;

U.RawKspaceData= permute(U.RawKspaceData,[1 4 2 3]); 
for d=1:dyn
    comb_raw(:,:,location+1:step)= squeeze(U.RawKspaceData(:,:,:,d));
    location=location+1000;
    step=step+1000;
end

U.RawKspaceData= permute(comb_raw,[1 3 4 2]);

% Add b0 eddy current phase correction of the raw data 
U.RawKspaceData = U.RawKspaceData.*exp(-1i*(b0_GIRF)); %  

% Also need after to apply the z-GIRF phase 
if sliceOrientation == 1 
    U.RawKspaceData = U.RawKspaceData.*exp(-1i*(trajx_GIRF));
elseif sliceOrientation == 2
    U.RawKspaceData = U.RawKspaceData.*exp(-1i*(trajz_GIRF));
elseif sliceOrientation == 3
    U.RawKspaceData = U.RawKspaceData.*exp(-1i*(trajy_GIRF));
else 
    disp('An error occured')
end

prompt = "Is there a B1 correction (1- Yes, 0- No)? ";
B1correct = input(prompt);

% Reconstruct multi-channel images to estimate CSM
nufft = bart(['nufft -a -d ',num2str(N),':',num2str(N),':1'],0.25*trajGIRF,permute(bsxfun(@times,dcf,U.RawKspaceData),[3 1 2 4])); %Phantom 0.25
imgForCSM = bart('fft 3',nufft);
csm   = bart('ecalib -m1',imgForCSM); %put a threshold ? 

%% Dictionary 

% loading dictionary from MATLAB file (5865 entries):
%load D_Phantom2025_invEff096_SPinf_norefocusingTEadj_576InvTime_1000RF_10mm_101iso_0.mat
%load D_IDX_SP_Phantom2025.mat % apparently this is also needed - these are likely the T1, T2 (and B1 (?)) values for each dictionary entry

folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/';
load(fullfile(folder, 'D_Phantom2025_invEff096_SPinf_norefocusingTEadj_576InvTime_1000RF_10mm_101iso_0.mat'));
load(fullfile(folder, 'D_IDX_SP_Phantom2025.mat'));


% alternatively: the dictionary with B1:
%load Dictionary_B1_SP.mat
load(fullfile(folder, 'Dictionary_B1_SP.mat'));


dict0 = dict0(1:1000*dyn,:);
[~,s,~]=svd((dict0),'econ');
result_s    = diag(s);
NRJ= zeros(100,1);
for R = 1:100
    NRJ(R,1) = sum(result_s(1:R))/sum(result_s);
end

%R = find(NRJ>=0.999,1);
R = 20;  % Try 10–30 to start small

% dict must be under the form timeFrame * parameters combinations
[D] = compression_mrf_dictionary(dict0,idx, R);
%%D = MRF_dictionary_umcu(idx', R, double(sum(dict,3))); % epg 

% First write all the cfl/hdr files in the correct dimensions for BART
writecfl('csm',single(csm));
writecfl('u',permute(D.u,[3 4 5 6 7 1 2]));
writecfl('traj',permute(.25 * trajGIRF,[1 2 5 4 6 3]))
writecfl('dcf',permute(sqrt(dcf),[5 1 3 4 6 2]));
writecfl('kdata',permute(U.RawKspaceData,[3 1 5 4 6 2]));

% LR inversion parameters
lambda = 0.05;
n_iter = 50;

% Reconstruct singular value images with BART
tic 
% for me this line did not work, I had to add [], but I might be using an
% old wrapper script ...
%recon = bart('bart pics -d1 -e -l1 -r ',num2str(lambda),' -i',num2str(n_iter),' -p ','dcf',' ',' -t ','traj',' ','-B ','u',' ','kdata',' ','csm');
cmd = sprintf('pics -d1 -e -l1 -r %f -i %d -p dcf -t traj -B u kdata csm', lambda, n_iter);
recon = bart(cmd);


writecfl('recon',single(recon));
toc 

svd_images= reshape(readcfl('recon'),N,N,R);
match_images = reshape(svd_images,N^2,R);

c   = zeros([size(match_images,1) 1],'single');
idx2 = zeros([size(match_images,1) 1],'single');

if B1correct ~= 1 && B1correct ~= 0
   disp("Wrong entry");
elseif B1correct == 1
    B1   = [0.8:0.02:0.9, 0.91:0.01:1.09, 1.10:0.02:1.2] ;
    % need to specify the folder for the B1 map...
    %folderB1 = 'D:\1_Data\2025\Article1\Phantom\Data\Session1\DICOM\';
    folderB1 = 'C:/data/MRFData_Lucy/DICOM/';
    mapB1 = double(dicomread([folderB1 'IM_0101']));
    info_mapB1 = dicominfo([folderB1 'IM_0101']);
    
    mapB1_FP = mapB1/info_mapB1.Private_2005_100e  -  info_mapB1.Private_2005_100d/info_mapB1.Private_2005_100e;
    mapB1_percent = mapB1_FP/100; % Map on which I want to apply the transformation
    
    % Interpolate to get the same size as the MRF images 
    sz_for_registration = size(mapB1_percent);
    x_before = 1:sz_for_registration(1);
    y_before = 1:sz_for_registration(2);
    
    F_mapB1percent = griddedInterpolant({x_before,y_before},double(mapB1_percent),'nearest');
    
    x_after = (0:2.0726:256)'; %From 256 to 124  %From 128 to 124:  1.03876
    y_after = (0:2.0726:256)';
   
    % x_after = (0:1.0827:432)'; %From 256 to 124  %From 128 to 124:  1.03876
    % y_after = (0:1.0827:432)';

    resized_B1_map_percent = (F_mapB1percent({x_after,y_after}));
    resized_B1_map_percent = fliplr(flip(resized_B1_map_percent));
    resized_B1_map_percent = fliplr(resized_B1_map_percent); %% yes for phantom
    
    %Dictionary match on singular value images
    
    resized_B1_map_percent = reshape(resized_B1_map_percent, N^2,1);
    
    B1D = double(D.lookup_table');
    nbB1 = size(B1D,2)/size(B1,2);
    
    diff = single(zeros(size(B1D(3,:))));
    indices = single(zeros(1,nbB1));
    
    for n = 1 : N^2
        diff = abs(B1D(3,:)-resized_B1_map_percent(n));
        b11 = min(abs(B1D(3,:)-resized_B1_map_percent(n)));
        indices(1,:)= find(diff == b11);
    
        [c(n,1),idx2(n,1)] = max(match_images(n,:) * conj(D.magnetization(:,indices)), [], 2);
        idx2(n,1) = indices (idx2(n,1));
    end
else
    % Dictionary match on singular value images
    for n = 1 : N^2
        [c(n,1),idx2(n,1)] = max(match_images(n,:) * conj(D.magnetization), [], 2);
    end
end 

Qmaps = reshape(D.lookup_table(idx2,:),[[N N 1], size(D.lookup_table,2)]);
Qmaps = cat(numel(size(Qmaps)),Qmaps,reshape(c ./ D.normalization(idx2).',[N N]));
% PD    = c ./ D.normalization(idx2).';
% PD    = reshape(PD, [N N]);

% Do some background masking
mask = sum(abs(reshape(svd_images,[[N N] R])),3);
mask(mask < .4* max(mask(:))) = 0;
mask (mask ~= 0) = 1;
Qmaps = bsxfun(@times,Qmaps,mask);

Qmaps = flip (Qmaps);
% figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 max(max(Qmaps(:,:,1,1)))])
% figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 max(max(Qmaps(:,:,1,2)))])
% max(max(Qmaps(:,:,1,1)))
% max(max(Qmaps(:,:,1,2)))

figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 3000]);
figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 1500]); 
% figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 3000]); colormap hot;
% figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 1500]); colormap turbo;