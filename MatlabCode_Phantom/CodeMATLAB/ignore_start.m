clear;
%% Magnetic Resonance Fingerprinting reconstruction code

% ... Using data from Philips scanner 

%% Initialize paths, check for bart toolbox and load data

% Required: BART (installed in WSL)
bart_loc = '/mnt/c/Users/lucya/MSC_PROJECT/bart/';

% matlab folder of BART contains a script to call BART in WSL
addpath(genpath('C:/Users/lucya/MSC_PROJECT/bart/matlab/'));
setenv('TOOLBOX_PATH',bart_loc);

% check if bart toolbox is installed
if ~isempty(which('bart')); else;disp('>> BART toolbox not found in path. Please check BART path.');end

% add the "utils" subfolder to the PATH:
addpath(genpath('C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/CodeMATLAB/utils'));

% folder that contains the raw data:
folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/';

%% Load data  
data = loadRawKspacePhilips([folder 'raw_001.list']);

start = 101;

% raw data and image dimensions (according to .list file):
Nkx = abs(data.kspace_properties.kx_range(2) - data.kspace_properties.kx_range(1)) + 1;
Nky = abs(data.kspace_properties.ky_range(2) - data.kspace_properties.ky_range(1)) + 1;
Nx = abs(data.kspace_properties.X_range(2) - data.kspace_properties.X_range(1)) + 1;
Ny = abs(data.kspace_properties.Y_range(2) - data.kspace_properties.Y_range(1)) + 1;
Nc  = numel(unique(data.chan));
dyn = 1;  
% dyn=numel(unique(data.dyn)); 
neco = data.kspace_properties.number_of_echoes;
slices=numel(unique(data.loca));

distanceSlice = 0;

% extract noise adjustment data:
noise = complex(zeros(numel(data.complexdata{find(cell2mat(data.typ) == 'NOI', 1)}), Nc));
for ii = 1:numel(data.complexdata)
    if strcmp(data.typ{ii}, 'NOI')
        noise(:, data.chan(ii) + 1) = data.complexdata{ii};
    end
end

% put raw data in more useful structure:
raw = complex(zeros(Nkx,Nc,Nky,dyn));
for ii = 1:(numel(data.complexdata))
    if strcmp(data.typ{ii}, 'STD')
        raw(:,  data.chan(ii)+1, data.ky(ii) +1 ,data.dyn(ii)+1) = data.complexdata{ii};
    end
end

%% Data preprocessing
U.RawKspaceData = squeeze(raw(:,:,:,1:dyn));
U.RawKspaceData = permute(U.RawKspaceData,[1 3 4 2]);
U.RawKspaceData = noise_prewhitening(U.RawKspaceData,noise);
kdim = c12d(size(U.RawKspaceData));

sliceOrientation = 1; % Coronal

kx_GIRF = fread(fopen([folder 'kx_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 
ky_GIRF = fread(fopen([folder 'ky_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 
kz_GIRF = fread(fopen([folder 'kz_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 
b0_GIRF = fread(fopen([folder 'b0_gve_001_girf.bin']),[Nkx*Nky,dyn],'float64'); 

b0_GIRF = b0_GIRF.*(2*pi*42.57746778); 
trajx_GIRF = reshape(kx_GIRF,[Nkx,Nky*dyn]);
trajy_GIRF = reshape(ky_GIRF,[Nkx,Nky*dyn]);
trajz_GIRF = reshape(kz_GIRF,[Nkx,Nky*dyn]);
b0_GIRF = reshape(b0_GIRF,[Nkx,Nky*dyn]);

trajx_GIRF = trajx_GIRF * 2 * pi * distanceSlice; 

N = 124;
comb_raw = zeros(kdim(1),Nc,kdim(2)*dyn);
location = 0;
step = 1000;

U.RawKspaceData = permute(U.RawKspaceData,[1 4 2 3]); 
for d = 1:dyn
    comb_raw(:,:,location+1:step) = squeeze(U.RawKspaceData(:,:,:,d));
    location = location+1000;
    step = step+1000;
end

U.RawKspaceData = permute(comb_raw,[1 3 4 2]);

U.RawKspaceData = U.RawKspaceData(:,start:1000*dyn,:,:);     % crops k-space data
b0_GIRF     = b0_GIRF(:,start:1000*dyn);      % crops the B0-induced phase error correction values
trajx_GIRF  = trajx_GIRF(:,start:1000*dyn);     % crops the x/y/z-component of the k-space trajectory (after GIRF correction)
trajy_GIRF  = trajy_GIRF(:,start:1000*dyn);   
trajz_GIRF  = trajz_GIRF(:,start:1000*dyn);


trajGIRF(1,:,:) = trajz_GIRF;
trajGIRF(2,:,:) = trajy_GIRF;
trajGIRF(3,:,:) = 0;
trajGIRF = repmat(trajGIRF,[1 1 1 1 kdim(6)]);
dcf = permute(radial_density(trajGIRF(:,:,:,:,1)),[1 2 3 5 4]);

U.RawKspaceData = U.RawKspaceData .* exp(-1i * b0_GIRF);
if sliceOrientation == 1 
    U.RawKspaceData = U.RawKspaceData .* exp(-1i * trajx_GIRF);
elseif sliceOrientation == 2
    U.RawKspaceData = U.RawKspaceData .* exp(-1i * trajz_GIRF);
elseif sliceOrientation == 3
    U.RawKspaceData = U.RawKspaceData .* exp(-1i * trajy_GIRF);
else 
    disp('An error occurred')
end

B1correct = 0;
nufft = bart(['nufft -a -d ',num2str(N),':',num2str(N),':1'],0.25*trajGIRF,permute(bsxfun(@times,dcf,U.RawKspaceData),[3 1 2 4]));
imgForCSM = bart('fft 3',nufft);
csm = bart('ecalib -m1',imgForCSM);


%% Dictionary 

% loading dictionary from MATLAB file (5865 entries):
% Each entry is time-series signal for specific combination of parameters
% load(fullfile(folder, 'D_Phantom2025_invEff096_SPinf_norefocusingTEadj_576InvTime_1000RF_10mm_101iso_0.mat'));
% load(fullfile(folder, 'D_IDX_SP_Phantom2025.mat'));     % apparently this is also needed - these are likely the T1, T2 (and B1 (?)) values for each dictionary entry
dict_folder = 'C:/Users/LucyA/MSC_PROJECT/Dictionary/dict';
load(fullfile(dict_folder, 'blochdict_10mm_701_short.mat'));

% alternatively: the dictionary with B1:
% load(fullfile(folder, 'Dictionary_B1_SP.mat'));


%dict0 = dict0(1:1000*dyn,:);    % ... First dyn*1000 time points
dict0 = dict0(start:1000*dyn,:);  % now only 900 timepoints
[~,s,~]=svd((dict0),'econ');    % ... SVD: captures dominant signal patterns
result_s    = diag(s);
NRJ= zeros(100,1);
for R = 1:100    % ... Energy retained in first R modes
    NRJ(R,1) = sum(result_s(1:R))/sum(result_s);
end

R = find(NRJ>=0.999,1);    % Pick # of modes to keep 99.9% of signal energy
disp(R)



% dict must be under the form timeFrame * parameters combinations
[D] = compression_mrf_dictionary(dict0,idx, R);       % ... Compressing dictionary
%%D = MRF_dictionary_umcu(idx', R, double(sum(dict,3))); % epg       % ???
% First write all the cfl/hdr files in the correct dimensions for BART
writecfl('csm',single(csm));
writecfl('u',permute(D.u,[3 4 5 6 7 1 2]));
writecfl('traj',permute(.25 * trajGIRF,[1 2 5 4 6 3]))
writecfl('dcf',permute(sqrt(dcf),[5 1 3 4 6 2]));
writecfl('kdata',permute(U.RawKspaceData,[3 1 5 4 6 2]));

% LR inversion parameters
lambda = 0.05;    % ...Regularisation weight
n_iter = 50;    % ... Iteration count 

% Reconstruct singular value images with BART
tic 
cmd = sprintf('pics -d1 -e -l1 -r %f -i %d -p dcf -t traj -B u kdata csm', lambda, n_iter);
recon = bart(cmd);   % ... Outputs complex valued image of shape [N, N, R]

writecfl('recon',single(recon));
toc 

svd_images= reshape(readcfl('recon'),N,N,R);
match_images = reshape(svd_images,N^2,R);

c   = zeros([size(match_images,1) 1],'single');
idx2 = zeros([size(match_images,1) 1],'single');

if B1correct ~= 1 && B1correct ~= 0
   disp("Wrong entry");
elseif B1correct == 1
    B1   = [0.8:0.02:0.9, 0.91:0.01:1.09, 1.10:0.02:1.2] ;    % ... Range of B1 scaling factors
    % need to specify the folder for the B1 map...
    folderB1 = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/DICOM/';
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

recon_folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/recon_results';
save(fullfile(recon_folder, '10mm_701_short_skip100.mat'), 'svd_images', 'D', 'dict0', 'idx', 'R', 'match_images', 'idx2', 'c', 'N');

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
