%   FLIM analysis code for light-sheet SPAD array (FLIMera HORIBA) data
%   for the manuscript "Light sheet autofluorescence lifetime imaging with a single photon avalanche diode array" 
%   submitted to the Journal of Biomedical Optics
%   Kayvan Samimi (2023)
%
%%

all_decays = [];
bad_decays = [];
clear NADH_fits NADH_fits_array


Shift_start = 0; % in steps of 0.1 time bin

time = (0:0.0412:20)';             % FLIMera time vector (0 to 20) - in nanoseconds [ns] with 41.1 ps time bins


%%  Apply pixel delay correction and view intensity image and sum decay
data_reshaped = permute(data_cube_final, [3,1,2]);
data_reshaped = reshape(data_reshaped, size(data_cube_final,3), []);

data_reshaped2 = zeros(size(data_reshaped));
for i=1:size(data_reshaped,2)
%     data_reshaped2(:,i) = circshift(data_reshaped(:,i),-lag_map(i));
%     data_reshaped2(:,i) = circshift(data_reshaped(:,i),-lag_map2(i));  % lag map for YG and BB beads
%     data_reshaped2(:,i) = circshift(data_reshaped(:,i),-lag_map3(i));  % lag map for PANC-1 cells
    data_reshaped2(:,i) = circshift(data_reshaped(:,i),-lag_map4(i));  % lag map for coumarin and YG beads
%     data_reshaped2(:,i) = circshift(data_reshaped(:,i),0);             % no lag correction
end
data_reshaped2 = reshape(data_reshaped2, [size(data_cube_final,3),size(data_cube_final,1),size(data_cube_final,2)]);
data_cube = permute(data_reshaped2, [2,3,1]);

multiframe_data = data_cube(:,:,3611:4096);   % last 486 time bins
figure, semilogy(squeeze(sum(multiframe_data, [1 2])));

size_data = size(multiframe_data);
multiframe_data_offsetrmvd = multiframe_data - 1.15*mean(multiframe_data(:,:,415:460),3);
multiframe_data_offsetrmvd(multiframe_data_offsetrmvd<0)=0;
decays_reshaped = reshape(multiframe_data_offsetrmvd, [size_data(1)*size_data(2) size_data(3)]);
intensity_reshaped2 = sum(decays_reshaped,2);   % unbinned intensity image
% figure, imagesc(sum(multiframe_data_offsetrmvd,3)); ax=gca; ax.DataAspectRatio=[2,1,1]; colormap gray; colorbar('FontSize', 20);
intensity2 = imrotate(sum(multiframe_data_offsetrmvd,3), -90);   % unbinned intensity image
figure, imagesc(intensity2); ax=gca; ax.DataAspectRatio=[1,2,1]; colormap gray; %colorbar('FontSize', 20);

%%  Pixel Binning
Bin = 1;    % Bin factor - same convention as SPCImage
clear binned_multiframe_data;
for timebin=1:size(multiframe_data,3)
    binned_multiframe_data(:,:,timebin) = conv2(multiframe_data(:,:,timebin), ones(2*Bin+1, 2*Bin+1), 'same');
end

multiframe_data = binned_multiframe_data;                 % B&H or FLIMera data
decays_reshaped = reshape(multiframe_data, [size_data(1)*size_data(2) size_data(3)]);
figure, imagesc(sum(multiframe_data,3)); ax=gca; ax.DataAspectRatio=[2,1,1]; colorbar('FontSize', 20);
%%
scrsz = get(0,'ScreenSize');
fig = figure('Position',[round(1.1*scrsz(3)) round(0.1*scrsz(4)) round(0.8*scrsz(3)) round(0.8*scrsz(4))]);
NADH_fits={};


%%  Fit entire image (multithreaded processing)

size_data = size(multiframe_data);
num_frames = size_data(end);    % number of time bins
% if (length(size_data)==4)
%     for frame_num=1:num_frames
%         multiframe_data_gray(:,:,frame_num) = rgb2gray(multiframe_data(:,:,1:3,frame_num));
%     end
% else
%     multiframe_data_gray = multiframe_data;
% end

phasor_freq = 0.050;            % in Giga Hertz [GHz]
Shift = Shift_start; % in steps of 0.1 time bin
SHG = FLIMeraIRF(3611:4096);     % IRF
% SHG = SHG - mean(SHG(1:10));
SHG = SHG - mean(SHG(415:460)); 
SHG(SHG<0)=0;
SHG(139:end)=0;
IRF = SHG/max(SHG(:));

IRFinterp = interp1((1:1:size(IRF,1))', IRF, (1:0.1:size(IRF,1))');  % 10 fold interpolation
IRFinterp = circshift(IRFinterp, round(10*Shift));
IRF = interp1((1:1:size(IRFinterp,1))', IRFinterp, (1:10:size(IRFinterp,1))');  % 1/10 fold interpolation


%%%%%%   Plot intensity image
intensity = sum(multiframe_data,3);
intensity_image = ind2rgb(uint16(4096*((intensity-min(intensity(:)))/max(max(intensity-min(intensity(:)))))),gray(4096));
scrsz = get(0,'ScreenSize');
% f1 = figure('Position',[round(0.25*scrsz(3)) round(0.1*scrsz(4)) round(0.84*scrsz(4)) round(0.8*scrsz(4))]);
f1 = figure('Position',[round(0.025*scrsz(3)) round(0.1*scrsz(4)) round(0.7*scrsz(4)) round(0.7*scrsz(4))]);
image(intensity_image); axis equal; axis tight;
%%%%%%

intensity_reshaped = sum(decays_reshaped,2);

[M, Phi, G, S] = phasor_calculator(phasor_freq, time, decays_reshaped-mean(decays_reshaped(:,415:460),2), IRF);       % shifted IRF

intensity_reshaped(:,2)=G;
intensity_reshaped(:,3)=S;


%%%%%%%%%%   Start parallel pool
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('local',14); % 14 for 14-core processor
end
%%%%%%%%%%
clear temp;
tic
parfor pix = 1:size(decays_reshaped,1)
% for pix = 1:size(decays_reshaped,1)
    decay = (decays_reshaped(pix,:))';
    L = size(decay,1);
    if (intensity_reshaped(pix,1)<=0)||(isnan(intensity_reshaped(pix,1)))
        temp(pix, :) = NaN(1,10);                 % parfor
        continue;
    end
    %   FLIMera fit
    Shift = 12; %start shift value
    [Tm, A0n, T1, A1n, T2, A2n, IRF, decay_orig_offset_rmvd, decay_orig_norm, decay_reconv_opt, Chi2, fit_opt_XY, f_opt] = Reconv_fit(time, decay, SHG, Shift, fig);    
        ind02 = find(decay_orig_offset_rmvd>=0.9*max(decay_orig_offset_rmvd),1,'first'); % determine search range for the shift parameter
        ind01 = find(decay_orig_offset_rmvd>=0.3*max(decay_orig_offset_rmvd),1,'first');
        ind00 = max(2*ind01 - ind02, 1);
        ind1 = ind00-1+find(diff(decay_orig_offset_rmvd(ind00:ind02+1))==max(diff(decay_orig_offset_rmvd(ind00:ind02+1))),1,'first');
        ind2 = find((SHG==max(SHG)),1,'first');
        % determine search range for the shift parameter
        Shift_low = max(round(0.9*ind01, 1),1)-ind2;
        Shift_high = min(round(1.1*ind01, 1),ind02)-ind2;
        Shift_mid = round(0.5*(Shift_low+Shift_high), 1);
        
        %%%%%%      Determine optimal shift value and perform 2-exp fit using it.
        [low_Tm, low_A0n, low_T1, low_A1n, low_T2, low_A2n, low_IRF, low_decay_orig_offset_rmvd, low_decay_orig_norm, ...
            low_decay_reconv_opt, low_Chi2, low_fit_opt_XY, low_f_opt] = Reconv_fit(time, decay, SHG, Shift_low, fig);
        [mid_Tm, mid_A0n, mid_T1, mid_A1n, mid_T2, mid_A2n, mid_IRF, mid_decay_orig_offset_rmvd, mid_decay_orig_norm, ...
            mid_decay_reconv_opt, mid_Chi2, mid_fit_opt_XY, mid_f_opt] = Reconv_fit(time, decay, SHG, Shift_mid, fig);
        [high_Tm, high_A0n, high_T1, high_A1n, high_T2, high_A2n, high_IRF, high_decay_orig_offset_rmvd, high_decay_orig_norm, ...
            high_decay_reconv_opt, high_Chi2, high_fit_opt_XY, high_f_opt] = Reconv_fit(time, decay, SHG, Shift_high, fig);

        num_iter=0;    
        while (abs(Shift_high-Shift_low)>0.2)&&(num_iter<=10)
            num_iter = num_iter + 1;
            [~,sort_ind] = sort([low_Chi2, mid_Chi2, high_Chi2]);
            if (sort_ind(1)==1)||(sum(sort_ind==[2,1,3])==3)
                [high_Tm, high_A0n, high_T1, high_A1n, high_T2, high_A2n, high_IRF, high_decay_orig_offset_rmvd, high_decay_orig_norm, ...
                    high_decay_reconv_opt, high_Chi2, high_fit_opt_XY, high_f_opt] =...
                deal(mid_Tm, mid_A0n, mid_T1, mid_A1n, mid_T2, mid_A2n, mid_IRF, mid_decay_orig_offset_rmvd, mid_decay_orig_norm, ...
                mid_decay_reconv_opt, mid_Chi2, mid_fit_opt_XY, mid_f_opt);
                Shift_high = Shift_mid;
                
                Shift_mid = round(mean([Shift_low, Shift_mid]), 1);
                [mid_Tm, mid_A0n, mid_T1, mid_A1n, mid_T2, mid_A2n, mid_IRF, mid_decay_orig_offset_rmvd, mid_decay_orig_norm, ...
                    mid_decay_reconv_opt, mid_Chi2, mid_fit_opt_XY, mid_f_opt] = Reconv_fit(time, decay, SHG, Shift_mid, fig);
            else    %other outcome: (sort_ind(1)==3)||(sum(sort_ind==[2,1,3])==3)
                [low_Tm, low_A0n, low_T1, low_A1n, low_T2, low_A2n, low_IRF, low_decay_orig_offset_rmvd, low_decay_orig_norm, ...
                    low_decay_reconv_opt, low_Chi2, low_fit_opt_XY, low_f_opt] =...
                deal(mid_Tm, mid_A0n, mid_T1, mid_A1n, mid_T2, mid_A2n, mid_IRF, mid_decay_orig_offset_rmvd, mid_decay_orig_norm, ...
                mid_decay_reconv_opt, mid_Chi2, mid_fit_opt_XY, mid_f_opt);
                Shift_low = Shift_mid;
                
                Shift_mid = round(mean([Shift_high, Shift_mid]), 1);
                [mid_Tm, mid_A0n, mid_T1, mid_A1n, mid_T2, mid_A2n, mid_IRF, mid_decay_orig_offset_rmvd, mid_decay_orig_norm, ...
                    mid_decay_reconv_opt, mid_Chi2, mid_fit_opt_XY, mid_f_opt] = Reconv_fit(time, decay, SHG, Shift_mid, fig);
            end
        end
        
        [~,sort_ind] = sort([low_Chi2, mid_Chi2, high_Chi2, Chi2]);
        switch sort_ind(1)
            case 1
                [Tm, A0n, T1, A1n, T2, A2n, IRF, decay_orig_offset_rmvd, decay_orig_norm, decay_reconv_opt, Chi2, fit_opt_XY, f_opt] =...
                deal(low_Tm, low_A0n, low_T1, low_A1n, low_T2, low_A2n, low_IRF, low_decay_orig_offset_rmvd, low_decay_orig_norm, ...
                low_decay_reconv_opt, low_Chi2, low_fit_opt_XY, low_f_opt);
                Shift = Shift_low;
            case 2
                [Tm, A0n, T1, A1n, T2, A2n, IRF, decay_orig_offset_rmvd, decay_orig_norm, decay_reconv_opt, Chi2, fit_opt_XY, f_opt] =...
                deal(mid_Tm, mid_A0n, mid_T1, mid_A1n, mid_T2, mid_A2n, mid_IRF, mid_decay_orig_offset_rmvd, mid_decay_orig_norm, ...
                mid_decay_reconv_opt, mid_Chi2, mid_fit_opt_XY, mid_f_opt);
                Shift = Shift_mid;
            case 3
                [Tm, A0n, T1, A1n, T2, A2n, IRF, decay_orig_offset_rmvd, decay_orig_norm, decay_reconv_opt, Chi2, fit_opt_XY, f_opt] =...
                deal(high_Tm, high_A0n, high_T1, high_A1n, high_T2, high_A2n, high_IRF, high_decay_orig_offset_rmvd, high_decay_orig_norm, ...
                high_decay_reconv_opt, high_Chi2, high_fit_opt_XY, high_f_opt);
                Shift = Shift_high;
            otherwise
                %keep the fit using user input Shift value.
        end
    
    % intensity_reshaped(pix, 4:12) = [Tm, A0n, T1, A1n, T2, A2n, Chi2_red, Chi2_cyan, f_cutoff_ind];       % for
%     temp(pix, :) = [Tm, A0n, T1, A1n, T2, A2n, Chi2_red, Chi2_cyan, f_cutoff_ind, G, S];                          % parfor
    [M, Phi, G, S] = phasor_calculator(phasor_freq, time, decay_orig_offset_rmvd', IRF);
    temp(pix, :) = [Tm, A0n, T1, A1n, T2, A2n, Shift, Chi2, G, S];                          % parfor
    %{
    figure(fig), plot(time, decay_orig_norm, 'LineWidth', 2, 'Color', 'b'); hold on; ax=gca;ax.FontSize=25;ax.LineWidth=2; ylim([-0.005 max(IRF(:))/50*sum(IRF)]); ylabel('Norm. Count', 'Interpreter', 'latex', 'FontSize', 25); xlabel('time [ns]', 'Interpreter', 'latex', 'FontSize', 25);
    figure(fig), plot(time, IRF/50*sum(IRF), 'LineWidth', 2, 'Color', 'g');        % plot IRF
    figure(fig), plot(time, decay_deconv_final_opt, 'LineWidth', 2, 'Color', 'm');
    figure(fig), plot(time, decay_reconv_opt(1:L), 'LineWidth', 2, 'Color', 'r');
    figure(fig), plot (fit_opt_XY(:,1), fit_opt_XY(:,2), 'Color', 'k');
    figure(fig), plot (time, decay_reconv2_opt, 'LineWidth', 1.5, 'Color', 'c');
        
    disp(char(strcat('Chi2_red ='," ",num2str(Chi2_red,'%05.4f'), "  ", 'Chi2_cyn ='," ",num2str(Chi2_cyan,'%05.4f'), "  ",...
        'A1% ='," ",num2str(100*A1n/(A1n+A2n),'%04.2f'), "  ", 'Tm ='," ",num2str(1000*Tm,'%04.0f'),'[ps]', "   ",...
        'T1 ='," ",num2str(1000*T1,'%04.0f'),'[ps]', "  ",  'T2 ='," ",num2str(1000*T2,'%04.0f'),'[ps]', "  ", 'f_cutoff_ind =',...
        " ",num2str(f_cutoff_ind), " ", 'Shift ='," ",num2str(Shift,'%5.1f'), " ", 'num_iter ='," ",num2str(num_iter), " ", 'count=',num2str(sum(decay)))));
    pause(0.01);
    hold off;
    %}
end
elapsedTime = toc
intensity_reshaped = cat(2, intensity_reshaped, temp);          % parfor


imgrgb = ind2rgb(uint16(1024*(reshape((intensity_reshaped(:,4)>0 & intensity_reshaped(:,4)<=5 & intensity_reshaped(:,1)>=200).*intensity_reshaped(:,4), [size_data(1),size_data(2)])-0.0)/(3.5-0.0)), jet(1024));

imgrgb_scaled = 3*repmat(reshape(intensity_reshaped2(:,1)/max(intensity_reshaped2(:,1)), [size_data(1),size_data(2)]), [1 1 3]).*imgrgb;

imgrgb_scaled = imrotate(imgrgb_scaled, -90);       %%%%%

fig2 = figure('Position',[round(0.525*scrsz(3)) round(0.1*scrsz(4)) round(0.7*scrsz(4)) round(0.7*scrsz(4))]);
% imshow(imgrgb_scaled); ax=gca; ax.DataAspectRatio=[2,1,1];
imshow(imgrgb_scaled); ax=gca; ax.DataAspectRatio=[1,2,1]; %%%%%
colormap jet;
caxis([0.0, 3.5]); colorbar('FontSize', 20);
%}

% figure, hist(intensity_reshaped(:,4),1000)
% figure, imagesc(reshape(intensity_reshaped(:,4), [size_data(1),size_data(2)])); ax=gca; ax.DataAspectRatio=[2,1,1]; colorbar;         % Tm
% figure, imagesc(imrotate(reshape(intensity_reshaped(:,4), [size_data(1),size_data(2)]),-90)); ax=gca; ax.DataAspectRatio=[1,2,1]; colorbar;         % Tm
% figure, imagesc(imrotate(sum(multiframe_data_offsetrmvd,3),-90)); ax=gca; ax.DataAspectRatio=[1,2,1]; colormap gray; colorbar;  % Intensity
% figure, imagesc(reshape(100*intensity_reshaped(:,7), [size_data(1),size_data(2)])); ax=gca; ax.DataAspectRatio=[2,1,1]; colorbar;     % Alpha1%
% figure, imagesc(imrotate(reshape(100*intensity_reshaped(:,7), [size_data(1),size_data(2)]),-90)); ax=gca; ax.DataAspectRatio=[1,2,1]; colorbar;         % Alpha1%
% figure, imagesc(reshape(1000*intensity_reshaped(:,6), [size_data(1),size_data(2)])); ax=gca; ax.DataAspectRatio=[2,1,1]; colorbar;    % T1
% figure, imagesc(reshape(1000*intensity_reshaped(:,8), [size_data(1),size_data(2)])); ax=gca; ax.DataAspectRatio=[2,1,1]; colorbar;    % T2
% figure, imagesc(imrotate(reshape(1000*intensity_reshaped(:,8), [size_data(1),size_data(2)]),-90)); ax=gca; ax.DataAspectRatio=[1,2,1]; colorbar;    % T2
% figure, imagesc(reshape(intensity_reshaped(:,11), [size_data(1),size_data(2)])); ax=gca; ax.DataAspectRatio=[2,1,1]; colorbar;        % Chi2

%%  Upscale FLIM resolution using CMOS camera intensity image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tm
FLIMera = imrotate(reshape(intensity_reshaped2(:,1), [size_data(1),size_data(2)]), -90);
figure, imagesc(double(FLIMera)); colormap gray; ax=gca; ax.DataAspectRatio=[1,2,1];

F = griddedInterpolant({1:192,1:128},FLIMera);
q1 = (1:6.5/9.2:192)';
q2 = (1:6.5/18.4:128)';
FLIMeraInterp = F({q1,q2});

% FLIMera_Tm = imrotate(reshape(intensity_reshaped(:,4), [size_data(1),size_data(2)]), -90);
tm_reshaped = (intensity_reshaped(:,7).*intensity_reshaped(:,6)+intensity_reshaped(:,9).*intensity_reshaped(:,8))./(intensity_reshaped(:,7)+intensity_reshaped(:,9));
tm_reshaped(isnan(tm_reshaped))=0;
FLIMera_Tm = imrotate(reshape(tm_reshaped, [size_data(1),size_data(2)]), -90);
ti_reshaped = (intensity_reshaped(:,7).*intensity_reshaped(:,6).^2+intensity_reshaped(:,9).*intensity_reshaped(:,8).^2)./(intensity_reshaped(:,7).*intensity_reshaped(:,6)+intensity_reshaped(:,9).*intensity_reshaped(:,8));
ti_reshaped(isnan(ti_reshaped))=0;
FLIMera_Tm = imrotate(reshape(ti_reshaped, [size_data(1),size_data(2)]), -90);
% figure, imagesc(double(FLIMera_Tm)); colormap jet; colorbar('FontSize', 20); ax=gca; ax.DataAspectRatio=[1,2,1];

F = griddedInterpolant({1:192,1:128},FLIMera_Tm);
q1 = (1:6.5/9.2:192)';
q2 = (1:6.5/18.4:128)';
FLIMera_TmInterp = F({q1,q2});
figure, imagesc(double(FLIMera_TmInterp)); colorbar('FontSize', 20); axis equal tight

imgrgb = ind2rgb(uint16(1024*(FLIMera_TmInterp-0.0)/(3.0-0.0)), jet(1024));   %%%%%%%%%%%% Tm
%%%% CMOS intensity-weighted
% imgrgb_scaled = repmat(CMOSreg/max(CMOSreg(:)), [1 1 3]).*imgrgb;
imgrgb_scaled1 = repmat(CMOSreg/23500, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled1); colormap jet; caxis([0.0, 3.0]); colorbar('FontSize', 20);
%%%% FLIMera intensity-wighted
% imgrgb_scaled = repmat(FLIMeraInterp/max(FLIMeraInterp(:)), [1 1 3]).*imgrgb;
imgrgb_scaled2 = repmat(FLIMeraInterp/3550, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled2); colormap jet; caxis([0.0, 3.0]); colorbar('FontSize', 20);
%%%% CMOS intensity-weighted & log scaled
% temp = log(double((CMOSreg-min(CMOSreg(:)))/(max(CMOSreg(:))-min(CMOSreg(:)))+1)); %temp = temp - min(temp(:));
CMOSregLog = 2*log(double((CMOSreg-1000)/(23500-1000)+1)); %temp = temp - min(temp(:));
% imgrgb_scaled = 1.25*repmat((temp-min(temp(:)))/(max(temp(:))-min(temp(:))), [1 1 3]).*imgrgb;
imgrgb_scaled3 = repmat(CMOSregLog, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled3); colormap jet; caxis([0.0, 3.0]); colorbar('FontSize', 20);

% figure, imshowpair(FLIMeraInterp,CMOSreg,'montage');
figure, imshowpair(FLIMeraInterp/3550,CMOSregLog,'montage','Scaling','none'); colorbar('FontSize', 20);
figure, imshowpair(imgrgb_scaled2,imgrgb_scaled3,'montage','Scaling','none'); colormap jet; caxis([0.0, 3.0]); colorbar('FontSize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% alpha1%

% FLIMera_a1 = imrotate(reshape(100*intensity_reshaped(:,7), [size_data(1),size_data(2)]), -90); %%%% alpha1%
FLIMera_a1 = imrotate(reshape(100*(intensity_reshaped(:,7)./(intensity_reshaped(:,7)+intensity_reshaped(:,9))), [size_data(1),size_data(2)]), -90); %%%% alpha1%
% figure, imagesc(double(FLIMera_a1)); colorbar; ax=gca; ax.DataAspectRatio=[1,2,1];

F = griddedInterpolant({1:192,1:128},FLIMera_a1);
q1 = (1:6.5/9.2:192)';
q2 = (1:6.5/18.4:128)';
FLIMera_a1Interp = F({q1,q2});
figure, imagesc(double(FLIMera_a1Interp)); colorbar; axis equal tight
imgrgb = ind2rgb(uint16(1024*(FLIMera_a1Interp-0)/(100-0)), jet(1024));   %%%%%%%%%%%% alpha1%
%%%% CMOS intensity-weighted
% imgrgb_scaled = repmat(CMOSreg/max(CMOSreg(:)), [1 1 3]).*imgrgb;
imgrgb_scaled1 = repmat(CMOSreg/23500, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled1); colormap jet; caxis([0, 100]); colorbar('FontSize', 20);
%%%% FLIMera intensity-wighted
% imgrgb_scaled = repmat(FLIMeraInterp/max(FLIMeraInterp(:)), [1 1 3]).*imgrgb;
imgrgb_scaled2 = repmat(FLIMeraInterp/3550, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled2); colormap jet; caxis([0, 100]); colorbar('FontSize', 20);
%%%% CMOS intensity-weighted & log scaled
% temp = log(double((CMOSreg-min(CMOSreg(:)))/(max(CMOSreg(:))-min(CMOSreg(:)))+1)); %temp = temp - min(temp(:));
CMOSregLog = 2*log(double((CMOSreg-1000)/(23500-1000)+1)); %temp = temp - min(temp(:));
% imgrgb_scaled = 1.25*repmat((temp-min(temp(:)))/(max(temp(:))-min(temp(:))), [1 1 3]).*imgrgb;
imgrgb_scaled3 = repmat(CMOSregLog, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled3); colormap jet; caxis([0, 100]); colorbar('FontSize', 20);

figure, imshowpair(imgrgb_scaled2,imgrgb_scaled3,'montage','Scaling','none'); colormap jet; caxis([0, 100]); colorbar('FontSize', 20);

%%
filename = 'Cyanide_Reconv_Analysis_time_0.mat';
save(filename, 'data_cube_Cy_time_0_min', 'intensity_reshaped', 'intensity_reshaped2', 'CMOSreg', 'CMOSregLog', 'FLIMera', 'FLIMeraInterp', 'FLIMera_Tm', 'FLIMera_TmInterp', 'FLIMera_a1', 'FLIMera_a1Interp');

%%
FLIMera_T1 = imrotate(reshape(1000*intensity_reshaped(:,6), [size_data(1),size_data(2)]), -90); %%%% T1
% figure, imagesc(double(FLIMera_T1)); colorbar; ax=gca; ax.DataAspectRatio=[1,2,1];

F = griddedInterpolant({1:192,1:128},FLIMera_T1);
q1 = (1:6.5/9.2:192)';
q2 = (1:6.5/18.4:128)';
FLIMera_T1Interp = F({q1,q2});
figure, imagesc(double(FLIMera_T1Interp)); colorbar; axis equal tight
imgrgb = ind2rgb(uint16(1024*(FLIMera_T1Interp-0)/(1000-0)), jet(1024));   %%%%%%%%%%%% T1
%%%% CMOS intensity-weighted
% imgrgb_scaled = repmat(CMOSreg/max(CMOSreg(:)), [1 1 3]).*imgrgb;
imgrgb_scaled1 = repmat(CMOSreg/23500, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled1); colormap jet; caxis([0, 1000]); colorbar('FontSize', 20);
%%%% FLIMera intensity-wighted
% imgrgb_scaled = repmat(FLIMeraInterp/max(FLIMeraInterp(:)), [1 1 3]).*imgrgb;
imgrgb_scaled2 = repmat(FLIMeraInterp/3550, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled2); colormap jet; caxis([0, 1000]); colorbar('FontSize', 20);
%%%% CMOS intensity-weighted & log scaled
% temp = log(double((CMOSreg-min(CMOSreg(:)))/(max(CMOSreg(:))-min(CMOSreg(:)))+1)); %temp = temp - min(temp(:));
CMOSregLog = 2*log(double((CMOSreg-1000)/(23500-1000)+1)); %temp = temp - min(temp(:));
% imgrgb_scaled = 1.25*repmat((temp-min(temp(:)))/(max(temp(:))-min(temp(:))), [1 1 3]).*imgrgb;
imgrgb_scaled3 = repmat(CMOSregLog, [1 1 3]).*imgrgb;
figure, imshow(imgrgb_scaled3); colormap jet; caxis([0, 1000]); colorbar('FontSize', 20);


figure, imshowpair(imgrgb_scaled2,imgrgb_scaled3,'montage','Scaling','none'); colormap jet; caxis([0, 1000]); colorbar('FontSize', 20);
