function [Tm, A0n, T1, A1n, T2, A2n, IRF, decay_orig_offset_rmvd, decay_orig_norm, decay_reconv_opt, Chi_Square_opt, fit_opt_XY, f_opt] = Reconv_fit(time, decay_orig, SHG, Shift, fig)

%%
%     function ymodel = reconv(time, a, b, c, d, irf)
%         L = length(time);
%         y = a*exp(b*time)+c*exp(d*time);
%         ymodel_full = conv(y, irf, 'full');
%         ymodel = ymodel_full(1:L);
%     end
%%
L = size(decay_orig,1);

offset_start_frac = 0.85;
offset_end_frac = 0.95;


offset = mean(decay_orig(round(offset_start_frac*L):round(offset_end_frac*L)));    % the background level is calculated from the early part of the decay histogram

% offset = min(decay_orig);    % minimal offset subtraction
%%%
offset = 0;
%%%


decay = decay_orig - offset;
decay(decay<0) = 0;
decay = decay_orig;

% decay(1:round(offset_end_frac*L)) = 0;
% decay(round(0.95*L):end) = 0;       % only for B&H decays that end at 9.5 ns

decay_orig_offset_rmvd = decay;
decay_orig_norm = decay_orig/sum(decay_orig);  % normalize the decay by its sum
offset_frac = offset/sum(decay);               % offset as a fraction of the decay sum (to be added to the reconvolved decay later)
decay = decay/sum(decay);                      % normalize the offset-removed decay by its sum

%%

IRF = SHG/max(SHG(:));
IRFinterp = interp1((1:1:size(IRF,1))', IRF, (1:0.1:size(IRF,1))');  % 10 fold interpolation of the IRF decay histogram
IRFinterp = circshift(IRFinterp, round(10*Shift));      % apply desired shift
IRF = interp1((1:1:size(IRFinterp,1))', IRFinterp, (1:10:size(IRFinterp,1))');  % 1/10 fold - bring back to original time resolution
%IRF = IRF/sqrt(sum(IRF.^2));



%%
warning('off');


Chi_Square_opt = Inf;
    

% b0 = [0.02;-2;0.02;-0.4; 0];    % starting parameter set                                                         %%%
% lb = [0; -100; 0; -100; 0];                                                                                        %%%
% ub = [1000; -0.1; 1000; -0.1; 1000];                                                                                %%%
% 
% fo = fitoptions('Method','NonlinearLeastSquares',...
%     'Algorithm', 'Levenberg-Marquardt',...
%     'Lower',lb,...
%     'Upper',ub,...
%     'StartPoint',b0);
% ft = fittype('reconv(time, a, b, c, d, e, irf)', 'independent', {'time', 'irf'}, 'dependent', {'decay'},'options',fo);
% 
% f2 = fit([time,IRF], decay, ft);
% b = [f2.e; f2.a; -f2.b; f2.c; -f2.d];        
% 
% decay_reconv = f2([time,IRF]); %decay_reconv = decay_reconv + offset_frac*sum(decay_reconv(1:L)); 


% b0 = [0.02;-2;0.00;-0.0001];    % starting parameter set for 1-exp fit                                        %%%
% lb = [0; -10; 0; -0.0001];                                                                                    %%%
% ub = [1000; -0.1; 0; -0.0001];                                                                                %%%

b0 = [0.02;-2;0.02;-0.4];    % starting parameter set for 2-exp fit                                           %%%
% lb = [0; -10; 0; -10];                                                                                        %%%
lb = [0; -4; 0; -4];                                                                                        %%%
% ub = [1000; -0.1; 1000; -0.1];                                                                                %%%
ub = [1000; -1.5; 1000; -0.1];                                                                                %%%
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Algorithm', 'Levenberg-Marquardt',...
    'Robust', 'Bisquare',...
    'Display','off',...
    'Lower',lb,...
    'Upper',ub,...
    'StartPoint',b0);
ft = fittype('reconv(time, a, b, c, d, irf)', 'independent', {'time', 'irf'}, 'dependent', {'decay'}, 'coefficients', {'a','b','c','d'}, 'options',fo);

% b0 = [0.02;-2;0.02;-0.4; 0];    % starting parameter set  shift                                                  %%%
% lb = [0; -10; 0; -10; 0];                                                                                        %%%
% ub = [1000; -0.1; 1000; -0.1; 50];                                                                                %%%
% fo = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',lb,...
%     'Upper',ub,...
%     'StartPoint',b0);
% ft = fittype('reconvshift(time, a, b, c, d, shift, irf)', 'independent', {'time', 'irf'}, 'dependent', {'decay'}, 'coefficients', {'a','b','c','d','shift'},'options', fo);

% b0 = [0.02;-2;0.00;-0.0001; 0];    % starting parameter set  offset   1-exp                                       %%%
% lb = [0; -10; 0; -0.0001; 0];                                                                                        %%%
% ub = [1000; -0.1; 0; -0.0001; 1000];                     

b0 = [0.02;-2;0.02;-0.4; 0];    % starting parameter set  offset   2-exp                                       %%%
lb = [0; -10; 0; -10; 0];                                                                                        %%%
ub = [1000; -0.1; 1000; -0.1; 1000];                                                                                %%%
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Algorithm', 'Levenberg-Marquardt',...
    'Robust', 'Bisquare',...
    'Display','off',...
    'Lower',lb,...
    'Upper',ub,...
    'StartPoint',b0);
ft = fittype('reconvshift(time, a, b, c, d, e, irf)', 'independent', {'time', 'irf'}, 'dependent', {'decay'}, 'coefficients', {'a','b','c','d','e'},'options', fo);
% 

f2 = fit([time,IRF], decay, ft);   %%%

% b = [offset_frac*sum(f2(time,IRF)); f2.a; -f2.b; f2.c; -f2.d];   
b = [f2.e; f2.a; -f2.b; f2.c; -f2.d];   

decay_reconv = f2([time,IRF]); decay_reconv = decay_reconv + offset_frac*sum(decay_reconv(1:L)); 
decay_reconv_norm = decay_reconv(1:L)/sum(decay_reconv(1:L));
    
Chi_Square_temp = (((decay_reconv_norm(1:L)-decay_orig_norm(1:L)).^2)...
    ./(decay_orig_norm(1:L)));
Chi_Square_red = (sum(decay_orig(:))/sum(decay_orig_norm(:)))*...
    sum(Chi_Square_temp(isfinite(Chi_Square_temp)))/L;

%%
Chi_Square_opt = Chi_Square_red;
b_opt = b;
decay_reconv_opt = decay_reconv_norm;
f_opt = f2;
fit_opt_XY = [time, f_opt(time,IRF)];

%%%%%%%
%%% Single-exponential decay fit
%{
T1 = 1/b_opt(3); A1n = b_opt(2)/(b_opt(1)+b_opt(2));
T2 = 0; A2n = 0;

A0n = b_opt(1)/(b_opt(1)+b_opt(2));
Tm = T1;
%}
%%% Bi-exponential decay fit
%
if b_opt(3) >= b_opt(5)
    T1 = 1/b_opt(3); A1n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4));
    T2 = 1/b_opt(5); A2n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4));
    
else
    T1 = 1/b_opt(5); A1n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4));
    T2 = 1/b_opt(3); A2n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4));
end
A0n = b_opt(1)/(b_opt(1)+b_opt(2)+b_opt(4));
Tm = A1n*T1+A2n*T2;
%}
%%%%%%%
%{
[~,sort_ind] = sort([b_opt(3), b_opt(5), b_opt(7)], 'descend');
switch sort_ind(1)
    case [1,2,3]
        T1 = 1/b_opt(3); A1n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T2 = 1/b_opt(5); A2n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T3 = 1/b_opt(7); A3n = b_opt(6)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
    case [1,3,2]
        T1 = 1/b_opt(3); A1n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T2 = 1/b_opt(7); A2n = b_opt(6)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T3 = 1/b_opt(5); A3n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
    case [2,1,3]
        T1 = 1/b_opt(5); A1n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T2 = 1/b_opt(3); A2n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T3 = 1/b_opt(7); A3n = b_opt(6)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
    case [2,3,1]
        T1 = 1/b_opt(5); A1n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T2 = 1/b_opt(7); A2n = b_opt(6)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T3 = 1/b_opt(3); A3n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
    case [3,1,2]
        T1 = 1/b_opt(7); A1n = b_opt(6)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T2 = 1/b_opt(3); A2n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T3 = 1/b_opt(5); A3n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
    case [3,2,1]
        T1 = 1/b_opt(7); A1n = b_opt(6)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T2 = 1/b_opt(5); A2n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
        T3 = 1/b_opt(3); A3n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
end
A0n = b_opt(1)/(b_opt(1)+b_opt(2)+b_opt(4)+b_opt(6));
Tm = (A1n*T1+A2n*T2)/(A1n+A2n);
%}
%%%%%%%
warning('on');

%% Plot all decays and fitted curves, i.e., deconvolved decay, and [reconvolutions of deconvolved decay with IRF]
%{
%figure(fig), plot(time, decay, 'LineWidth', 2, 'Color', 'b'); hold on; ax=gca;ax.FontSize=25;ax.LineWidth=2; ylim([-0.005 0.075]); ylabel('Norm. Count', 'Interpreter', 'latex', 'FontSize', 25); xlabel('time [ns]', 'Interpreter', 'latex', 'FontSize', 25);
figure(fig), plot(time, decay, 'LineWidth', 2, 'Color', 'b'); hold on; ax=gca;ax.FontSize=25;ax.LineWidth=2; ylim([-0.005 max(IRF(:))/50*sum(IRF)]); ylabel('Norm. Count', 'Interpreter', 'latex', 'FontSize', 25); xlabel('time [ns]', 'Interpreter', 'latex', 'FontSize', 25);
figure(fig), plot(time, IRF/80*sum(IRF), 'LineWidth', 2, 'Color', 'g');        % plot IRF
figure(fig), plot(time, decay_deconv_final_opt, 'LineWidth', 2, 'Color', 'm');
figure(fig), plot(time, decay_reconv_opt(1:L), 'LineWidth', 2, 'Color', 'r');
% figure(fig), plot (time(2:end-round(0.24*L)-round(Shift)), modelfun(b_opt,time(2:end-round(0.24*L)-round(Shift))), 'Color', 'm');
figure(fig), plot (time(2:end), modelfun(b_opt,time(2:end)), 'Color', 'm');
% figure(fig), plot (time, decay_reconv2_opt(1:L)/sum(decay_reconv2_opt(round(0.164*L):L-round(0.085*L))), 'LineWidth', 1.5, 'Color', 'c');
figure(fig), plot (time, decay_reconv2_opt(1:L)/sum(decay_reconv2_opt), 'LineWidth', 1.5, 'Color', 'c');

disp(' ');
disp(char(strcat('T0 =', {' '},num2str(inf),'    ns', '    A0% =', {' '},num2str(100*A0n))));
disp(char(strcat('T1 =', {' '},num2str(T1, '%05.4f'),' ns', '    A1% =', {' '},num2str(100*A1n, '%05.4f'))));
disp(char(strcat('T2 =', {' '},num2str(T2, '%05.4f'),' ns', '    A2% =', {' '},num2str(100*A2n, '%05.4f'))));
%disp(char(strcat('p-value =', {' '}, num2str(max(mdl.Coefficients{:,4})))));
disp(char(strcat('Tm =', {' '},num2str((A1n*T1+A2n*T2), '%05.4f'),' ns')));
disp(char(strcat('Chi2_red  =', {' '},num2str(Chi_Square_red_opt))));
disp(char(strcat('Chi2_cyan =', {' '},num2str(Chi_Square_cyan_opt))));
%}
end
