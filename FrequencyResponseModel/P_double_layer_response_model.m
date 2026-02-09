clc
clear;
close all;
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultTextFontSize',12);
set(0,'DefaultTextFontName','Times');
set(0,'DefaultLineLineWidth',0.5);
Color=[1 0 0; 0 0.5 0; 0 0 0.8];
%% Graph rule 
xlabel_name='\itf \rm[Hz]';
ylabel_name_Gain='Gain [dB]';
ylabel_name_Phase='Phase [deg]';
x_lim=[100 10000];
y_lim_gain=[-3.5 0.5];
y_lim_phase=[-40 5];
%% Data import and sorting
Data = xlsread('Measurement Data.xlsx','A1:C20');
Data(1, :) = [];
[~, idx] = sort(Data(:,1));  
Data = Data(idx,:);           
outdir='fig_Typical_Modle\';
title1='Gain';
title2='Phase';
mkdir(outdir); 
figsize=[300 100 400 350];
marker_size=8;
%% Initial values 
% Luminescent lifetime
tau =[7.04e-6];
%Diffusivity coefficient
D1_ini=1.96d-7;                                        
% Frequency
f=Data(:,1); 

%% Parse gain and phase data
[m,n]=size(Data);  
Data(:,1)=[];
for i=1:(n-1)/2
    Gain(:,i)=Data(:,(i-1)*2+1);
    Phase(:,i)=Data(:,(i-1)*2+2);

end 

%% Typical single (Kameda parameters)
% Diffusivity coefficient
design_variable_up = 1.96d-7; % Upper layer
design_variable_bottom = design_variable_up;
% Coating thickness
dl1=3.0d-6;  % Upper layer
dl2=0;     % Bottom layer 
% Hidden factor
hf_ini=0;  

% Predict gain and phase based on optimal coefficient
for i=1:(n-1)/2
  Result_est_min_single(:,:,i)=F_transfer_function_doublelayer(f,dl1,dl2,design_variable_up(1,i),...
         design_variable_bottom(1,i),hf_ini,tau(1,i));
end

%% Typical double (Pandey and Gregory parameter)
% Diffusivity coefficient
design_variable_up = 9.78d-8; % Upper layer
design_variable_bottom = 9.78d-11; % Bottom layer
% Coating thickness
dl1=3.0d-6;  % Upper layer
dl2=37.0d-6;  % Bottom layer 
% Hidden factor
hf_ini=1.67d6;  

% Predict gain and phase based on optimal coefficient
for i=1:(n-1)/2
  Result_est_min_double(:,:,i)=F_transfer_function_doublelayer(f,dl1,dl2,design_variable_up(1,i),...
         design_variable_bottom(1,i),hf_ini,tau(1,i));
end


%% figure output
% Gain attenuation
for i=1:(n-1)/2
    figure(i)
    hAx=axes;
    plot(f(:,1), Gain(:,i), 'x', ...
     'LineStyle','none', ...
     'Color', Color(1,:,:), ...
     'MarkerSize', marker_size);
    hold on;
    plot(f,Result_est_min_single(:,1,i),'Color',Color(2,:,:))
    plot(f,Result_est_min_double(:,1,i),'Color',Color(3,:,:))
    legend('Measurement', 'Single layer model', 'Double layer model', ...
           'Location', 'southwest');
    hAx.XScale='log';
    set(gcf,'Position',figsize);
    grid on;
    xlabel(xlabel_name)
    ylabel(ylabel_name_Gain)
    xticks([100 1000 10000])
    xticklabels({'100','1000','10000'})
    xlim(x_lim)
    ylim(y_lim_gain)
    saveas(gcf,[ outdir title1 num2str(i) '.fig']);
    saveas(gcf,[ outdir title1 num2str(i)  '.tif']);
    saveas(gcf,[ outdir title1 num2str(i)  '.svg']);

% Phase delay
    figure(i*10)
    hAx=axes;
    plot(f(:,1), Phase(:,i), 'x', ...
     'LineStyle','none', ...
     'Color', Color(1,:,:), ...
     'MarkerSize', marker_size);
    hold on;
    plot(f,Result_est_min_single(:,2,i),'Color',Color(2,:,:))
    plot(f,Result_est_min_double(:,2,i),'Color',Color(3,:,:))
    legend('Measurement', 'Single layer model', 'Double layer model', ...
           'Location', 'southwest');
    hAx.XScale='log';  
    set(gcf,'Position',figsize);
    grid on;
    xticks([100 1000 10000])
    xticklabels({'100','1000','10000'})
    xlabel(xlabel_name)
    ylabel(ylabel_name_Phase)
    xlim(x_lim)
    ylim(y_lim_phase)
    saveas(gcf,[ outdir title2 num2str(i)  '.fig']);
    saveas(gcf,[ outdir title2 num2str(i)  '.tif']);
    saveas(gcf,[ outdir title2 num2str(i)  '.svg']);
end