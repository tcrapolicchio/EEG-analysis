clear all
close all
clc

%% Loading and cutting signals
load('chanlocs.mat'); 

num_subjects = 6; % Number of subjects
%deriv_names = {'FP1', 'FP2', 'F3', 'F4', 'F7', 'F8', 'T3', 'T4', 'C3', 'C4', 'T5', 'T6', 'P3', 'P4', 'O1', 'O2', 'FZ', 'CZ', 'PZ'};

all_labels = {chanlocs.labels}; 
deriv_names = all_labels(1:19); 
deriv_names{1} = 'FP1'; 
deriv_names{2} = 'FP2'; 
deriv_names{17} = 'FZ'; 
deriv_names{18} = 'CZ'; 
deriv_names{19} = 'PZ'; 


num_derivations = length(deriv_names); % Number of derivations
fs = 500; % Sampling Frequency
t = 0:1/fs:60 - 1/fs; 

% Loading data
k = [1, 2, 3, 5, 7, 8];
for subj = 1:num_subjects
    data_rest(subj) = load(sprintf('Subject%02d_1', k(subj)));
    data_task(subj) = load(sprintf('Subject%02d_2', k(subj))); 
    
    % Selecting last minute of rest and first minute of task
    for deriv = 1:num_derivations
       data_rest(subj).(deriv_names{deriv}) = data_rest(subj).(deriv_names{deriv})(60001:90000); 
       data_task(subj).(deriv_names{deriv}) = data_task(subj).(deriv_names{deriv})(1:30000); 
    end
end

%% FILTERING: notch at 50Hz already present

lowCut = 0.5;
highCut = 45;

% Bandpass Filter
for subj = 1:num_subjects
    for deriv = 1:num_derivations
        [b_bp, a_bp] = butter(4, [lowCut highCut] / (fs / 2));
        data_rest(subj).(deriv_names{deriv}) = filtfilt(b_bp, a_bp, data_rest(subj).(deriv_names{deriv}));
        data_task(subj).(deriv_names{deriv}) = filtfilt(b_bp, a_bp, data_task(subj).(deriv_names{deriv}));

        % Normalization
        data_rest(subj).(deriv_names{deriv}) = (data_rest(subj).(deriv_names{deriv}) - mean(data_rest(subj).(deriv_names{deriv}))) / std(data_rest(subj).(deriv_names{deriv}));
        data_task(subj).(deriv_names{deriv}) = (data_task(subj).(deriv_names{deriv}) - mean(data_task(subj).(deriv_names{deriv}))) / std(data_task(subj).(deriv_names{deriv}));
    end
end 

%% PLOT EEG
% Plot one derivation of an EEG rest as an example
figure(1)
plot(t, data_rest(1).P3, color = 'blue')
xlabel('Time [s]')
ylabel('\mu V')
title('EEG Subject 1, derivation P3, Rest')

% Plot one derivation of an EEG task as an example
figure(2)
plot(t, data_task(1).P3, color = 'red')
xlabel('Time [s]')
ylabel('\mu V')
title('EEG Subject 1, derivation P3, Task')

%% PLOT FFT
fft_rest = fft(data_rest(1).P3);
f = linspace(0, fs, length(fft_rest)); 
figure(3)
plot(f, abs(fft_rest), color = 'blue')
xlim([0 fs/2])
xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('FFT Subject 1, derivation P3, Rest')

fft_task = fft(data_task(1).P3);
f = linspace(0, fs, length(fft_task)); 
figure(4)
plot(f, abs(fft_task), color = 'red')
xlim([0 fs/2])
xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('FFT Subject 1, derivation P3, Task')

%% Parameters definition
window_size = 10 * fs; 
overlap = 0.1 * fs; 

%% Spectrogram plot
spectrogram(data_rest(1).P3, window_size, overlap, [], fs); 
colorbar;              
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% PSD calculation

% ---------- Frequency Ranges selection ---------- %

bandwidths = [0.5 4; 4 8; 8 12; 12 30; 30 45]; 
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};

% bandwidths = [4.1 5.8; 5.9 7.4; 13 19.9; 20 25; 30 45]; %Theta1, Theta2, Beta1, Beta2
% bands = {'theta1', 'theta2', 'beta1', 'beta2'};

% ---------- Frequency Ranges selection ---------- %

num_bands = length(bands);

power_rest = zeros(num_subjects, num_derivations, num_bands); 
power_task = zeros(num_subjects, num_derivations, num_bands);

f_rest = zeros(4097, num_derivations, num_subjects);
f_task = zeros(4097, num_derivations, num_subjects);

pxx_rest = zeros(4097, num_derivations, num_subjects);
pxx_task = zeros(4097, num_derivations, num_subjects);
for subj = 1:num_subjects   
    for deriv = 1:num_derivations
        % PSD calculation
        EEG_rest_vector = data_rest(subj).(deriv_names{deriv}); 
        EEG_task_vector = data_task(subj).(deriv_names{deriv});

        [pxx_rest(:,deriv,subj), f_rest(:,deriv,subj)] = pwelch(EEG_rest_vector, hamming(window_size), overlap, [], fs); % PSD rest
        [pxx_task(:,deriv,subj), f_task(:,deriv,subj)] = pwelch(EEG_task_vector, hamming(window_size), overlap, [], fs); % PSD task

        % Check validity pxx 
        if isempty(pxx_rest(:,deriv,subj)) || isempty(f_rest(:,deriv,subj))
            error('PSD not valid for subject %d, deriv %d.', subj, deriv);
        end

        power_tot_rest(subj, deriv) = trapz(f_rest(:, deriv, subj), pxx_rest(:, deriv, subj));
        power_tot_task(subj, deriv) = trapz(f_task(:, deriv, subj), pxx_task(:, deriv, subj));

        % Power calculation for each band
        for band = 1: num_bands 
            try
                power_rest(subj, deriv, band) = bandpower(squeeze(pxx_rest(:,deriv,subj)), f_rest(:,deriv,subj), bandwidths(band, :), 'psd');
                power_task(subj, deriv, band) = bandpower(squeeze(pxx_task(:,deriv,subj)), f_task(:,deriv,subj), bandwidths(band, :), 'psd');
 
            catch ME
                disp(['Error in band ' num2str(bandwidths(band, :)) ': ' ME.message]);
            end
        end

    end
end

% Normalization for maximum value
for subj = 1:num_subjects
    max_rest(subj) = max(power_tot_rest(subj, :)); 
    max_task(subj) = max(power_tot_task(subj, :)); 

    max_tot(subj) = max(max_rest(subj), max_task(subj));
    power_rest(subj, :, :) = power_rest(subj, :, :) ./ max_tot(subj); 
    power_task(subj, :, :) = power_task(subj, :, :) ./ max_tot(subj); 
end


num_derivations = size(power_rest, 2); 
rest_tot = zeros(num_derivations, num_bands); 
task_tot = zeros(num_derivations, num_bands);

for band = 1:num_bands
    for deriv = 1:num_derivations
        rest_tot(deriv, band) = mean(power_rest(:, deriv, band));
        task_tot(deriv, band) = mean(power_task(:, deriv, band));
    end

    rest_tot(:, band) = rest_tot(:, band)';
    task_tot(:, band) = task_tot(:, band)';
end

figure(5);
colormap(parula);

number = [1 2; 3 4; 5 6; 7 8; 9 10];
for band = 1:num_bands
    max_val = max([max(rest_tot(:, band)), max(task_tot(:, band))]);
    min_val = min([min(rest_tot(:, band)), min(task_tot(:, band))]);

    % First subplot: REST
    subplot(num_bands, 2, number(band, 1));
    topoplot(rest_tot(:, band), chanlocs, 'maplimits', [min_val, max_val], 'electrodes', 'on');
    colorbar;
    title(sprintf('REST - %s', bands{band}));        

    % Second subplot: TASK
    subplot(num_bands, 2, number(band, 2));
    topoplot(task_tot(:, band), chanlocs, 'maplimits', [min_val, max_val], 'electrodes', 'on');
    colorbar;
    title(sprintf('TASK - %s', bands{band}));
end



%% WATERFALL PLOT

% Concatenating EEG signals
data_tot = repmat(struct(), num_subjects, 1);

for subj = 1:num_subjects
    for deriv = 1:num_derivations
        data1_task = data_task(subj).(deriv_names{deriv});
        data1_rest = data_rest(subj).(deriv_names{deriv});
        
        data_tot(subj).(deriv_names{deriv}) = [data1_rest', data1_task'];
    end
end

t2 = 0:1/fs:120 - 1/fs;

%% Analysis Theta band 

window_length = fs * 3;  
shift = window_length * 0.5;  
nfft = 2048;  
noverlap = fs * 0.5;  

deriv_names2 = {'FZ', 'F3', 'F4'}; 

% Plot
figure(6);
sgtitle('Waterfall Plots for FZ, F3, F4 across Subjects') 
for subj = 1:num_subjects
    for rep = 1:length(deriv_names2)
        i = 1;
        counter = 1;
        CSA{rep} = [];
        time{rep} = [];
    
        % Length specific derivation
        N = length(data_tot(subj).(deriv_names2{rep}));
    
        while i + window_length < N 
            EEG_segment = detrend(data_tot(subj).(deriv_names2{rep})(i:i + window_length));
            [PSD{rep}, fw{rep}] = pwelch(EEG_segment, hamming(fs), noverlap, nfft, fs);
            PSD_theta{rep} = PSD{rep}(fw{rep} < 8 & fw{rep} > 4);
    
            CSA{rep}(counter, :) = PSD_theta{rep};
            time{rep}(counter) = t2(i + window_length / 2) / 60; % in minutes
    
            i = i + shift;
            counter = counter + 1;
        end
    
        % Plot
        freq = fw{rep}(fw{rep} < 8 & fw{rep} > 4); 
        [Time, Frequency] = meshgrid(time{rep}, freq);
        
        % Subplot for each subject and derivation
        subplot(num_subjects, length(deriv_names2), (subj - 1) * length(deriv_names2) + rep);
        waterfall(Time', Frequency', CSA{rep});
    
        xlabel('Time [m]');
        ylabel('Theta Band [Hz]');
        zlabel('\muV^{2}/Hz');
        title([deriv_names2{rep} ' - Subject ' num2str(k(subj))]);
    
        xline(t(end) / 60); % Showing end of rest condition
    end
end

%% Analysis Beta band 

deriv_names2 = {'O1', 'O2'};

figure;
sgtitle('Waterfall Plots for O1, O2 across Subjects') 
for subj = 1:num_subjects
    for rep = 1:length(deriv_names2)
        i = 1;
        counter = 1;
        CSA{rep} = [];
        time{rep} = [];
    
        % Length of specific derivation
        N = length(data_tot(subj).(deriv_names2{rep}));
    
        while i + window_length < N 
            EEG_segment = detrend(data_tot(subj).(deriv_names2{rep})(i:i + window_length));
            [PSD{rep}, fw{rep}] = pwelch(EEG_segment, hamming(fs), noverlap, nfft, fs);
            PSD_beta{rep} = PSD{rep}(fw{rep} < 30 & fw{rep} > 12);
    
            CSA{rep}(counter, :) = PSD_beta{rep};
            time{rep}(counter) = t2(i + window_length / 2) / 60; % in minutes
    
            i = i + shift;
            counter = counter + 1;
        end
    
        % Plot
        freq = fw{rep}(fw{rep} < 30 & fw{rep} > 12); 
        [Time, Frequency] = meshgrid(time{rep}, freq);
        
        % Subplot for each subject and derivation
        subplot(num_subjects, length(deriv_names2), (subj - 1) * length(deriv_names2) + rep);
        waterfall(Time', Frequency', CSA{rep});
    
        xlabel('Time [m]');
        ylabel('Beta Band [Hz]');
        zlabel('\muV^{2}/Hz');
        title([deriv_names2{rep} ' - Subject ' num2str(k(subj))]);
    
        xline(t(end) / 60); %  Showing end of rest condition
    end
end