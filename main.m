%% SET UP
close all
clear
clc

chanlocs = load('chanlocs.mat'); % Load channel location data
chanlocs = chanlocs.chanlocs;
num_channels = length(chanlocs);

bands = struct('theta1', [4.1, 5.8], 'theta2', [5.9, 7.4], 'beta1', [13, 19.9], 'beta2', [20, 25]);
bandNames = fieldnames(bands);

lowCut= 2;
highCut = 30;

fs = 500;

for k = [1, 2, 3, 5, 7, 8]
    % Load EEG data
    subject_rest = load(sprintf('Subject%02d_1.mat', k));
    subject_task = load(sprintf('Subject%02d_2.mat', k));
    
    % ---------- Selezionare la derivazione di interesse ---------- %
    deriv_name = 'F7'; % Modifica questo valore per selezionare un altro canale
    % ---------- Selezionare la derivazione di interesse ---------- %
    
    [EEG_rest, EEG_task] = plot_EEG_periodogram(subject_rest, subject_task, deriv_name, fs, k);

    % ANALIZZANDO I PERIODOGRAMMI NON è NECESSARIO APPLICARE FILTRO NOTCH
    % APPLICO FILTRO PASSA-BANDA PER STUDIARE LE FREQ DI INTERESSE e
    % NORMALIZZO (-mean /std)
    
    subject_rest_filtered = struct();
    subject_task_filtered = struct();
    fieldNames = fieldnames(subject_rest);
    for i = 1:length(fieldNames)
        deriv_name = fieldNames{i};      
        signal_rest = subject_rest.(deriv_name)(60001:90000);
        signal_task = subject_task.(deriv_name)(1:30000);
    
        filtered_signal_rest = process_sign(signal_rest, lowCut, highCut, fs);
        filtered_signal_task = process_sign(signal_task, lowCut, highCut, fs);
        
        subject_rest_filtered.(deriv_name) = filtered_signal_rest;
        subject_task_filtered.(deriv_name) = filtered_signal_task;
    end

    % potenze per tutte nelle 4 bande per i 21 canali
    power_values_rest = zeros(num_channels, length(bandNames));
    power_values_task = zeros(num_channels, length(bandNames));
    
    % Calcolo delle potenze per ogni canale e banda
    for i = 1:num_channels
        deriv_name = chanlocs(i).labels;
        if isfield(subject_rest_filtered, deriv_name)
            % Carico il segnale
            signal_rest = subject_rest_filtered.(deriv_name);
            signal_task = subject_task_filtered.(deriv_name);
            
            % Calcolo della PSD utilizzando il metodo di Welch
            [PSD_rest, f] = pwelch(signal_rest, [], [], [], fs);
            [PSD_task, ~] = pwelch(signal_task, [], [], [], fs);
            
            % Calcolo della potenza per ciascuna banda
            for j = 1:length(bandNames)
                band = bands.(bandNames{j});
                band_indices = f >= band(1) & f <= band(2); % indici degli elementi appartenenti alla j-esima banda
                
                % Verifica che ci siano dati nella banda
                if any(band_indices)
                    % Calcolo della potenza integrando la PSD nella banda
                    power_rest = trapz(f(band_indices), PSD_rest(band_indices));
                    power_task = trapz(f(band_indices), PSD_task(band_indices));
                else
                    power_rest = NaN;
                    power_task = NaN;
                end
                
                % Memorizzare le potenze nel vettore allineato ai canali
                power_values_rest(i, j) = power_rest;
                power_values_task(i, j) = power_task;
            end
        else
            % Canale non è presente -> NaN
            power_values_rest(i, :) = NaN;
            power_values_task(i, :) = NaN;
        end
    end
    
    % Creazione della figura per il subject k
    figure('Name', sprintf('Subject %d - Power Density Maps', k), 'NumberTitle', 'off');
    
    % Impostare il numero di righe e colonne per i subplot
    num_rows = 2; % Rest e Task
    num_cols = length(bandNames); % Numero di bande
    
    % Iterare su ciascuna banda per creare i subplot
    for j = 1:length(bandNames)
        band_label = bandNames{j};
        
        data_rest = power_values_rest(:, j);
        data_task = power_values_task(:, j);
        
        % Subplot rest period
        subplot(num_rows, num_cols, j);
        topoplot(data_rest, chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
        title(sprintf('Rest - %s', band_label));
        colorbar;
        
        % Subplot task period
        subplot(num_rows, num_cols, j + length(bandNames));
        topoplot(data_task, chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
        title(sprintf('Task - %s', band_label));
        colorbar;
    end
    
    mainTitle = sgtitle(sprintf('Subject %d - Power Density Maps', k));
    mainTitle.FontWeight = 'bold';

end

%%


% fs = 500;  % Frequenza di campionamento
window_length = fs * 1;       % Lunghezza finestra di 1 secondo
shift = window_length * 0.5;  % Sovrapposizione del 50%
order = 15;                   % Ordine del modello
nfft = 2048;                  % Numero di punti FFT
%lowCut = 2;
%highCut = 30;

subjects = [1, 2, 3, 5, 7, 8];

% ---------- Selezionare la derivazione di interesse ---------- %
derivazione_interesse = 'C3';  % Modifica questo valore per selezionare un altro canale
% ---------- Selezionare la derivazione di interesse ---------- %

% Creare la figura per tutti i soggetti
figure('Name', sprintf('Waterfall Plot - Derivazione %s', derivazione_interesse), 'NumberTitle', 'off','Position', [100, 100, 1200, 800]);

% Calcolare il numero di subplot necessari
num_subjects = length(subjects);
num_rows = 2;
num_cols = num_subjects;

% Contatore per i subplot
subplot_idx = 1;

for idx = 1:num_subjects
    k = subjects(idx);
    
    % Caricare i dati EEG per rest e task
    subject_rest = load(sprintf('Subject%02d_1.mat', k));
    subject_task = load(sprintf('Subject%02d_2.mat', k));
    

    if isfield(subject_rest, derivazione_interesse) && isfield(subject_task, derivazione_interesse)
        % Estrarre il segnale per la derivazione di interesse
        EEG_rest = subject_rest.(derivazione_interesse)(60001:90000);
        EEG_task = subject_task.(derivazione_interesse)(1:30000);
        
        % Filtrare i segnali
        EEG_rest = process_sign(EEG_rest, lowCut, highCut, fs);
        EEG_task = process_sign(EEG_task, lowCut, highCut, fs);
              
        % Calcolare CSA per rest
        [CSA_rest, time_rest, freq_rest] = compute_CSA(EEG_rest, fs, window_length, shift, order, nfft);
        
        % Calcolare CSA per task
        [CSA_task, time_task, freq_task] = compute_CSA(EEG_task, fs, window_length, shift, order, nfft);
        
        % Plot per la condizione di riposo (rest)
        subplot(num_rows, num_cols, idx);
        [Time_rest, Frequency_rest] = meshgrid(time_rest / 60, freq_rest);
        waterfall(Time_rest', Frequency_rest', CSA_rest);
        xlabel('Time [m]');
        ylabel('Frequency [Hz]');
        zlabel('PSD [\muV^{2}/Hz]');
        title(sprintf("Rest - Subject %d - Derivation 'C3'", k));
        grid on;
        hold on;
        
        
        % Plot per la condizione di task
        subplot(num_rows, num_cols, idx + num_subjects);
        [Time_task, Frequency_task] = meshgrid(time_task / 60, freq_task);
        waterfall(Time_task', Frequency_task', CSA_task);
        xlabel('Time [m]');
        ylabel('Frequency [Hz]');
        zlabel('PSD [\muV^{2}/Hz]');
        title(sprintf('Task - S%d', k));
        grid on;
        hold on;
       
    else
        fprintf('Derivazione %s non trovata per il soggetto %d.\n', derivazione_interesse, k);
    end
end


%%

function [EEG_rest, EEG_task] = plot_EEG_periodogram(subj_rest, subj_task, deriv_name, fs, k)
    EEG_rest = subj_rest.(deriv_name)(60001:90000); % analizzo i primi 30s
    EEG_task = subj_task.(deriv_name)(1:30000); %analizzo gli ultimi 30s
    
    t = (0:length(EEG_rest)-1) / fs;
    
    % Periodogramma
    N = length(EEG_rest);
    
    [PSDx_rest, freq] = periodogram(EEG_rest, rectwin(N), N, fs);
    [PSDx_task, ~] = periodogram(EEG_task, rectwin(N), N, fs);


    figure('Name', sprintf(['EEG - Soggetto ' num2str(k), ', Derivazione ' deriv_name]), 'NumberTitle', 'off');
    % Primo sottoplot: segnale nel dominio del tempo (fase di riposo e task)
    subplot(1, 2, 1);
    plot(t, EEG_rest, 'b');
    xlabel('[s]');
    ylabel('[\mu V]');
    legend('Riposo');

    subplot(1, 2, 2);
    plot(t, EEG_task, 'r');
    xlabel('[s]');
    ylabel('[\mu V]');
    mainTitle = sgtitle(['Segnale EEG - Soggetto ' num2str(k) ', Derivazione ' deriv_name]);
    mainTitle.FontWeight = 'bold';
    legend('Task');
    
    figure('Name', sprintf(['Periodogramma - Soggetto ' num2str(k), ', Derivazione ' deriv_name]), 'NumberTitle', 'off');
    % Secondo sottoplot: periodogramma (fase di riposo e task)
    subplot(1, 2, 1);
    plot(freq, PSDx_rest, 'b');
    xlabel('[Hz]');
    ylabel('[\mu V^{2}/Hz]');
    legend('Riposo');

    subplot(1, 2, 2);
    plot(freq, PSDx_task, 'r');
    xlabel('[Hz]');
    ylabel('[\mu V^{2}/Hz]');
    mainTitle = sgtitle(['Periodogramma - Soggetto ' num2str(k) ', Derivazione ' deriv_name]);
    mainTitle.FontWeight = 'bold';
    legend('Task');

end

function eeg_filt = process_sign(eeg, lowCut, highCut, fs)
    % Bandpass Filter
    [b_bp, a_bp] = butter(4, [lowCut highCut] / (fs / 2));
    eeg_filt = filtfilt(b_bp, a_bp, eeg);

    % Normalization
    eeg_filt = (eeg_filt - mean(eeg_filt)) / std(eeg_filt);
end

% Funzione per calcolare la CSA
function [CSA, time, freq] = compute_CSA(EEG_signal, fs, window_length, shift, order, nfft)
    N = length(EEG_signal);
    i = 1;
    counter = 1;
    CSA = [];
    time = [];
    while i + window_length <= N
        % Selezionare il segmento EEG
        EEG_segment = detrend(EEG_signal(i:i+window_length-1));        
                
        % Stimare la PSD
        [PSD, f] = pyulear(EEG_segment, order, nfft, fs);
        
        % Selezionare solo fino a 30 Hz
        PSD = PSD(f < 30);
        f_segment = f(f < 30);

        % Accumulare i dati
        CSA(counter, :) = PSD';
        time(counter) = (i + window_length / 2) / fs;  % Tempo in secondi

        i = i + shift;
        counter = counter + 1;
    end
    freq = f_segment;
end