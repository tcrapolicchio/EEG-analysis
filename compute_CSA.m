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