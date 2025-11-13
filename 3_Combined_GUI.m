function CombinedAudioGUIfinal_11
    fig = figure('Name', 'Combined Audio Processing GUI', ...
                 'NumberTitle', 'on', ...
                 'Position', [50, 50, 1400, 750], ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'none', ...
                 'Resize', 'on');

    tabgroup = uitabgroup(fig);
    tab1 = uitab(tabgroup, 'Title', 'Audio File Mode');
    tab2 = uitab(tabgroup, 'Title', 'Voice Recorder Mode');
    tab3 = uitab(tabgroup, 'Title', 'Audio Watermarking');
   
    
    panel1 = uipanel('Parent', tab1, ...
                     'Units', 'normalized', ...
                     'Position', [0 0 1 1], ...
                     'BackgroundColor', [1 0.5 0]);

    panel2 = uipanel('Parent', tab2, ...
                     'Units', 'normalized', ...
                     'Position', [0 0 1 1], ...
                     'BackgroundColor', [0 0.447 0.741]);
                  
    panel3 = uipanel('Parent', tab3, ...
                     'Units', 'normalized', ...
                     'Position', [0 0 1 1], ...
                     'BackgroundColor', [0 0.447 0.741]);


    AudioEffectAndNoiseRemovalGUI(tab1);
    VoiceRecorderGUI(tab2);
    AudiowatermarkingGUI(tab3);
end

function AudioEffectAndNoiseRemovalGUI(parent)
    fig = uipanel(parent, 'Units', 'normalized', 'Position', [0 0 1 1],'BackgroundColor', [0 0 0]);

    % File Selection
    uicontrol(fig, 'Style', 'text', 'String', 'File Name', 'Position', [500, 680, 60, 25], 'HorizontalAlignment', 'left','BackgroundColor', [0 0 0], 'ForegroundColor', [1 1 1]);
    fileNameInput = uicontrol(fig, 'Style', 'edit', 'Position', [580, 680, 300, 25],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Browse', 'Position', [900, 680, 90, 25],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1], 'Callback', @browseFile);

    % Effects Section
    effects = {'Gain', 'Trim', 'Fade In', 'Fade Out', 'Reverse', 'Speed Up/Down', 'Normalize', 'Compressor', 'Distortion', 'Echo', 'Reverb', 'Delay', 'Chorus', 'Equalizer', 'Pitch Shifting'};
    numEffects = length(effects);
    numColumns = 5;
    buttonWidth = 65;
    buttonHeight = 30;
    spacingX = 10;
    spacingY = 6;
    startX = 50;
    startY = 650;

    for i = 1:numEffects
        col = mod(i-1, numColumns);
        row = floor((i-1) / numColumns);
        xPos = startX + col * (buttonWidth + spacingX);
        yPos = startY - row * (buttonHeight + spacingY);
        uicontrol(fig, 'Style', 'pushbutton', 'String', effects{i}, 'Position', [xPos, yPos, buttonWidth, buttonHeight],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1], 'Callback', {@applyEffect, effects{i}});
    end

    % Sliders for effect parameters
    uicontrol(fig, 'Style', 'text', 'String', 'Gain', 'Position', [50, 400, 60, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    gainSlider = uicontrol(fig, 'Style', 'slider', 'Min', 0, 'Max', 3, 'Value', 1, 'Position', [120, 400, 200, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    uicontrol(fig, 'Style', 'text', 'String', 'Fade In Duration (s)', 'Position', [50, 360, 120, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    fadeInSlider = uicontrol(fig, 'Style', 'slider', 'Min', 0, 'Max', 10, 'Value', 3, 'Position', [180, 360, 200, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    uicontrol(fig, 'Style', 'text', 'String', 'Fade Out Duration (s)', 'Position', [50, 320, 120, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    fadeOutSlider = uicontrol(fig, 'Style', 'slider', 'Min', 0, 'Max', 10, 'Value', 3, 'Position', [180, 320, 200, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    uicontrol(fig, 'Style', 'text', 'String', 'Speed Factor', 'Position', [50, 280, 120, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    speedSlider = uicontrol(fig, 'Style', 'slider', 'Min', 0.5, 'Max', 2, 'Value', 1, 'Position', [180, 280, 200, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    uicontrol(fig, 'Style', 'text', 'String', 'Pitch Shift (semitones)', 'Position', [50, 240, 150, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    pitchSlider = uicontrol(fig, 'Style', 'slider', 'Min', -12, 'Max', 12, 'Value', 0, 'Position', [180, 240, 200, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    uicontrol(fig, 'Style', 'text', 'String', 'Trim Duration (s)', 'Position', [50, 200, 120, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    trimSlider = uicontrol(fig, 'Style', 'slider', 'Min', 0, 'Max', 10, 'Value', 5, 'Position', [180, 200, 200, 20],'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    % Waveform Displays
    axesInput = axes('Parent', fig, 'Units', 'pixels', 'Position', [550, 550, 600, 100],'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1]);
    title(axesInput, 'Input Signal', 'Color', [1 1 1]);
    xlabel(axesInput, 'Time (s)', 'Color', [1 1 1]);
    ylabel(axesInput, 'Amplitude', 'Color', [1 1 1]);

    axesOutput = axes('Parent', fig, 'Units', 'pixels', 'Position', [550, 400, 600, 100],'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1]);
    title(axesOutput, 'Output Signal', 'Color', [1 1 1]);
    xlabel(axesOutput, 'Time (s)', 'Color', [1 1 1]);
    ylabel(axesOutput, 'Amplitude', 'Color', [1 1 1]);

    axesNoisy = axes('Parent', fig, 'Units', 'pixels', 'Position', [550, 250, 600, 100],'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1]);
    title(axesNoisy, 'Noisy Signal', 'Color', [1 1 1]);
    xlabel(axesNoisy, 'Time (s)', 'Color', [1 1 1]);
    ylabel(axesNoisy, 'Amplitude', 'Color', [1 1 1]);

    axesFiltered = axes('Parent', fig, 'Units', 'pixels', 'Position', [550, 100, 600, 100],'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1]);
    title(axesFiltered, 'Filtered Signal', 'Color', [1 1 1]);
    xlabel(axesFiltered, 'Time (s)', 'Color', [1 1 1]);
    ylabel(axesFiltered, 'Amplitude', 'Color', [1 1 1]);

    % Control Buttons
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Show Original Graph', 'Position', [50, 520, 120, 30], 'Callback', @showOriginalGraph,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play Original Audio', 'Position', [200, 520, 120, 30], 'Callback', @playOriginalAudio,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Pause Audio', 'Position', [50, 480, 120, 30], 'Callback', @pauseAudio,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play Processed Audio', 'Position', [200, 480, 120, 30], 'Callback', @playProcessedAudio,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Reset', 'Position', [500, 20, 80, 30], 'Callback', @resetAudio,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save Processed Audio', 'Position', [600, 20, 150, 30], 'Callback', @saveProcessedAudio,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Add Noise', 'Position', [50, 150, 120, 30], 'Callback', @selectNoiseType,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Apply Filter', 'Position', [200, 150, 120, 30], 'Callback', @selectFilterType,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play Noisy', 'Position', [50, 100, 120, 30], 'Callback', @playNoisy,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play Filtered', 'Position', [200, 100, 120, 30], 'Callback', @playFiltered,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop Noise Sig.', 'Position', [50, 50, 120, 30], 'Callback', @stopAudio,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save Filtered', 'Position', [200, 50, 120, 30], 'Callback', @saveFiltered,'BackgroundColor', [1 0.5 0], 'ForegroundColor', [1 1 1]);

    % Variables to store audio data and playback state
    audioData = struct('y', [], 'fs', [], 'processed', [], 'noisy', [], 'filtered', []);
    audioState = struct('player', [], 'paused', false);

    % Callback Functions
    function browseFile(~, ~)
        [fileName, filePath] = uigetfile({'*.wav;*.mp3', 'Audio Files (*.wav, *.mp3)'});
        if fileName
            fileFullPath = fullfile(filePath, fileName);
            set(fileNameInput, 'String', fileFullPath);
            [audioData.y, audioData.fs] = audioread(fileFullPath);
            audioData.processed = audioData.y;
            audioData.noisy = audioData.y;
            audioData.filtered = audioData.y;
            plotWaveform(audioData.y, audioData.fs, axesInput, 'b');
        end
    end

    function showOriginalGraph(~, ~)
        if ~isempty(audioData.y)
            plotWaveform(audioData.y, audioData.fs, axesInput, 'b');
        else
            errordlg('Please upload an audio file first.');
        end
    end

    function playOriginalAudio(~, ~)
        if ~isempty(audioData.y)
            stopAudio();
            audioState.player = audioplayer(audioData.y, audioData.fs);
            audioState.paused = false;
            play(audioState.player);
        else
            errordlg('Please upload an audio file first.');
        end
    end

    function playProcessedAudio(~, ~)
        if ~isempty(audioData.processed)
            stopAudio();
            audioState.player = audioplayer(audioData.processed, audioData.fs);
            audioState.paused = false;
            play(audioState.player);
        else
            errordlg('No processed audio to play.');
        end
    end

    function pauseAudio(~, ~)
        if ~isempty(audioState.player)
            if isplaying(audioState.player)
                pause(audioState.player);
                audioState.paused = true;
            elseif audioState.paused
                resume(audioState.player);
                audioState.paused = false;
            end
        else
            errordlg('No audio is currently playing.');
        end
    end

    function resetAudio(~, ~)
        if ~isempty(audioData.y)
            audioData.processed = audioData.y;
            audioData.noisy = audioData.y;
            audioData.filtered = audioData.y;
            cla(axesInput); cla(axesOutput); cla(axesNoisy); cla(axesFiltered); % Clear all axes
            drawnow; % Force immediate update of the figure
            plotWaveform(audioData.y, audioData.fs, axesInput, 'b');
            plotWaveform(audioData.processed, audioData.fs, axesOutput, 'g');
            plotWaveform(audioData.noisy, audioData.fs, axesNoisy, 'r');
            plotWaveform(audioData.filtered, audioData.fs, axesFiltered, 'm');
        else
            errordlg('Please upload an audio file first.');
        end
    end

    function saveProcessedAudio(~, ~)
        if ~isempty(audioData.processed)
            maxAmplitude = max(abs(audioData.processed));
            if maxAmplitude > 1
                audioData.processed = audioData.processed / maxAmplitude;
            end
            [fileName, filePath] = uiputfile('*.wav', 'Save Processed Audio');
            if fileName
                audiowrite(fullfile(filePath, fileName), audioData.processed, audioData.fs);
            end
        else
            errordlg('No processed audio to save.');
        end
    end

    function stopAudio(~, ~)
        if ~isempty(audioState.player)
            stop(audioState.player);
            audioState.paused = false;
        end
    end

    function plotWaveform(audio, fs, ax, color)
        t = (0:length(audio)-1) / fs;
        axes(ax);
        plot(t, audio, 'Color', color, 'LineWidth', 1.5);
        grid on;
        xlim([0 t(end)]);
        ylim([-1 1]);
    end

    function applyEffect(~, ~, effectName)
        if ~isempty(audioData.y) %is it loading
            gain = get(gainSlider, 'Value'); %sliders
            fadeInDuration = get(fadeInSlider, 'Value');
            fadeOutDuration = get(fadeOutSlider, 'Value');
            speedFactor = get(speedSlider, 'Value');
            pitchShiftAmount = get(pitchSlider, 'Value');
            trimDuration = get(trimSlider, 'Value');

            switch effectName
                case 'Gain'
                    audioData.processed = audioData.y * gain;
                case 'Trim'
                    trimSamples = round(trimDuration * audioData.fs);
                    audioData.processed = audioData.y(trimSamples:end, :);
                case 'Fade In'
                    fadeSamples = min(length(audioData.y), fadeInDuration * audioData.fs);
                    fadeCurve = linspace(0, 1, fadeSamples)';
                    audioData.processed = audioData.y;
                    for channel = 1:size(audioData.y, 2)
                        audioData.processed(1:fadeSamples, channel) = audioData.y(1:fadeSamples, channel) .* fadeCurve;
                    end
                case 'Fade Out'
                    fadeSamples = min(length(audioData.y), fadeOutDuration * audioData.fs);
                    fadeCurve = linspace(1, 0, fadeSamples)';
                    audioData.processed = audioData.y;
                    for channel = 1:size(audioData.y, 2)
                        audioData.processed(end-fadeSamples+1:end, channel) = audioData.y(end-fadeSamples+1:end, channel) .* fadeCurve;
                    end
                case 'Speed Up/Down'
                    audioData.processed = resample(audioData.y, round(speedFactor * audioData.fs), audioData.fs);
                case 'Reverse'
                    audioData.processed = flipud(audioData.y);
                case 'Normalize'
                    maxAmplitude = max(abs(audioData.y));
                    audioData.processed = audioData.y / maxAmplitude;
                case 'Compressor'
                    threshold = 0.5;
                    ratio = 4;
                    audioData.processed = min(audioData.y, threshold + (audioData.y - threshold) / ratio);
                case 'Distortion'
                    gain = 2;
                    audioData.processed = max(-1, min(1, gain * audioData.y));
                case 'Echo'
                    delay = 0.3               
                    ;
                    gain = 0.5;
                    audioData.processed = echoEffect(audioData.y, audioData.fs, delay, gain);
                case 'Reverb'
                    reverbTime = 0.5;
                    audioData.processed = reverbEffect(audioData.y, audioData.fs, reverbTime);
                case 'Delay'
                    delay = 0.3;
                    gain = 0.5;
                    audioData.processed = echoEffect(audioData.y, audioData.fs, delay, gain);
                case 'Chorus'
                    modulationDepth = 0.005;
                    modulationRate = 1.5;
                    audioData.processed = chorusEffect(audioData.y, audioData.fs, modulationDepth, modulationRate);
                case 'Equalizer'
                    audioData.processed = equalizerEffect(audioData.y, audioData.fs);
                case 'Pitch Shifting'
                    audioData.processed = pitchShiftEffect(audioData.y, audioData.fs, pitchShiftAmount);
                otherwise
                    warndlg('Effect not implemented.');
            end
            plotWaveform(audioData.processed, audioData.fs, axesOutput, 'g');
        else
            errordlg('Please upload an audio file first.'); %file isnt loaded
        end
    end

    function yOut = pitchShiftEffect(y, fs, semitones)
        factor = 2^(semitones / 12);
        yOut = resample(y, round(factor * fs), fs);
    end

    function yOut = echoEffect(y, fs, delay, gain)
        delaySamples = round(delay * fs);
        yOut = y;
        yOut(delaySamples+1:end, :) = yOut(delaySamples+1:end, :) + gain * y(1:end-delaySamples, :);
    end

    function yOut = reverbEffect(y, fs, reverbTime)
        alpha = 0.5;
        delaySamples = round(reverbTime * fs);
        yOut = y;
        for i = delaySamples+1:length(y)
            yOut(i) = yOut(i) + alpha * y(i - delaySamples);
        end
    end

    function yOut = equalizerEffect(y, fs)
        yOut = y;
    end

    function yOut = chorusEffect(y, fs, depth, rate)
        t = (0:length(y)-1) / fs;
        modSignal = depth * sin(2 * pi * rate * t);
        yOut = y;
        for i = 1:length(y)
            delaySamples = round(modSignal(i) * fs);
            if i > delaySamples
                yOut(i) = yOut(i) + y(i - delaySamples);
            end
        end
    end

    function selectNoiseType(~, ~)
        if ~isempty(audioData.y)
            noiseOptions = {'White Gaussian Noise', 'Salt & Pepper Noise', 'Spike Impulse Noise'};
            [indx, tf] = listdlg('PromptString', 'Select a noise type:', 'SelectionMode', 'single', 'ListString', noiseOptions);
            if tf == 1
                addNoise(noiseOptions{indx});
            end
        else
            errordlg('Please upload an audio file first.');
        end
    end

    function addNoise(selectedNoise)
        switch selectedNoise
            case 'White Gaussian Noise'
                noise = 0.02 * randn(size(audioData.y));
                audioData.noisy = audioData.y + noise;
            case 'Salt & Pepper Noise'
                audioData.noisy = imnoise(audioData.y, 'salt & pepper', 0.02);
            case 'Spike Impulse Noise'
                noise = zeros(size(audioData.y));
                numSpikes = max(5, round(length(audioData.y) / 1000)); % Approx 5-10 spikes
                spikeIndices = randi([1 length(audioData.y)], 1, numSpikes);
                noise(spikeIndices) = 0.5 * (2 * rand(1, numSpikes) - 1); % Random amplitude between -0.5 and 0.5
                audioData.noisy = audioData.y + noise;
            otherwise
                audioData.noisy = audioData.y;
        end
        plotWaveform(audioData.noisy, audioData.fs, axesNoisy, 'r');
    end

    function selectFilterType(~, ~)
        if ~isempty(audioData.noisy)
            filterOptions = {'Median Filter', 'Adaptive Filter (LMS)', 'Wiener Filter'};
            [indx, tf] = listdlg('PromptString', 'Select a filter type:', 'SelectionMode', 'single', 'ListString', filterOptions);
            if tf == 1
                applyFilter(filterOptions{indx});
            end
        else
            errordlg('Please add noise first.');
        end
    end

    function applyFilter(selectedFilter)
        if isempty(audioData.noisy)
            errordlg('Add noise first before filtering.');
            return;
        end
        switch selectedFilter
            case 'Median Filter'
                y = audioData.noisy;
                for ch = 1:size(y, 2)
                    y(:, ch) = medfilt1(y(:, ch), 5); % Window size of 5, adjustable
                end
                audioData.filtered = y;
                plotWaveform(y, audioData.fs, axesFiltered, 'm');
            case 'Adaptive Filter (LMS)'
                y = audioData.noisy;
                mu = 0.01; % Step size, adjustable
                order = 10; % Filter order, adjustable
                filtered = zeros(size(y));
                for ch = 1:size(y, 2)
                    w = zeros(order, 1); % Initial weights
                    for n = order:length(y)
                        x = y(max(1, n-order+1):n, ch); % Input vector
                        e = y(n, ch) - w' * x; % Error
                        w = w + mu * e * x; % Weight update
                        filtered(n, ch) = e; % Output is error (noise estimate subtracted)
                    end
                end
                audioData.filtered = filtered;
                plotWaveform(filtered, audioData.fs, axesFiltered, 'm');
            case 'Wiener Filter'
                y = audioData.noisy;
                filtered = zeros(size(y));
                for ch = 1:size(y, 2)
                    filtered(:, ch) = wienerFilter(y(:, ch), audioData.fs);
                end
                audioData.filtered = filtered;
                plotWaveform(filtered, audioData.fs, axesFiltered, 'm');
        end
    end

    function yOut = wienerFilter(y, fs)
        % Wiener Filter implementation for 1D audio signal
        % Estimate noise power spectrum using a segment of the signal
        noiseSegmentLength = round(0.1 * fs); % Use first 0.1 seconds for noise estimation
        if length(y) < noiseSegmentLength
            noiseSegmentLength = length(y);
        end
        noise = y(1:noiseSegmentLength); % Assume initial segment is noise
        noisePower = mean(noise.^2); % Noise power estimation
        
        % Signal power estimation
        signalPower = mean(y.^2);
        
        % Wiener filter gain
        H = max(0, (signalPower - noisePower) / signalPower);
        
        % Apply Wiener filter in frequency domain
        Y = fft(y);
        Y_filtered = Y * H;
        yOut = real(ifft(Y_filtered));
        
        % Ensure output length matches input
        yOut = yOut(1:length(y));
    end

    function playNoisy(~, ~)
        if ~isempty(audioData.noisy)
            stopAudio();
            audioState.player = audioplayer(audioData.noisy, audioData.fs);
            audioState.paused = false;
            play(audioState.player);
        else
            errordlg('No noisy audio to play.');
        end
    end

    function playFiltered(~, ~)
        if ~isempty(audioData.filtered)
            stopAudio();
            audioState.player = audioplayer(audioData.filtered, audioData.fs);
            audioState.paused = false;
            play(audioState.player);
        else
            errordlg('No filtered audio to play.');
        end
    end

    function saveFiltered(~, ~)
        if ~isempty(audioData.filtered)
            maxAmplitude = max(abs(audioData.filtered));
            if maxAmplitude > 1
                audioData.filtered = audioData.filtered / maxAmplitude;
            end
            [fileName, filePath] = uiputfile('*.wav', 'Save Filtered Audio');
            if fileName
                audiowrite(fullfile(filePath, fileName), audioData.filtered, audioData.fs);
            end
        else
            errordlg('No filtered audio to save.');
        end
    end

    function yOut = spectralSubtract(y, noise)
        y = y(:);
        noise = noise(:);
        N = length(y);
        Y = fft(y);
        Nfft = fft(noise, N);
        magY = abs(Y);
        magN = abs(Nfft);
        phaseY = angle(Y);
        cleanMag = max(magY - magN, 0);
        Yclean = cleanMag .* exp(1i * phaseY);
        yOut = real(ifft(Yclean));
    end
end

function VoiceRecorderGUI(parent)
    fig = uipanel(parent, 'Units', 'normalized', 'Position', [0 0 1 1],'BackgroundColor', [0 0 0]);

    global player Fs recorded_audio noisy_audio filtered_audio;
    Fs = 16000; % Adjusted to 16 kHz for voice
    duration = 5;
    bgBlack = [0 0 0];
    textWhite = [1 1 1];
    boxOrange = [1 0.5 0];

    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Set Duration', ...
        'Position', [60, 600, 180, 30], 'Callback', @setDuration ,'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Record Audio', ...
        'Position', [60, 560, 180, 30], 'Callback', @recordAudio, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play Original', ...
        'Position', [60, 520, 180, 30], 'Callback', @playOriginal, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    noisePopup = uicontrol(fig, 'Style', 'popupmenu', 'String', {'Select Noise', 'Spike Impulse Noise', 'White Gaussian Noise (WGN)'}, ...
        'Position', [60, 480, 180, 30], 'Callback', @addNoise, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    filterPopup = uicontrol(fig, 'Style', 'popupmenu', 'String', {'Low Pass Filter', 'High Pass Filter', 'Median Filter', 'Adaptive Filter (LMS)', 'Wiener Filter'}, ...
        'Position', [60, 440, 180, 30], 'Callback', @updateFilterType, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Set Cutoff', ...
        'Position', [60, 400, 180, 30], 'Callback', @setCutoff, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Apply Filter', ...
        'Position', [60, 360, 180, 30], 'Callback', @applyFilter, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play Noisy', ...
        'Position', [60, 320, 180, 30], 'Callback', @playNoisy, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play Filtered', ...
        'Position', [60, 280, 180, 30], 'Callback', @playFiltered, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save WAV', ...
        'Position', [60, 240, 180, 30], 'Callback', @saveAudio, 'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);

    % Waveform and Spectrogram Displays (3 rows, 2 columns with gap)
    axes1_time = axes('Parent', fig, 'Units', 'pixels', 'Position', [300, 500, 450, 100], 'Color', bgBlack, ...
        'XColor', textWhite, 'YColor', textWhite);
    title(axes1_time, 'Original Signal (Amplitude vs Samples)', 'Color', textWhite);
    xlabel(axes1_time, 'Samples', 'Color', textWhite);
    ylabel(axes1_time, 'Amplitude', 'Color', textWhite);

    axes1_spec = axes('Parent', fig, 'Units', 'pixels', 'Position', [800, 500, 450, 100], 'Color', bgBlack, ...
        'XColor', textWhite, 'YColor', textWhite);
    title(axes1_spec, 'Original Signal (Frequency vs Samples)', 'Color', textWhite);
    xlabel(axes1_spec, 'Samples', 'Color', textWhite);
    ylabel(axes1_spec, 'Frequency (Hz)', 'Color', textWhite);

    axes2_time = axes('Parent', fig, 'Units', 'pixels', 'Position', [300, 300, 450, 100], 'Color', bgBlack, ...
        'XColor', textWhite, 'YColor', textWhite);
    title(axes2_time, 'Noisy Signal (Amplitude vs Samples)', 'Color', textWhite);
    xlabel(axes2_time, 'Samples', 'Color', textWhite);
    ylabel(axes2_time, 'Amplitude', 'Color', textWhite);

    axes2_spec = axes('Parent', fig, 'Units', 'pixels', 'Position', [800, 300, 450, 100], 'Color', bgBlack, ...
        'XColor', textWhite, 'YColor', textWhite);
    title(axes2_spec, 'Noisy Signal (Frequency vs Samples)', 'Color', textWhite);
    xlabel(axes2_spec, 'Samples', 'Color', textWhite);
    ylabel(axes2_spec, 'Frequency (Hz)', 'Color', textWhite);

    axes3_time = axes('Parent', fig, 'Units', 'pixels', 'Position', [300, 100, 450, 100], 'Color', bgBlack, ...
        'XColor', textWhite, 'YColor', textWhite);
    title(axes3_time, 'Filtered Signal (Amplitude vs Samples)', 'Color', textWhite);
    xlabel(axes3_time, 'Samples', 'Color', textWhite);
    ylabel(axes3_time, 'Amplitude', 'Color', textWhite);

    axes3_spec = axes('Parent', fig, 'Units', 'pixels', 'Position', [800, 100, 450, 100], 'Color', bgBlack, ...
        'XColor', textWhite, 'YColor', textWhite);
    title(axes3_spec, 'Filtered Signal (Frequency vs Samples)', 'Color', textWhite);
    xlabel(axes3_spec, 'Samples', 'Color', textWhite);
    ylabel(axes3_spec, 'Frequency (Hz)', 'Color', textWhite);

    % Stop and Reset Buttons
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop', ...
        'Position', [60, 60, 180, 30], 'Callback', @stopPlayback,'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Reset', ...
        'Position', [60, 20, 180, 30], 'Callback', @reset,'BackgroundColor', boxOrange, 'ForegroundColor', textWhite);

    filtered_audio = [];
    filterType = 'Low Pass Filter';
    cutoffFreq = 1000; % Default cutoff

    function setDuration(~, ~)
        newDuration = inputdlg('Enter recording duration (seconds):', 'Set Duration', [1 30], {num2str(duration)});
        if ~isempty(newDuration)
            duration = str2double(newDuration{1});
            if isnan(duration) || duration <= 0
                duration = 5;
                warndlg('Invalid duration. Set to default 5 seconds.');
            end
        end
    end

    function recordAudio(~, ~)
        recObj = audiorecorder(Fs, 16, 1);% sampling freq,16 bit depth,1 input microphone
        recordblocking(recObj, duration);%recobj= audio record
        recorded_audio = getaudiodata(recObj);
        samples = 0:length(recorded_audio)-1;
        axes(axes1_time);
        plot(samples, recorded_audio, 'b', 'LineWidth', 1.5);
        xlim([0 length(recorded_audio)-1]);
        title(axes1_time, 'Original Recorded Signal (Amplitude vs Samples)', 'Color', textWhite);
        [s, f, t_spec] = spectrogram(recorded_audio, 256, 200, 256, Fs, 'yaxis');% FFT size.overlap.window size,samp freq,freq in y axis
        axes(axes1_spec);
        imagesc(1:size(s,2), f, 10*log10(abs(s)));%t vs f ke image akare dekhay
        axis xy;
        colormap jet;%jjet clr map
        xlim([1 size(s,2)]);
        colorbar;%clr scale bar
        title(axes1_spec, 'Original Signal (Frequency vs Samples)', 'Color', textWhite);
        noisy_audio = recorded_audio; % Initialize noisy_audio with original
        filtered_audio = [];
        cla(axes2_time);
        cla(axes2_spec);%cla-age graph muche ney
        cla(axes3_time);%noisy sig er time domain graph clear
        cla(axes3_spec);%noisy sig er spectrum clear
    end

    function playOriginal(~, ~)%callback function first one ignores object 2nd one ignores data
        if isempty(recorded_audio)%checking is the file empty
            errordlg('No original audio to play.');%will show error dialog box
            return;
        end
        stopPlayback();%previous playing audio will stop
        player = audioplayer(recorded_audio, Fs);
        play(player);
    end

    function addNoise(hObj, ~)%hObj= eta popup menu jeta noise select korar jonno user use kre
        if isempty(recorded_audio)
            errordlg('Please record audio first.'); return;
        end
        noiseType = hObj.String{hObj.Value};%hObj.string = show all option list
        %hObj.value= User kon option beche niyeche tar index
        if strcmp(noiseType, 'Select Noise')
            return;
        end
        switch noiseType
            case 'Spike Impulse Noise'
                noise = zeros(size(recorded_audio));
                numSpikes = max(5, round(length(recorded_audio) / 1000)); % Approx 5-10 spikes
                spikeIndices = randi([1 length(recorded_audio)], 1, numSpikes);%Random index selecr kora
                noise(spikeIndices) = 0.5 * (2 * rand(1, numSpikes) - 1); % Random amplitude between -0.5 and 0.5
                noisy_audio = recorded_audio + noise;
                samples = 0:length(noisy_audio)-1;%sample index making
                axes(axes2_time);
                plot(samples, noisy_audio, 'r', 'LineWidth', 1.5);% making noisy signal red
                xlim([0 length(noisy_audio)-1]);
                title(axes2_time, 'Noisy Signal with Spike Impulse Noise (Amplitude vs Samples)', 'Color', textWhite);
                [s, f, t_spec] = spectrogram(noisy_audio, 256, 200, 256, Fs, 'yaxis');
                axes(axes2_spec);
                imagesc(1:size(s,2), f, 10*log10(abs(s)));
                axis xy;
                colormap jet;
                xlim([1 size(s,2)]);
                colorbar;
                title(axes2_spec, 'Noisy Signal (Frequency vs Samples)');
            case 'White Gaussian Noise (WGN)'
                noise = 0.02 * randn(size(recorded_audio));%making gaussian noise 
                noisy_audio = recorded_audio + noise;
                samples = 0:length(noisy_audio)-1;
                axes(axes2_time);
                plot(samples, noisy_audio, 'r', 'LineWidth', 1.5);
                xlim([0 length(noisy_audio)-1]);
                title(axes2_time, 'Noisy Signal with White Gaussian Noise (Amplitude vs Samples)', 'Color', textWhite);
                [s, f, ~] = spectrogram(noisy_audio, 256, 200, 256, Fs, 'yaxis');
                axes(axes2_spec);
                imagesc(1:size(s,2), f, 10*log10(abs(s)));
                axis xy;
                colormap jet;
                xlim([1 size(s,2)]);
                colorbar;
                title(axes2_spec, 'Noisy Signal (Frequency vs Samples)');
        end
        filtered_audio = []; % Clear filtered audio when noise is added
        cla(axes3_time); cla(axes3_spec);
    end

    function updateFilterType(hObj, ~)
        filterType = hObj.String{hObj.Value};
    end

    function setCutoff(~, ~)
        freq = inputdlg('Enter cutoff frequency (Hz):', 'Set Cutoff', [1 30], {num2str(cutoffFreq)});
        %Ekta dialog box open kora hobe jekhane user cutoff freq input dibe
        if ~isempty(freq)
            cutoffFreq = str2double(freq{1});%input string to numeric convert
            if isnan(cutoffFreq) || cutoffFreq <= 0 || cutoffFreq > Fs/2
                cutoffFreq = 1000;
                warndlg('Invalid frequency. Set to default 1000 Hz.');
            end
        end
    end

    function applyFilter(~, ~)
        if isempty(noisy_audio)
            errordlg('Please add noise first.'); return;
        end
        switch filterType
            case 'Low Pass Filter'
                [b, a] = butter(6, cutoffFreq/(Fs/2), 'low');
                filtered_audio = filter(b, a, noisy_audio);
            case 'High Pass Filter'
                [b, a] = butter(6, cutoffFreq/(Fs/2), 'high');
                filtered_audio = filter(b, a, noisy_audio);
            case 'Median Filter'
                filtered_audio = medfilt1(noisy_audio, 3); % Sliding window size of 3
            case 'Adaptive Filter (LMS)'
                mu = 0.01; % Step size, adjustable
                order = 10; % Filter order, adjustable
                filtered_audio = zeros(size(noisy_audio));
                for ch = 1:size(noisy_audio, 2)
                    w = zeros(order, 1); % Initial weights
                    for n = order:length(noisy_audio)
                        x = noisy_audio(max(1, n-order+1):n, ch); % Input vector
                        e = noisy_audio(n, ch) - w' * x; % Error
                        w = w + mu * e * x; % Weight update
                        filtered_audio(n, ch) = e; % Output is error (noise estimate subtracted)
                    end
                end
            case 'Wiener Filter'
                filtered_audio = wiener2(noisy_audio, [3 3]); % 3x3 window for 1D signal approximation
        end
        samples = 0:length(filtered_audio)-1;
        axes(axes3_time);
        plot(samples, filtered_audio, 'g', 'LineWidth', 1.5);
        xlim([0 length(filtered_audio)-1]);
        title(axes3_time, 'Filtered Signal (Amplitude vs Samples)', 'Color', textWhite);
        [s, f, t_spec] = spectrogram(filtered_audio, 256, 200, 256, Fs, 'yaxis');
        axes(axes3_spec);
        imagesc(1:size(s,2), f, 10*log10(abs(s)));
        axis xy;
        colormap jet;
        xlim([1 size(s,2)]);
        colorbar;
        title(axes3_spec, 'Filtered Signal (Frequency vs Samples)', 'Color', textWhite);
    end

    function playNoisy(~, ~)
        if isempty(noisy_audio)
            errordlg('No noisy audio to play.'); return;
        end
        stopPlayback();
        player = audioplayer(noisy_audio, Fs);
        play(player);
    end

    function playFiltered(~, ~)
        if isempty(filtered_audio)
            errordlg('No filtered audio to play.'); return;
        end
        stopPlayback();
        player = audioplayer(filtered_audio, Fs);
        play(player);
    end

    function saveAudio(~, ~)
        if isempty(filtered_audio)
            errordlg('No filtered audio to save.'); return;
        end
        [file, path] = uiputfile('*.wav', 'Save Filtered Audio');
        if file
            audiowrite(fullfile(path, file), filtered_audio, Fs);
            msgbox('Filtered audio saved.');
        end
    end

    function stopPlayback(~, ~)
        if ~isempty(player) && isplaying(player)
            stop(player);
        end
    end

    function reset(~, ~)
        recorded_audio = [];
        noisy_audio = [];
        filtered_audio = [];
        cla(axes1_time); cla(axes1_spec);
        cla(axes2_time); cla(axes2_spec);
        cla(axes3_time); cla(axes3_spec);
        title(axes1_time, 'Original Signal (Amplitude vs Samples)', 'Color', textWhite);
        title(axes1_spec, 'Original Signal (Frequency vs Samples)', 'Color', textWhite);
        title(axes2_time, 'Noisy Signal (Amplitude vs Samples)', 'Color', textWhite);
        title(axes2_spec, 'Noisy Signal (Frequency vs Samples)', 'Color', textWhite);
        title(axes3_time, 'Filtered Signal (Amplitude vs Samples)', 'Color', textWhite);
        title(axes3_spec, 'Filtered Signal (Frequency vs Samples)', 'Color', textWhite);
    end
end
function AudiowatermarkingGUI(parent)

    data = struct();
    data.audio = [];
    data.fs = [];
    data.watermarked = [];
    data.key = [];
    data.msg = '';
    data.fileName = '';
    data.savedPath = '';
bg=[0 0 0];
textColor=[1 1 1];
orange=[1 0.5 0];


 f = uipanel(parent, 'Units', 'normalized', 'Position', [0 0 1 1],'BackgroundColor', [0 0 0]);
  
  % Transmitter section header
   
    uicontrol(f,'Style','text','String','Transmitter Section','FontWeight','bold','FontSize',11, ...
        'Units','normalized','Position',[0.03 0.90 0.3 0.05],'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);

    % Audio file edit + browse
    uicontrol(f,'Style','text','String','Audio File:', 'Units','normalized','Position',[0.03 0.84 0.10 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    audioPathBox = uicontrol(f,'Style','edit','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.14 0.84 0.50 0.04],'String','');
    uicontrol(f,'Style','pushbutton','String','Browse','Units','normalized','Position',[0.66 0.84 0.07 0.04],...
        'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Callback',@browseAudio);

    % Load Audio button 
    uicontrol(f,'Style','pushbutton','String','Load Audio','Units','normalized','Position',[0.74 0.84 0.10 0.04],...
        'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Callback',@loadAudio);

    % Reload button 
    uicontrol(f,'Style','pushbutton','String','Reload','Units','normalized','Position',[0.85 0.84 0.11 0.04],...
        'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Callback',@reloadGUI);

    % Sampling frequency display (transmitter)
    uicontrol(f,'Style','text','String','Sampling Frequency (Hz):','Units','normalized','Position',[0.03 0.79 0.20 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    fsBox = uicontrol(f,'Style','edit','Units','normalized','Position',[0.25 0.79 0.12 0.04],...
        'String','','Enable','inactive','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);

    % Message to embed
    uicontrol(f,'Style','text','String','Message:','Units','normalized','Position',[0.03 0.74 0.08 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    msgBox = uicontrol(f,'Style','edit','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.12 0.74 0.55 0.04],'String','');

    % Secret key (transmitter)
    uicontrol(f,'Style','text','String','Secret Key:','Units','normalized','Position',[0.03 0.69 0.10 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    keyBox = uicontrol(f,'Style','edit','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.14 0.69 0.12 0.04],'String','');

    % Buttons row: Embed, Play Original, Play Watermarked, Save Watermarked
    uicontrol(f,'Style','pushbutton','String','Embed Watermark','Units','normalized','Position',[0.03 0.62 0.22 0.05],...
        'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Callback',@embedWatermark);
    uicontrol(f,'Style','pushbutton','String','Play Original','Units','normalized','Position',[0.27 0.62 0.18 0.05],...
        'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Callback',@playOriginal);
    uicontrol(f,'Style','pushbutton','String','Play Watermarked','Units','normalized','Position',[0.47 0.62 0.18 0.05],...
        'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Callback',@playWatermarked);
    uicontrol(f,'Style','pushbutton','String','Save Watermarked','Units','normalized','Position',[0.67 0.62 0.20 0.05],...
        'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Callback',@saveWatermarked);

    % Saved file path display (transmitter)
    uicontrol(f,'Style','text','String','Saved File Path:','Units','normalized','Position',[0.03 0.57 0.13 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    savedPathBox = uicontrol(f,'Style','edit','Units','normalized','Position',[0.17 0.57 0.77 0.04],...
        'String','','Enable','inactive','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);

    % Receiver section header
    
    uicontrol(f,'Style','text','String','Receiver Section','FontWeight','bold','FontSize',11, ...
        'Units','normalized','Position',[0.03 0.51 0.3 0.05],'ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);

    % Watermarked audio path (receiver) + Browse
    uicontrol(f,'Style','text','String','Watermarked Audio:','Units','normalized','Position',[0.03 0.46 0.18 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    recvAudioBox = uicontrol(f,'Style','edit','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.22 0.46 0.60 0.04],'String','');
    uicontrol(f,'Style','pushbutton','String','Browse','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.84 0.46 0.08 0.04],...
        'Callback',@browseWatermarked);

    % Receiver sampling frequency 
    uicontrol(f,'Style','text','String','Sampling Frequency (Hz):','Units','normalized','Position',[0.03 0.41 0.2 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    recvFsBox = uicontrol(f,'Style','edit','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.24 0.41 0.12 0.04],'String','');

    % Receiver secret key 
    uicontrol(f,'Style','text','String','Secret Key:','Units','normalized','Position',[0.40 0.41 0.10 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    recvSeedBox = uicontrol(f,'Style','edit','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.51 0.41 0.12 0.04],'String','');

    % Extract button
    uicontrol(f,'Style','pushbutton','String','Extract Message','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.03 0.35 0.22 0.05],...
        'Callback',@extractMessage);

    % Extracted message label + box
    uicontrol(f,'Style','text','String','Extracted Message:','Units','normalized','Position',[0.27 0.35 0.18 0.03],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);
    extractedMsgBox = uicontrol(f,'Style','edit','Units','normalized','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0],'Position',[0.46 0.33 0.52 0.07],'Max',2,'String','');

   
    statusBox = uicontrol(f,'Style','text','String','Status: Ready','FontWeight','bold',...
        'Units','normalized','Position',[0.03 0.28 0.94 0.035],...
        'HorizontalAlignment','left','ForegroundColor',[1 1 1],'BackgroundColor',[1 0.5 0]);

 
    axAmp = axes('Parent',f,'Units','normalized','Color', [0 0 0], ...
        'XColor',[1 1 1], 'YColor', [1 1 1],'Position',[0.06 0.045 0.45 0.22], 'NextPlot','add', 'Box','on');
    title(axAmp,'Amplitude (FFT) Spectrum','Color', [1 1 1]);
    xlabel(axAmp,'Frequency (Hz)','Color', textColor);
    ylabel(axAmp,'Magnitude','Color', textColor);



    function browseAudio(~,~)
        [f,p] = uigetfile({'*.wav;*.mp3','Audio Files (*.wav, *.mp3)';'*.*','All Files'}, 'Select audio file');
        if isequal(f,0); return; end
        set(audioPathBox,'String',fullfile(p,f));
    end

    function loadAudio(~,~)
        audioFileStr = get(audioPathBox,'String');
        if isempty(audioFileStr)
            set(statusBox,'String','Status: Please specify an audio filename.'); return;
        end
        if exist(audioFileStr,'file') ~= 2
            set(statusBox,'String','Status: Audio file not found.'); return;
        end
        try
            [y,fs] = audioread(audioFileStr);
        catch
            set(statusBox,'String','Status: Error reading audio file.'); return;
        end
        if size(y,2) > 1
            y = mean(y,2);
        end
        y = y(:);
        y = y(1:min(length(y), fs*10)); 
        if max(abs(y))>0
            y = y / (max(abs(y)) + eps);
        end
        data.audio = y;
        data.fs = fs;
        data.fileName = audioFileStr;
        set(fsBox,'String',num2str(fs));
        set(statusBox,'String',['Status: Loaded "' audioFileStr '" (' num2str(fs) ' Hz)']);
        cla(axAmp); 
        drawnow;
    end

    function embedWatermark(~,~)
        if isempty(data.audio)
            set(statusBox,'String','Status: No audio loaded.'); return;
        end
        data.msg = get(msgBox,'String');
        if isempty(data.msg)
            set(statusBox,'String','Status: Enter a message to embed.'); return;
        end
        keyStr = get(keyBox,'String');
        seed = str2double(keyStr);
        if isnan(seed)
            set(statusBox,'String','Status: Secret key must be numeric.'); return;
        end
        data.key = seed;
        alpha = 0.5;
        rep = 5;

        % bits
        bits_local = [];
        for k = 1:length(data.msg)
            ch = uint8(data.msg(k));
            for b = 8:-1:1
                bits_local = [bits_local, bitget(ch,b)];
            end
        end
        data.bits = bits_local;
        data.nbits = length(bits_local);

        data.halfN = floor(length(data.audio)/2);
        data.total_pairs = data.nbits * rep;

       
        if data.total_pairs > data.halfN
           
            data.total_pairs = data.halfN;
           
        end

        rng(data.key);
        data.bit_positions = randperm(data.halfN, data.total_pairs);

        % FFT & embedding (phase-based)
        X = fft(data.audio);
        magX = abs(X);
        phaseX = angle(X);

        for k = 1:data.total_pairs
            idx = data.bit_positions(k);
            bit_val = data.bits(ceil(k/rep));
            if bit_val == 1
                phaseX(idx) = alpha;
                phaseX(length(X)-idx+2) = -alpha;
            else
                phaseX(idx) = -alpha;
                phaseX(length(X)-idx+2) = alpha;
            end
        end

        X_new = magX .* exp(1i*phaseX);
        water = real(ifft(X_new));
        if max(abs(water))>0
            water = water / (max(abs(water)) + eps);
        end
        data.watermarked = water;

        %  saving filename in current folder
        [pname, fname, ext] = fileparts(data.fileName);
        if isempty(fname), fname = 'audio'; end
        defaultSaved = fullfile(pwd, [fname '_watermarked.wav']);
        try
            audiowrite(defaultSaved, data.watermarked, data.fs);
            data.savedPath = defaultSaved;
            set(savedPathBox,'String',defaultSaved);
            set(statusBox,'String',['Status: Watermark embedded and saved as "' defaultSaved '"']);
        catch
            set(statusBox,'String','Status: Error saving default watermarked file.');
        end

        % plot amplitude  spectra
        Nloc = length(data.audio);
        f_axis = (0:floor(Nloc/2)-1)*(data.fs/Nloc);
        O = fft(data.audio);
        W = fft(data.watermarked);

        % amplitude spectra
        cla(axAmp);
        plot(axAmp, f_axis, abs(O(1:floor(Nloc/2)))); hold(axAmp,'on');
        plot(axAmp, f_axis, abs(W(1:floor(Nloc/2))), 'r'); hold(axAmp,'off');
        legend(axAmp,'Original','Watermarked');
        title(axAmp,'Amplitude (FFT) Spectrum');
        xlabel(axAmp,'Frequency (Hz)'); ylabel(axAmp,'Magnitude');

        drawnow;
    end

    function saveWatermarked(~,~)
        if isempty(data.watermarked)
            set(statusBox,'String','Status: No watermarked audio to save. Embed first.'); return;
        end
        [f,p] = uiputfile('watermarked_audio.wav','Save watermarked as');
        if isequal(f,0); return; end
        try
            audiowrite(fullfile(p,f), data.watermarked, data.fs);
            data.savedPath = fullfile(p,f);
            set(savedPathBox,'String',data.savedPath);
            set(statusBox,'String',['Status: Watermarked saved as "' fullfile(p,f) '"']);
        catch
            set(statusBox,'String','Status: Error saving watermarked audio.');
        end
    end

    function playOriginal(~,~)
        if isempty(data.audio)
            set(statusBox,'String','Status: No audio loaded.'); return;
        end
        set(statusBox,'String','Status: Playing original audio...');
        try
            sound(data.audio, data.fs);
        catch
            set(statusBox,'String','Status: Playback error (original).');
        end
    end

    function playWatermarked(~,~)
        if isempty(data.watermarked)
            set(statusBox,'String','Status: No watermarked audio. Embed first.'); return;
        end
        set(statusBox,'String','Status: Playing watermarked audio...');
        try
            sound(data.watermarked, data.fs);
        catch
            set(statusBox,'String','Status: Playback error (watermarked).');
        end
    end

    function browseWatermarked(~,~)
        [f,p] = uigetfile({'*.wav;*.mp3','Audio Files (*.wav,*.mp3)';'*.*','All Files'}, 'Select watermarked file');
        if isequal(f,0); return; end
        full = fullfile(p,f);
        set(recvAudioBox,'String',full);
        try
            [y, fs] = audioread(full);
            if size(y,2)>1, y = mean(y,2); end
            data.watermarked = y(:);
            data.fs = fs;
            set(statusBox,'String',['Status: Loaded "' full '" for extraction']);
            if ~isempty(data.audio)
                Nloc = length(data.audio);
                f_axis = (0:floor(Nloc/2)-1)*(data.fs/Nloc);
                O = fft(data.audio);
                W = fft(data.watermarked);
                cla(axAmp);
                plot(axAmp, f_axis, abs(O(1:floor(Nloc/2)))); hold(axAmp,'on');
                plot(axAmp, f_axis, abs(W(1:floor(Nloc/2))), 'r'); hold(axAmp,'off');
                legend(axAmp,'Original','Watermarked');
                
                 drawnow;
            end
        catch
            set(statusBox,'String','Status: Error reading watermarked file.');
        end
    end

    function extractMessage(~,~)
        if isempty(data.watermarked)
            set(statusBox,'String','Status: No watermarked audio present. Embed or load first.'); return;
        end
        recv_fs = str2double(get(recvFsBox,'String'));
        recv_seed = str2double(get(recvSeedBox,'String'));
        if isnan(recv_fs) || isnan(recv_seed)
            set(statusBox,'String','Status: Receiver sampling freq and seed must be numeric.'); return;
        end

       
        if isempty(data.fs)
            set(statusBox,'String','Status: Original sampling frequency unknown.'); return;
        end
        fsFactor = abs(recv_fs - data.fs)/data.fs;

       
        if recv_seed ~= data.key
            set(statusBox,'String','Status: Wrong Key ? Extraction failed.'); 
     
        end
        if fsFactor > 0
            set(statusBox,'String','Status: Wrong Fs ? Extracted message likely distorted.');
        else
            set(statusBox,'String','Status: Access Granted. Extracting message...');
        end

        rng_seed_for_positions = recv_seed + round(fsFactor * 1e6);
        rng(rng_seed_for_positions);

        recv = data.watermarked(:); 
        R = fft(recv);
        phaseR = angle(R);

        halfN_local = floor(length(recv)/2);
        total_pairs_rx = data.nbits * 5; % rep = 5
        if total_pairs_rx > halfN_local
     
            total_pairs_rx = min(total_pairs_rx, halfN_local);
        end
        bit_positions_rx = randperm(halfN_local, total_pairs_rx);

        pair_bits = zeros(1,total_pairs_rx);
        for k = 1:total_pairs_rx
            idx = bit_positions_rx(k);
            if idx <= length(phaseR)
                pair_bits(k) = phaseR(idx) > 0;
            else
                pair_bits(k) = 0;
            end
        end

        % Majority voting
        nbits_local = data.nbits;
        recovered_bits = zeros(1, nbits_local);
        for i = 1:nbits_local
            idxs = (i-1)*5 + (1:5);
            if max(idxs) <= length(pair_bits)
                recovered_bits(i) = round(mean(pair_bits(idxs)));
            else
             
                recovered_bits(i) = 0;
            end
        end

        % Bits to message
        nchars = floor(nbits_local/8);
        recovered_msg = '';
        for c = 1:nchars
            val = 0;
            for b = 1:8
                val = val + recovered_bits((c-1)*8 + b) * 2^(8-b);
            end
            recovered_msg = [recovered_msg char(val)];
        end

        set(extractedMsgBox,'String',recovered_msg);

        
        if recv_seed ~= data.key
            set(statusBox,'String','Status: Wrong Key ? Extraction failed.');
        elseif fsFactor > 0
            set(statusBox,'String','Status: Wrong Fs ? Extracted message likely distorted.');
        else
            set(statusBox,'String',['Status: Extracted Message: "' recovered_msg '"']);
        end
    end

    function reloadGUI(~,~)
        % Reset GUI 
        set(audioPathBox,'String','');
        set(fsBox,'String','');
        set(msgBox,'String','');
        set(keyBox,'String','');
        set(savedPathBox,'String','');
        set(recvAudioBox,'String','');
        set(recvFsBox,'String','');
        set(recvSeedBox,'String','');
        set(extractedMsgBox,'String','');
        set(statusBox,'String','Status: Ready');
        data = struct('audio',[],'fs',[],'watermarked',[],'key',[],'msg','','fileName','','savedPath','');
        cla(axAmp);  
        drawnow;
    end

end
