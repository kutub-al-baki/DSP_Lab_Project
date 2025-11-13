function CombinedAudioGUIfinal
    fig = figure('Name', 'Combined Audio Processing GUI', ...
                 'NumberTitle', 'on', ...
                 'Position', [50, 50, 1400, 750], ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'none', ...
                 'Resize', 'on');

    tabgroup = uitabgroup(fig);
    tab1 = uitab(tabgroup, 'Title', 'Audio File Mode');
    tab2 = uitab(tabgroup, 'Title', 'Voice Recorder Mode');
   
    
    panel1 = uipanel('Parent', tab1, ...
                     'Units', 'normalized', ...
                     'Position', [0 0 1 1], ...
                     'BackgroundColor', [1 0.5 0]);

    panel2 = uipanel('Parent', tab2, ...
                     'Units', 'normalized', ...
                     'Position', [0 0 1 1], ...
                     'BackgroundColor', [0 0.447 0.741]);


    AudioEffectAndNoiseRemovalGUI(tab1);
    VoiceRecorderGUI(tab2);
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
                    delay = 0.3;
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