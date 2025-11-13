function gui_phase_watermarking_final_resizable_4
% PHASE-BASED AUDIO WATERMARKING GUI (Resizable + Phase plots + Reload)
% This file preserves your working logic and restores missing callbacks,
% while adding the reload button and phase-spectrum plotting.
% Compatible with MATLAB ~2015a (uses uicontrol & axes).

    % -------------------------
    % Initialize data
    % -------------------------
    data = struct();
    data.audio = [];
    data.fs = [];
    data.watermarked = [];
    data.key = [];
    data.msg = '';
    data.fileName = '';
    data.savedPath = '';

    % -------------------------
    % Create main figure
    % -------------------------
    f = figure('Name','Phase-Based Audio Watermarking',...
               'Units','normalized',...
               'Position',[0.15 0.10 0.70 0.78],...
               'NumberTitle','off','MenuBar','none','Resize','on',...
               'Color',[0.95 0.95 0.95]);

    % -------------------------
    % Transmitter section header
    % -------------------------
    uicontrol('Style','text','String','Transmitter Section','FontWeight','bold','FontSize',11, ...
        'Units','normalized','Position',[0.03 0.88 0.3 0.05],'BackgroundColor',[0.95 0.95 0.95]);

    % Audio file edit + browse
    uicontrol('Style','text','String','Audio File:', 'Units','normalized','Position',[0.03 0.84 0.10 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    audioPathBox = uicontrol('Style','edit','Units','normalized','Position',[0.13 0.84 0.50 0.04],'String','');
    uicontrol('Style','pushbutton','String','Browse','Units','normalized','Position',[0.64 0.84 0.07 0.04],...
        'Callback',@browseAudio);

    % Load Audio button (keeps original semantics)
    uicontrol('Style','pushbutton','String','Load Audio','Units','normalized','Position',[0.73 0.84 0.10 0.04],...
        'Callback',@loadAudio);

    % Reload button - top-right
    uicontrol('Style','pushbutton','String','Reload','Units','normalized','Position',[0.85 0.84 0.11 0.04],...
        'BackgroundColor',[0.90 0.93 0.98],'Callback',@reloadGUI);

    % Sampling frequency display (transmitter)
    uicontrol('Style','text','String','Sampling Frequency (Hz):','Units','normalized','Position',[0.03 0.79 0.20 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    fsBox = uicontrol('Style','edit','Units','normalized','Position',[0.25 0.79 0.12 0.04],...
        'String','','Enable','inactive','BackgroundColor',[0.92 0.92 0.92]);

    % Message to embed
    uicontrol('Style','text','String','Message:','Units','normalized','Position',[0.03 0.74 0.08 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    msgBox = uicontrol('Style','edit','Units','normalized','Position',[0.12 0.74 0.55 0.04],'String','');

    % Secret key (transmitter) - initially blank
    uicontrol('Style','text','String','Secret Key:','Units','normalized','Position',[0.03 0.69 0.10 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    keyBox = uicontrol('Style','edit','Units','normalized','Position',[0.13 0.69 0.12 0.04],'String','');

    % Buttons row: Embed, Play Original, Play Watermarked, Save Watermarked
    uicontrol('Style','pushbutton','String','Embed Watermark','Units','normalized','Position',[0.03 0.62 0.22 0.05],...
        'Callback',@embedWatermark);
    uicontrol('Style','pushbutton','String','Play Original','Units','normalized','Position',[0.27 0.62 0.18 0.05],...
        'Callback',@playOriginal);
    uicontrol('Style','pushbutton','String','Play Watermarked','Units','normalized','Position',[0.47 0.62 0.18 0.05],...
        'Callback',@playWatermarked);
    uicontrol('Style','pushbutton','String','Save Watermarked','Units','normalized','Position',[0.67 0.62 0.20 0.05],...
        'Callback',@saveWatermarked);

    % Saved file path display (transmitter)
    uicontrol('Style','text','String','Saved File Path:','Units','normalized','Position',[0.03 0.57 0.13 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    savedPathBox = uicontrol('Style','edit','Units','normalized','Position',[0.17 0.57 0.77 0.04],...
        'String','','Enable','inactive','BackgroundColor',[0.92 0.92 0.92]);

    % -------------------------
    % Receiver section header
    % -------------------------
    uicontrol('Style','text','String','Receiver Section','FontWeight','bold','FontSize',11, ...
        'Units','normalized','Position',[0.03 0.50 0.3 0.05],'BackgroundColor',[0.95 0.95 0.95]);

    % Watermarked audio path (receiver) + Browse
    uicontrol('Style','text','String','Watermarked Audio:','Units','normalized','Position',[0.03 0.46 0.18 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    recvAudioBox = uicontrol('Style','edit','Units','normalized','Position',[0.20 0.46 0.60 0.04],'String','');
    uicontrol('Style','pushbutton','String','Browse','Units','normalized','Position',[0.82 0.46 0.08 0.04],...
        'Callback',@browseWatermarked);

    % Receiver sampling frequency box - initially blank
    uicontrol('Style','text','String','Sampling Frequency (Hz):','Units','normalized','Position',[0.03 0.41 0.2 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    recvFsBox = uicontrol('Style','edit','Units','normalized','Position',[0.24 0.41 0.12 0.04],'String','');

    % Receiver secret key - initially blank
    uicontrol('Style','text','String','Secret Key:','Units','normalized','Position',[0.40 0.41 0.10 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    recvSeedBox = uicontrol('Style','edit','Units','normalized','Position',[0.50 0.41 0.12 0.04],'String','');

    % Extract button
    uicontrol('Style','pushbutton','String','Extract Message','Units','normalized','Position',[0.03 0.35 0.22 0.05],...
        'Callback',@extractMessage);

    % Extracted message label + box
    uicontrol('Style','text','String','Extracted Message:','Units','normalized','Position',[0.27 0.35 0.18 0.03],...
        'HorizontalAlignment','left','BackgroundColor',[0.95 0.95 0.95]);
    extractedMsgBox = uicontrol('Style','edit','Units','normalized','Position',[0.45 0.33 0.52 0.07],'Max',2,'String','');

    % Status box (single-line)
    statusBox = uicontrol('Style','text','String','Status: Ready','FontWeight','bold',...
        'Units','normalized','Position',[0.03 0.28 0.94 0.035],...
        'HorizontalAlignment','left','BackgroundColor',[0.88 0.95 0.88]);

    % -------------------------
    % Spectra axes: left amplitude, right phase
    % -------------------------
    axAmp = axes('Units','normalized','Position',[0.07 0.03 0.45 0.22]);
    title(axAmp,'Amplitude (FFT) Spectrum');
    xlabel(axAmp,'Frequency (Hz)');
    ylabel(axAmp,'Magnitude');

    axPhase = axes('Units','normalized','Position',[0.55 0.03 0.40 0.22]);
    title(axPhase,'Phase Spectrum');
    xlabel(axPhase,'Frequency (Hz)');
    ylabel(axPhase,'Phase (radians)');

    % -------------------------------------------------------------------------
    % CALLBACKS (preserve original logic; restored missing functions and added
    % spectrum plotting and reload functionality)
    % -------------------------------------------------------------------------

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
        y = y(1:min(length(y), fs*10)); % 10 sec limit
        if max(abs(y))>0
            y = y / max(abs(y));
        end
        data.audio = y;
        data.fs = fs;
        data.fileName = audioFileStr;
        set(fsBox,'String',num2str(fs));
        set(statusBox,'String',['Status: Loaded "' audioFileStr '" (' num2str(fs) ' Hz)']);
        cla(axAmp); cla(axPhase);
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
            water = water / max(abs(water));
        end
        data.watermarked = water;

        % default saved filename in current folder
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

        % plot amplitude & phase spectra
        Nloc = length(data.audio);
        f_axis = (0:Nloc/2-1)*(data.fs/Nloc);
        O = fft(data.audio);
        W = fft(data.watermarked);

        % amplitude (magnitude) on left
        axes(axAmp);
        cla(axAmp);
        plot(f_axis, abs(O(1:Nloc/2))); hold(axAmp,'on');
        plot(f_axis, abs(W(1:Nloc/2)), 'r'); hold(axAmp,'off');
        legend(axAmp,'Original','Watermarked');
        title(axAmp,'Amplitude (FFT) Spectrum');
        xlabel(axAmp,'Frequency (Hz)'); ylabel(axAmp,'Magnitude');

        % phase on right
        axes(axPhase);
        cla(axPhase);
        plot(f_axis, unwrap(angle(O(1:Nloc/2)))); hold(axPhase,'on');
        plot(f_axis, unwrap(angle(W(1:Nloc/2))), 'r'); hold(axPhase,'off');
        legend(axPhase,'Original','Watermarked');
        title(axPhase,'Phase Spectrum');
        xlabel(axPhase,'Frequency (Hz)'); ylabel(axPhase,'Phase (radians)');
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
            % Update amplitude & phase plots for the loaded file (if original present)
            if ~isempty(data.audio)
                Nloc = length(data.audio);
                f_axis = (0:Nloc/2-1)*(data.fs/Nloc);
                O = fft(data.audio);
                W = fft(data.watermarked);
                axes(axAmp); cla(axAmp);
                plot(f_axis, abs(O(1:Nloc/2))); hold(axAmp,'on');
                plot(f_axis, abs(W(1:Nloc/2)), 'r'); hold(axAmp,'off');
                legend(axAmp,'Original','Watermarked');
                axes(axPhase); cla(axPhase);
                plot(f_axis, unwrap(angle(O(1:Nloc/2)))); hold(axPhase,'on');
                plot(f_axis, unwrap(angle(W(1:Nloc/2))), 'r'); hold(axPhase,'off');
                legend(axPhase,'Original','Watermarked');
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

        % compute fs mismatch factor
        if isempty(data.fs)
            set(statusBox,'String','Status: Original sampling frequency unknown.'); return;
        end
        fsFactor = abs(recv_fs - data.fs)/data.fs;

        % indicate mismatch if any
        if recv_seed ~= data.key
            set(statusBox,'String','Status: Wrong Key ? Extraction failed.'); 
            % still attempt extraction (will be wrong)
        end
        if fsFactor > 0
            % notify user — extraction will likely be distorted
            set(statusBox,'String','Status: Wrong Fs ? Extracted message likely distorted.');
        else
            set(statusBox,'String','Status: Access Granted. Extracting message...');
        end

        % Use receiver seed but perturb with fsFactor to desynchronize when Fs is wrong
        rng_seed_for_positions = recv_seed + round(fsFactor * 1e6);
        rng(rng_seed_for_positions);

        recv = data.watermarked(:); % extraction from watermarked signal (no noise)
        R = fft(recv);
        phaseR = angle(R);

        halfN_local = floor(length(recv)/2);
        total_pairs_rx = data.nbits * 5; % rep = 5
        % Guard: if randperm size requested exceeds halfN_local, fallback
        if total_pairs_rx > halfN_local
            % fallback to using available indices (this case shouldn't normally happen)
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
                % if out of bounds, set to zero
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

        % Final status message (already set above for mismatches, refine)
        if recv_seed ~= data.key
            set(statusBox,'String','Status: Wrong Key ? Extraction failed.');
        elseif fsFactor > 0
            set(statusBox,'String','Status: Wrong Fs ? Extracted message likely distorted.');
        else
            set(statusBox,'String',['Status: Extracted Message: "' recovered_msg '"']);
        end
    end

    function reloadGUI(~,~)
        % Reset GUI fields, stored data, and plots
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
        cla(axAmp); cla(axPhase);
    end

end
