% QE_EOG()
% This function computes onset and offset of the QE period, given EOG data
%
% INPUT: 
%
% timeSeries: a 1D vector of EOG values across time, or an array where the
% first dimension is time and remaining dimensions are trials/repetitions.
% The time series can represent any EOG channel (e.g., horizontal,
% vertical).
%
% timeVector: a 1D vector of time points in seconds. This vector must be of
% the same length of 'timeSeries'. Values must range from negative to
% positive. Time = 0 seconds is the reference point (e.g., movement
% initiation)
%
% algorithmChoice: structure with the following fields
%  - name: either 'dispersion' or 'velocity'
%  - winlen: filter length in samples (points). To express seconds, divide
%  by the sampling rate. This value is used for a median filter if using
%  'dispersion' or for a SG 1st degree differentiation filter if using
%  'velocity'
%  - polynDeg: degree of the polynomial [only used with the 'velocity'
%  algorithm, irrelevant when using the 'dispersion' algorithm]
%  - threshold: choose one numeric value in the same units as the
%  'timeSeries'. For example, a value of 3 would be a good choice if the
%  timeSeries is in degrees of visual angle. Alternatively, use 'auto' or
%  leave empty to automatically estimate threshold based on the most stable
%  window in the signal (3× standard deviation of minimum variance window).
%
% showPlot: logical true or false / A potentially useful figure
%
%
%
% OUTPUT:
%
% [QEonset, QEoffset]: vectors (1 x nTrials) of onset and offset times in seconds
%
% code written by Germano Gallicchio (Bangor University) in May 2022
% germano.gallicchio@gmail.com

function [QEonset, QEoffset] = QE_EOG(timeSeries, timeVector, algorithmChoice, showPlot)

                            
% sanity checks on the timeVector
if sum(timeVector==0)~=1;  error('The time vector must have a 0'); end
if timeVector(1)>=0; error('The time vector must start with a negative value'); end
if timeVector(end)<=0; error('The time vector must end with a positive value'); end
if mod(numel(timeSeries), length(timeVector))~=0; error('timeSeries length must be a multiple of timeVector length'); end

% extract the sampling rate from timeVector
dt = mean(diff(timeVector));    % time between consecutive samples
EOG_srate = 1/dt;               % sampling rate

% identify movement initiation
start_pnt = find(timeVector==0);

% reshape input to support multiple trials/repetitions in higher dimensions
nSamples = length(timeVector);
if mod(numel(timeSeries), nSamples)~=0
    error('timeSeries length must be a multiple of timeVector length');
end
timeSeries_mat = reshape(timeSeries, nSamples, []); % columns = trials
nTrials = size(timeSeries_mat,2);

QEonset  = nan(1, nTrials);
QEoffset = nan(1, nTrials);

if showPlot && nTrials>1
    warning('showPlot:multiTrial','Multiple trials provided; plotting only the first trial.');
end

% base threshold (may be empty or scalar or vector)
if isfield(algorithmChoice,'threshold')
    base_threshold = algorithmChoice.threshold;
else
    base_threshold = [];
end

for trialIdx = 1:nTrials
    ts = timeSeries_mat(:,trialIdx);

    % determine threshold for this trial (scalar or per-trial vector or auto)
    thr_val = base_threshold;
    if isnumeric(thr_val) && numel(thr_val)>1
        if numel(thr_val)~=nTrials
            error('If threshold is a numeric vector, it must have one value per trial.');
        end
        thr_val = thr_val(trialIdx);
    end

    % Automatic threshold estimation if not provided
    if ~isfield(algorithmChoice, 'threshold') || isempty(thr_val) || (ischar(thr_val) && strcmp(thr_val, 'auto'))
        % Find the most stable window in the signal (minimum robust variance)
        win_duration = 0.200; % 200ms window for stability assessment
        win_samples = round(win_duration * EOG_srate);
        
        if strcmp(algorithmChoice.name, 'velocity')
            % For velocity: compute velocity first, then find stable window
            vel_temp = [0; diff(ts)]/dt;
            
            % Sliding window robustness metric: squared scaled MAD (more robust than variance)
            num_windows = length(vel_temp) - win_samples + 1;
            window_mad2 = zeros(num_windows, 1);
            for i = 1:num_windows
                sigma_hat = mad(vel_temp(i:i+win_samples-1), 1); % scaled MAD ~ std
                window_mad2(i) = sigma_hat.^2;                   % variance proxy
            end
            
            % Find most stable window (minimum robust variance proxy)
            [~, min_idx] = min(window_mad2);
            stable_window = vel_temp(min_idx:min_idx+win_samples-1);
            
            % Threshold = 3× std of most stable window
            thr_val = 3 * std(stable_window);
        else
            % For dispersion: find most stable window in position using robust metric
            num_windows = length(ts) - win_samples + 1;
            window_mad2 = zeros(num_windows, 1);
            for i = 1:num_windows
                sigma_hat = mad(ts(i:i+win_samples-1), 1); % scaled MAD ~ std
                window_mad2(i) = sigma_hat.^2;             % variance proxy
            end
            
            % Find most stable window (minimum robust variance proxy)
            [~, min_idx] = min(window_mad2);
            stable_window = ts(min_idx:min_idx+win_samples-1);
            
            % Threshold = 3× std of most stable window
            thr_val = 3 * std(stable_window);
        end
        
        if showPlot && trialIdx==1
            fprintf('Auto-estimated threshold (trial %d): %.2f\n', trialIdx, thr_val);
        end
    end

    switch algorithmChoice.name
        
        case 'dispersion'
            
            % apply the median filter
            ts_filt = movmedian(ts,algorithmChoice.winlen);
            
            % center the time series to the point with time = 0 s
            ts_centered = ts_filt - ts_filt(start_pnt);
            
            % boundaries
            timeSeries_upper = ts_centered + thr_val;
            timeSeries_lower = ts_centered - thr_val;
            
            % 100 ms rule
            % In lay terms: make it harder for the upper series to cross y=0 by
            % using a moving max window; if it still crosses, we know there is a
            % full window under y=0. Same for the lower series with a moving min.
            % The window is [0 k] before time=0 and [k 0] after time=0, k = 100 ms.
            timeWin_sec = 0.100;                    % in s
            timeWin_pnt = timeWin_sec * EOG_srate;  % in pnt
            timeSeries_upper_pre  = movmax(timeSeries_upper, [0 timeWin_pnt], 'Endpoints','shrink');
            timeSeries_upper_post = movmax(timeSeries_upper, [timeWin_pnt 0], 'Endpoints','shrink');
            timeSeries_upper = [timeSeries_upper_pre(timeVector<=0); timeSeries_upper_post(timeVector>0)];
            timeSeries_lower_pre  = movmin(timeSeries_lower,  [0 timeWin_pnt], 'Endpoints','shrink');
            timeSeries_lower_post = movmin(timeSeries_lower,  [timeWin_pnt 0], 'Endpoints','shrink');
            timeSeries_lower = [timeSeries_lower_pre(timeVector<=0); timeSeries_lower_post(timeVector>0)];
            
            % QEonset
            qe_onset_pnt = find([1:length(timeVector)]'<=start_pnt & (timeSeries_lower>0 | timeSeries_upper<0),1,'last');
            if qe_onset_pnt~=start_pnt
                qe_onset_pnt = qe_onset_pnt + 1;
            end
            if isempty(qe_onset_pnt)
                qe_onset_pnt = 1;
            end
            qe_onset_sec = timeVector(qe_onset_pnt);
            
            % QEoffset
            qe_offset_pnt = find([1:length(timeVector)]'>=start_pnt & (timeSeries_lower>0 | timeSeries_upper<0),1,'first');
            if qe_offset_pnt~=start_pnt
                qe_offset_pnt = qe_offset_pnt - 1;
            end
            if isempty(qe_offset_pnt)
                qe_offset_pnt = length(timeVector);
            end
            qe_offset_sec = timeVector(qe_offset_pnt);
            
            % sanity check: figure (first trial only if multiple)
            plotFlag = showPlot && trialIdx==1;
            if plotFlag
                figure; clf;
                tl = tiledlayout(2,1,'TileSpacing','compact');
                
                % Top tile: dispersion profiles
                nexttile;
                y_all = [ts_centered(:); timeSeries_upper(:); timeSeries_lower(:)];
                y_pad = 0.05 * range(y_all);
                if ~isfinite(y_pad); y_pad = 0; end
                y_lim_disp = [min(y_all)-y_pad, max(y_all)+y_pad];
                
                h_signal = plot(timeVector,ts_centered, 'Color', [0 0 0],'LineWidth',1.5); hold on
                h_upper  = plot(timeVector,timeSeries_upper, 'Color', [1 0 0],'LineWidth',1);
                h_lower  = plot(timeVector,timeSeries_lower, 'Color', [0 0 1],'LineWidth',1);
                line([0 0], y_lim_disp,'Color', [0 0 0 0.2],'LineStyle','--')
                line(timeVector([1 end]), [0 0],'Color', [0 0 0 0.2],'LineStyle','--')
                onset_line  = line(repmat(timeVector(qe_onset_pnt),1,2),y_lim_disp, 'Color',[1 1 0],'LineWidth',2,'LineStyle','-');
                offset_line = line(repmat(timeVector(qe_offset_pnt),1,2),y_lim_disp,'Color',[0 1 1],'LineWidth',2,'LineStyle','-');
                
                ylabel('Dispersion');
                set(gca,'XLim',[timeVector(1) timeVector(end)], 'YLim', y_lim_disp);
                grid on; grid minor;
                legend([h_signal h_upper h_lower onset_line offset_line],["Filtered position" "Upper bound" "Lower bound" "QE onset" "QE offset"],'Location','best')
                
                % Bottom tile: positional EOG
                nexttile;
                y_pos_norm = normalize(ts,'range');
                y_pos_lim = [min(y_pos_norm)-0.05*range(y_pos_norm), max(y_pos_norm)+0.05*range(y_pos_norm)];
                
                h_pos = plot(timeVector,y_pos_norm, 'Color', [0 0 0],'LineWidth',1.5); hold on
                line([0 0], y_pos_lim,'Color', [0 0 0 0.2],'LineStyle','--')
                line(timeVector([1 end]), [0.5 0.5],'Color', [0 0 0 0.2],'LineStyle','--')
                onset_line2  = line(repmat(timeVector(qe_onset_pnt),1,2),y_pos_lim, 'Color',[1 1 0],'LineWidth',2,'LineStyle','-');
                offset_line2 = line(repmat(timeVector(qe_offset_pnt),1,2),y_pos_lim,'Color',[0 1 1],'LineWidth',2,'LineStyle','-');
                
                ylabel('Position (normalized)');
                xlabel('Time (s)');
                set(gca,'XLim',[timeVector(1) timeVector(end)], 'YLim', y_pos_lim);
                grid on; grid minor;
                legend([h_pos onset_line2 offset_line2],["EOG position" "QE onset" "QE offset"],'Location','best')
                
                title(tl, sprintf('Dispersion QE Analysis (trial %d)', trialIdx));
            end
            
        
        case 'velocity'

            % apply Savitzky-Golay differentiation filter
            deriv_order = 1;
            sg_framelen = algorithmChoice.winlen;
            sg_order    = algorithmChoice.polynDeg;
            [~, sg_coeffs] = sgolay(sg_order, sg_framelen);
            timeSeries_D = conv(ts, factorial(deriv_order)/(-dt)^deriv_order * sg_coeffs(:,deriv_order+1), 'same');

            % center the time series to the point with time = 0s
            timeSeries_D = timeSeries_D - timeSeries_D(start_pnt);
            
            % boundaries
            timeSeries_D_upper = timeSeries_D + thr_val;
            timeSeries_D_lower = timeSeries_D - thr_val;
            

            % QEonset
            qe_onset_pnt = find([1:length(timeVector)]'<=start_pnt & (timeSeries_D_lower>0 | timeSeries_D_upper<0),1,'last');
            if qe_onset_pnt~=start_pnt
                qe_onset_pnt = qe_onset_pnt + 1;
            end
            if isempty(qe_onset_pnt)
                qe_onset_pnt = 1;
            end
            qe_onset_sec = timeVector(qe_onset_pnt);
            
            % QEoffset
            qe_offset_pnt = find([1:length(timeVector)]'>=start_pnt & (timeSeries_D_lower>0 | timeSeries_D_upper<0),1,'first');
            if qe_offset_pnt~=start_pnt
                qe_offset_pnt = qe_offset_pnt - 1;
            end
            if isempty(qe_offset_pnt)
                qe_offset_pnt = length(timeVector);
            end
            qe_offset_sec = timeVector(qe_offset_pnt);
            
            
            % sanity check: figure (first trial only if multiple)
            plotFlag = showPlot && trialIdx==1;
            if plotFlag
                figure; clf;
                tl = tiledlayout(2,1,'TileSpacing','compact');
                
                % Top tile: velocity profiles
                nexttile;
                y_all = [timeSeries_D(:); timeSeries_D_upper(:); timeSeries_D_lower(:)];
                y_pad = 0.05 * range(y_all);
                if ~isfinite(y_pad); y_pad = 0; end
                y_lim_vel = [min(y_all)-y_pad, max(y_all)+y_pad];
                
                h_vel       = plot(timeVector,timeSeries_D, 'Color', [0 1 0],'LineWidth',1.5); hold on
                h_vel_upper = plot(timeVector,timeSeries_D_upper, 'Color', [1 0 0],'LineWidth',1);
                h_vel_lower = plot(timeVector,timeSeries_D_lower, 'Color', [0 0 1],'LineWidth',1);
                line([0 0], y_lim_vel,'Color', [0 0 0 0.2],'LineStyle','--')
                line(timeVector([1 end]), [0 0],'Color', [0 0 0 0.2],'LineStyle','--')
                onset_line  = line(repmat(timeVector(qe_onset_pnt),1,2),y_lim_vel, 'Color',[1 1 0],'LineWidth',2,'LineStyle','-');
                offset_line = line(repmat(timeVector(qe_offset_pnt),1,2),y_lim_vel,'Color',[0 1 1],'LineWidth',2,'LineStyle','-');
                
                ylabel('Velocity');
                set(gca,'XLim',[timeVector(1) timeVector(end)], 'YLim', y_lim_vel);
                grid on; grid minor;
                legend([h_vel h_vel_upper h_vel_lower onset_line offset_line],["Velocity" "Upper bound" "Lower bound" "QE onset" "QE offset"],'Location','best')
                
                % Bottom tile: positional EOG
                nexttile;
                y_pos_norm = normalize(ts,'range');
                y_pos_lim = [min(y_pos_norm)-0.05*range(y_pos_norm), max(y_pos_norm)+0.05*range(y_pos_norm)];
                
                h_pos = plot(timeVector,y_pos_norm, 'Color', [0 0 0],'LineWidth',1.5); hold on
                line([0 0], y_pos_lim,'Color', [0 0 0 0.2],'LineStyle','--')
                line(timeVector([1 end]), [0.5 0.5],'Color', [0 0 0 0.2],'LineStyle','--')
                onset_line2  = line(repmat(timeVector(qe_onset_pnt),1,2),y_pos_lim, 'Color',[1 1 0],'LineWidth',2,'LineStyle','-');
                offset_line2 = line(repmat(timeVector(qe_offset_pnt),1,2),y_pos_lim,'Color',[0 1 1],'LineWidth',2,'LineStyle','-');
                
                ylabel('Position (normalized)');
                xlabel('Time (s)');
                set(gca,'XLim',[timeVector(1) timeVector(end)], 'YLim', y_pos_lim);
                grid on; grid minor;
                legend([h_pos onset_line2 offset_line2],["EOG position" "QE onset" "QE offset"],'Location','best')
                
                title(tl, sprintf('Velocity QE Analysis (trial %d)', trialIdx));
            end
            
                    
        otherwise
            error('select a valid algorithm')
    end

    % store outputs
    QEonset(trialIdx)  = qe_onset_sec;
    QEoffset(trialIdx) = qe_offset_sec;

    % sanity checks (for debugging)
    if QEonset(trialIdx)>0; error('QEonset cannot be positive (i.e., after time=0)'); end
    if QEoffset(trialIdx)<0; error('QEoffset cannot be negative (i.e., before time=0)'); end
end







