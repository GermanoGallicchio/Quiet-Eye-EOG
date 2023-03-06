% QE_EOG()
% This function computes onset and offset of the QE period, given EOG data
%
% INPUT: 
%
% timeSeries: a 1D vector of EOG values across time. The time series
% can represent any EOG channel (e.g., horizontal, vertical).
%
% timeVector: a 1D vector of time points in seconds. This vector must be of
% the same length of 'timeSeries'. Values must range from negative to
% positive. Time = 0 seconds is the reference point (e.g., movement
% initiation)
%
% algorithmChoice: structure with the following fields
%  - name: either 'dispersion' or 'velocity'
%  - winlen: filter order in time points (in seconds = time points divided
%  by sampling rate). This value is used for a median filter if using
%  'dispersion' or for a SG 1st degree differentiation filter if using
%  'velocity'
%  - polynDeg: degree of the polynomial [only used with the 'velocity'
%  algorithm, irrelevant when using the 'dispersion' algorithm]
%  - threshold: choose one numeric value in the same units as the
%  'timeSeries'. For example, a value of 3 would be a good choice if the
%  timeSeries is in degrees of visual angle.
%
% doYouWantTheFigure: choose "yes" or "no" / A potentially useful figure
%
%
%
% OUTPUT:
%
% [QEonset, QEoffset]: a 1D vector of two time points in seconds
%
% code written by Germano Gallicchio (Bangor University) in May 2022
% germano.gallicchio@gmail.com

function [QEonset, QEoffset] = QE_EOG(timeSeries, timeVector, algorithmChoice, doYouWantTheFigure)

                            
% sanity checks on the timeVector
if sum(timeVector==0)~=1;  error('The time vector must have a 0'); end
if timeVector(1)>=0; error('The time vector must start with a negative value'); end
if timeVector(end)<=0; error('The time vector must end with a positive value'); end
if length(timeVector)~=length(timeSeries); error('The time vector must be of the same length as the time series'); end

% extract the sampling rate from timeVector
dt = mean(diff(timeVector));    % time between consecutive samples
EOG_srate = 1/dt;               % sampling rate

% identify movement initiation
start_pnt = find(timeVector==0);

switch algorithmChoice.name
    
    case 'dispersion'
        
        % apply the median filter
        timeSeries = movmedian(timeSeries,algorithmChoice.winlen);
        
        % center the time series to the point with time = 0 s
        timeSeries = timeSeries - timeSeries(start_pnt);
        
        % finding the timeseries exceeding a certain threshold is equivalent to
        % computing two boundaries time series (upper and lower, by
        % subracting/adding the threshold to the dataseries) and finding
        % when one of the two boundaries crosses the y=0 line
        threshold = algorithmChoice.threshold;
        timeSeries_upper = timeSeries + threshold;
        timeSeries_lower = timeSeries - threshold;
        
        % 100 ms rule 
        % In lay terms: Make it more difficult for the upper time series to
        % cross the y=0 line by using a moving-max window. In this way, if
        % the upper time series crosses the y=0 line, we have found a
        % window completely under the y=0 line. The same principle applies
        % to the lower time series using a moving-min window. 
        % The window is defined as [0 k] before time=0 and [k 0] after
        % time=0, where k is 100 ms
        timeWin_sec = 0.100;                    % in s
        timeWin_pnt = timeWin_sec * EOG_srate;  % in pnt
        timeSeries_upper_pre  = movmax(timeSeries_upper, [0 timeWin_pnt], 'Endpoints','shrink');
        timeSeries_upper_post = movmax(timeSeries_upper, [timeWin_pnt 0], 'Endpoints','shrink');
        timeSeries_upper = [timeSeries_upper_pre(timeVector<=0); timeSeries_upper_post(timeVector>0)];
        timeSeries_lower_pre  = movmin(timeSeries_lower,  [0 timeWin_pnt], 'Endpoints','shrink');
        timeSeries_lower_post = movmin(timeSeries_lower,  [timeWin_pnt 0], 'Endpoints','shrink');
        timeSeries_lower = [timeSeries_lower_pre(timeVector<=0); timeSeries_lower_post(timeVector>0)];
        
        % QEonset find the nearest value within the threshold, in the time preceding movement initiation
        qe_onset_pnt = find([1:length(timeVector)]'<=start_pnt & (timeSeries_lower>0 | timeSeries_upper<0),1,'last');
        qe_onset_pnt = qe_onset_pnt + 1; % +1 because above we found the first value exceeding the threshold but we want the last *within* the threshold
        if isempty(qe_onset_pnt)
            qe_onset_pnt = 1; % beginning of the epoch
        end
        qe_onset_sec = timeVector(qe_onset_pnt);
        
        % QEoffset find the nearest value within the threshold, in the time following movement initiation
        qe_offset_pnt = find([1:length(timeVector)]'>=start_pnt & (timeSeries_lower>0 | timeSeries_upper<0),1,'first');
        qe_offset_pnt = qe_offset_pnt - 1; % -1 because above we found the first value exceeding the threshold but we want the last *within* the threshold
        if isempty(qe_offset_pnt)
            qe_offset_pnt = length(timeVector); % end of the epoch
        end
        qe_offset_sec = timeVector(qe_offset_pnt);
        
        % sanity check: figure
        if strcmp(doYouWantTheFigure,"yes")
            figure(10); clf;
            p1 = plot(timeVector,timeSeries, 'Color', [0 0 0 1]); hold on
            p2 = plot(timeVector,timeSeries_upper, 'Color', [1 0 0 0.5]);
            p3 = plot(timeVector,timeSeries_lower, 'Color', [0 0 1 0.5]);
            line([0 0], get(gca,'YLim'),'Color', [0 0 0 0.1])
            line(get(gca,'XLim'), [0 0],'Color', [0 0 0 0.1])
            ylabel(['EOG']);
            xlabel('Time (s)');
            set(gca,'XLim',[timeVector(1) timeVector(end)]);
            l1 = line(repmat(timeVector(qe_onset_pnt),1,2),get(gca,'YLim'), 'Color',[1 1 0 1],'LineStyle','-');
            l2 = line(repmat(timeVector(qe_offset_pnt),1,2),get(gca,'YLim'),'Color',[0 1 1 1],'LineStyle','-');
            legend([l1 l2],["QEonset" "QEoffset"])
        end
        
        
        
    case 'velocity'

        % apply Savitzky-Golay differentiation filter
        p = 1; % derivative number (1=first derivative for velocity)
        SGframelen = algorithmChoice.winlen;
        SGorder    = algorithmChoice.polynDeg;
        [~,g] = sgolay(SGorder,SGframelen);
        timeSeries_D = conv(timeSeries, factorial(p)/(-dt)^p * g(:,p+1), 'same');

        % center the time series to the point with time = 0s
        timeSeries_D = timeSeries_D - timeSeries_D(start_pnt);
        
        % compute upper and lower boundaries by subracting/adding the
        % threshold to the dataseries
        threshold = algorithmChoice.threshold;
        timeSeries_D_upper = timeSeries_D + threshold;
        timeSeries_D_lower = timeSeries_D - threshold;
        

        % QEonset find the nearest value within the threshold, in the time preceding movement initiation
        qe_onset_pnt = find([1:length(timeVector)]'<=start_pnt & (timeSeries_D_lower>0 | timeSeries_D_upper<0),1,'last');
        if qe_onset_pnt~=start_pnt
            qe_onset_pnt = qe_onset_pnt + 1; % +1 because we want the last value *within* the threshold
        end
        if isempty(qe_onset_pnt)
            qe_onset_pnt = 1; % beginning of the epoch
        end
        qe_onset_sec = timeVector(qe_onset_pnt);
        
        % QEoffset find the nearest value within the threshold, in the time following movement initiation
        qe_offset_pnt = find([1:length(timeVector)]'>=start_pnt & (timeSeries_D_lower>0 | timeSeries_D_upper<0),1,'first');
        if qe_offset_pnt~=start_pnt
            qe_offset_pnt = qe_offset_pnt - 1; % -1 because we want the last value *within* the threshold
        end
        if isempty(qe_offset_pnt)
            qe_offset_pnt = length(timeVector); % end of the epoch
        end
        qe_offset_sec = timeVector(qe_offset_pnt);
        
        
        % sanity check: figure
        if strcmp(doYouWantTheFigure,"yes")
            figure(10); clf;
            p2 = plot(timeVector,timeSeries_D, 'Color', [0 1 0 0.1]); hold on
            p3 = plot(timeVector,timeSeries_D_upper, 'Color', [1 0 0 0.1]);
            p4 = plot(timeVector,timeSeries_D_lower, 'Color', [0 0 1 0.1]);
            p1 = plot(timeVector,normalize(timeSeries,'range',get(gca,'YLim')), 'Color', [0 0 0 1]); hold on   % this is normalized only for qualitative evaluation of the EOG position waveform in relation to the QEonset and QEoffset lines. otherwise the two time series (position and velocity) would be on very different y-axis scales that would make it very hard to inspect visually
            line([0 0], get(gca,'YLim'),'Color', [0 0 0 0.1])
            line(get(gca,'XLim'), [0 0],'Color', [0 0 0 0.1])
            ylabel('EOG velocity');
            xlabel('Time (s)');
            set(gca,'XLim',[timeVector(1) timeVector(end)]);
            l1 = line(repmat(timeVector(qe_onset_pnt),1,2),get(gca,'YLim'), 'Color',[1 1 0 1],'LineStyle','-');
            l2 = line(repmat(timeVector(qe_offset_pnt),1,2),get(gca,'YLim'),'Color',[0 1 1 1],'LineStyle','-');
            legend([l1 l2],["QEonset" "QEoffset"])
        end
        
                
    otherwise
        error('select a valid algorithm')
end



% QE period
QEonset  = qe_onset_sec;
QEoffset = qe_offset_sec;

% sanity checks (for debugging)
if QEonset>0; error('QEonset cannot be positive (i.e., after time=0)'); end
if QEoffset<0; error('QEonset cannot be negative (i.e., before time=0)'); end







