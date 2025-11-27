%% STX_FUSFP_chronic-chip
% This code is for extraction and processing of FP data generated during FUS stimulation
% Clear workspace
clear; clc; close all;
%% Input
% fileName
fileName = 'SD_5_R5C4-250919-132833';

% Update this with your actual CSV data path
DATA_PATH = ['D:\UKF 2025\LisaR\' fileName '\Data']; % csv data from pMAT
dFF_Table = readtable([DATA_PATH,'\DeltaF_Trace .csv']);
Event_Table = readtable([DATA_PATH,'\PC0_EventData .csv']); % control first stim time - can be off if chip test was performed during recording

% Keep only rows where event is 'PC0_'-change if TTL input is not PC0
Event_Table = Event_Table(strcmp(Event_Table.Event, 'PC0_'), :);


%% read Event
Event_Table = sortrows(Event_Table, 'Onset');
% Determine recording start time from first TTL pulse
if isempty(Event_Table)
    error('No TTL events (PC0_) found in Event_Table.');
end

onset = Event_Table.Onset;
offset = Event_Table.Offset;
% Compute time gap between consecutive events
gap = onset(2:end) - offset(1:end-1);
% Define threshold for a condition switch (e.g. >40 s)
condition_breaks = find(gap > 40);

% edges: segment boundaries (use 0 for before first and height for after last)
edges = [0; condition_breaks; height(Event_Table)];
nSegments = numel(edges) - 1;

% expected layout
nSites_expected = 21;
reps_per_condition = 3;
if nSegments ~= nSites_expected * reps_per_condition
    warning('Detected %d segments, expected %d (= %d x %d).', ...
        nSegments, nSites_expected * reps_per_condition, ...
        nSites_expected, reps_per_condition);
end

% Initialize labeling columns
Event_Table.Site = zeros(height(Event_Table),1);
Event_Table.Repetition = zeros(height(Event_Table),1);

valid_seg_counter = 0; % only increments when segment is valid

for seg = 1:nSegments
    idx = (edges(seg) + 1) : edges(seg + 1);
    nEvents_seg = numel(idx);

    % Skip short/incomplete segments (?2 events)
    if nEvents_seg <= 2
        continue;
    end

    % Count only valid segments
    valid_seg_counter = valid_seg_counter + 1;

    % Assign Site (1?21) and Repetition (1?3)
    site = ceil(valid_seg_counter / reps_per_condition);
    rep  = mod(valid_seg_counter - 1, reps_per_condition) + 1;

    Event_Table.Site(idx) = site;
    Event_Table.Repetition(idx) = rep;
end
                        
%% General plot
% Ensure event table sorted by onset
Event_Table = sortrows(Event_Table, 'Onset');

% Plot continuous dFF trace
figure('Color','w'); hold on;
plot(dFF_Table.Time, dFF_Table.DeltaF_F, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('\DeltaF/F');
title('Continuous \DeltaF/F with Condition Backgrounds');
grid on;

% Get y-range for shading
ymin = min(dFF_Table.DeltaF_F);
ymax = max(dFF_Table.DeltaF_F);
yrange = ymax - ymin;
if yrange == 0
    yrange = abs(ymax)*0.1 + 1;
    ymin = ymin - yrange/2;
    ymax = ymax + yrange/2;
    yrange = ymax - ymin;
end

% Condition info
conditions = unique(Event_Table.Site);
nCond = numel(conditions);
colors = lines(nCond);  % color per condition
alphaVal = 0.15;        % transparency

% For each condition: merge all repetitions, then shade + vertical text
for ic = 1:nCond
    c = conditions(ic);
    rows = Event_Table(Event_Table.Site == c, :);
    start_t = min(rows.Onset);
    stop_t  = max(rows.Offset);
    
    % Background shaded region
    fill([start_t stop_t stop_t start_t], ...
         [ymin ymin ymax ymax], ...
         colors(ic,:), 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'Clipping', 'off');
    
    % Vertical text in the middle
    xmid = (start_t + stop_t)/2;
    text(xmid, (ymin + ymax)/2, sprintf('Cond %d', c), ...
         'Rotation', 90, ...  % vertical
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', 8, ...
         'FontWeight', 'bold', ...
         'Color', [0 0 0], ...
         'BackgroundColor', 'none', ...
         'Clipping', 'off');
end

% Redraw trace on top so it?s not hidden by shading
uistack(findobj(gca, 'Type', 'line'), 'top');

ylim([ymin ymax]);
hold off;

%% Parameters (edit if your protocol differs)
preDur  = 5;   % seconds before stim onset used for baseline (pre)
stimDur = 5;   % seconds of stimulation
postDur = 30;  % seconds after stim end for post window
minValidPreSamples = 5;  % minimal samples required in pre window

%% Group events into trials by Site+Repetition
% We will treat each unique (Site,Repetition) as one trial. Use onset = min(Onset) of rows in that group
if ~all(ismember({'Site','Repetition','Onset','Offset'}, Event_Table.Properties.VariableNames))
    error('Event_Table must contain columns: Site, Repetition, Onset, Offset.');
end

groups = findgroups(Event_Table.Site, Event_Table.Repetition);
nGroups = max(groups);

trialList = struct('Site',[],'Rep',[],'onset',[],'offset',[],'rows',[]);
trialCounter = 0;
for g = 1:nGroups
    idx = find(groups==g);
    if isempty(idx), continue; end
    % skip groups with too few events
    if numel(idx) < 1, continue; end
    % ===  skip items where Site==0 OR Repetition==0 ===
    thisSite = Event_Table.Site(idx(1));
    thisRep  = Event_Table.Repetition(idx(1));
    if thisSite == 0 || thisRep == 0
        continue;   % skip this group entirely
    end
    
    trialCounter = trialCounter + 1;
    trialList(trialCounter).Site = Event_Table.Site(idx(1));
    trialList(trialCounter).Rep  = Event_Table.Repetition(idx(1));
    trialList(trialCounter).onset  = min(Event_Table.Onset(idx));
    trialList(trialCounter).offset = max(Event_Table.Offset(idx));
    trialList(trialCounter).rows   = idx;
end
nTrials = numel(trialList);
fprintf('Found %d trials (unique Site+Rep groups).\n', nTrials);

%% Build aligned trial matrices
% Pre-compute nominal sample counts
Fs = 1000; % Hz
preN  = round(preDur  * Fs);
stimN = round(stimDur * Fs);
postN = round(postDur * Fs);
totalN = preN + stimN + postN;

% relative time vector (centered on stim onset)
t_rel = ((-preN):(stimN+postN-1))' / Fs; % column vector length totalN

% allocate: trials x repetitions x time

DFFTrials = nan(nTrials, numel(t_rel));
ZTrials    = nan(nTrials, numel(t_rel));
TrialInfo = table((1:nTrials)', [trialList.Site]', [trialList.Rep]', ...
    [trialList.onset]', [trialList.offset]', 'VariableNames', ...
    {'TrialIdx','Site','Rep','Onset','Offset'});

% For overlap checks: collect all trial onsets for later comparison
allOnsets = [trialList.onset];
%% === Prepare DFFTrials (trial x time) ===
t = dFF_Table.Time;
dFF = dFF_Table.DeltaF_F;

for i = 1:nTrials
    % trial information
    on = trialList(i).onset;
    rep = trialList(i).Rep;
    
    % convert durations to number of samples
    preN  = round(preDur  * Fs);
    stimN = round(stimDur * Fs);
    postN = round(postDur * Fs);
    
    % idx in the signal vector
    [~, idx_on] = min(abs(t - on));       % index of the closest time point
    idx_pre_start = idx_on - preN;                % start of pre window
    idx_stim_end  = idx_on + stimN - 1;          % end of stim window
    idx_post_end  = idx_on + stimN + postN - 1;  % end of post window
    
    
    % extract segment
    seg = dFF(idx_pre_start : idx_post_end);  % raw ?F/F for pre+stim+post
    DFFTrials(i, :) = seg;
    
end
fprintf('Finished preparing raw trial-aligned DFFTrials.\n');

%% === Prepare ZTrials (trial x  time) ===
% hard coding
base_idx = find(t_rel >= -2 & t_rel <= -1);

for i = 1:nTrials
    seg = DFFTrials(i,:);   %
    % baseline = [-2 -1] s window
    base = seg(base_idx);
    
    mu = mean(base, 'omitnan');
    sd = std(base, 0, 'omitnan');
    
    if sd == 0 || isnan(sd)
        warning('Trial %d rep %d: baseline SD invalid. Skipped.', i, r);
        continue;
    end
    
    % Z-score
    ZTrials(i, :) = (seg - mu) ./ sd;
    
end
fprintf('Finished computing ZTrials using [-2 -1]s baseline.\n');


%% plot separately
plotSite = 12;     % <---- YOU SET THIS
idx = find(TrialInfo.Site == plotSite);
if isempty(idx)
    error('No trials found for Site %d', plotSite);
end

% Combine the 3 repetitions for each trial (flatten rep dimension)
dffData = DFFTrials(idx, :);    % (trials x time)
zData   = ZTrials(idx, :);

% Compute mean & SEM
dffMean = mean(dffData, 1, 'omitnan');
dffSEM  = std( dffData, 0, 1, 'omitnan') ./ sqrt(size(dffData,1));

zMean = mean(zData, 1, 'omitnan');
zSEM  = std( zData, 0, 1, 'omitnan') ./ sqrt(size(zData,1));

    % ----------- PLOT ?F/F -----------
    figure; hold on;
    title(sprintf('Site %d - ?F/F', plotSite));

    % individual grey traces
    plot(t_rel, dffData', 'Color', [0.7 0.7 0.7 0.4]);

    % SEM shading
    fill([t_rel; flipud(t_rel)], ...
        [dffMean + dffSEM, fliplr(dffMean - dffSEM)], ...
        [0.4 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % mean line
    plot(t_rel, dffMean, 'Color', [0 0.2 0.8], 'LineWidth', 2);

    xlabel('Time (s)');
    ylabel('?F/F');
    xline(0, '--r', 'Stim On');
    grid on; box on;


    % ----------- PLOT Z-Score -----------
    figure; hold on;
    title(sprintf('Site %d - Z-Score', plotSite));

    % individual grey traces
    plot(t_rel, zData', 'Color', [0.7 0.7 0.7 0.4]);

    % SEM shading
    fill([t_rel; flipud(t_rel)], ...
        [zMean + zSEM, fliplr(zMean - zSEM)], ...
        [0.4 0.9 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % mean line
    plot(t_rel, zMean, 'Color', [0 0.6 0], 'LineWidth', 2);

    xlabel('Time (s)');
    ylabel('Z-Score');
    xline(0, '--r', 'Stim On');
    grid on; box on;

%% pause
%% all plots
sitesToPlot = unique(TrialInfo.Site);   % or manually define, e.g. [6 12 15]
nSites = numel(sitesToPlot);
% dff
figure; 
for i = 1:nSites

    plotSite = sitesToPlot(i);

    % --- Find trials for this site ---
    idx = find(TrialInfo.Site == plotSite);
    if isempty(idx)
        warning('No trials found for Site %d', plotSite);
        continue;
    end

    % --- Extract data ---
    dffData = DFFTrials(idx, :);   % (trials x time)
    dffMean = mean(dffData, 1, 'omitnan');
    dffSEM  = std(dffData, 0, 1, 'omitnan') ./ sqrt(size(dffData,1));

    % ---- Create subplot ----
    subplot(ceil(nSites/3), 3, i);  % 3 columns layout
    hold on; box on;

    % individual grey traces
    plot(t_rel, dffData', 'Color', [0.7 0.7 0.7 0.4]);

    % SEM shading
    fill([t_rel; flipud(t_rel)], ...
         [dffMean + dffSEM, fliplr(dffMean - dffSEM)], ...
         [0.4 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % mean line
    plot(t_rel, dffMean, 'Color', [0 0.2 0.8], 'LineWidth', 2);

    title(sprintf('Site %d', plotSite));
    xlabel('Time (s)');
    ylabel('?F/F');
    xline(0, '--r');
end
sgtitle('?F/F for All Sites');


% Zscore
figure;
for i = 1:nSites
    
    plotSite = sitesToPlot(i);
    
    % --- Find trials for this site ---
    idx = find(TrialInfo.Site == plotSite);
    if isempty(idx)
        warning('No trials found for Site %d', plotSite);
        continue;
    end
    
    % --- Extract data ---
    zData = ZTrials(idx, :);
    zMean = mean(zData, 1, 'omitnan');
    zSEM  = std(zData, 0, 1, 'omitnan') ./ sqrt(size(zData,1));
    
    
    % ---- Create subplot ----
    subplot(ceil(nSites/3), 3, i);  % 3 columns layout
    hold on; box on;
    
    % individual grey traces
    plot(t_rel, zData', 'Color', [0.7 0.7 0.7 0.4]);
    
    % SEM shading
    fill([t_rel; flipud(t_rel)], ...
        [zMean + zSEM, fliplr(zMean - zSEM)], ...
        [0.4 0.9 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % mean line
    plot(t_rel, zMean, 'Color', [0 0.6 0], 'LineWidth', 2);
    
    title(sprintf('Site %d', plotSite));
    xlabel('Time (s)');
    ylabel('ZScore');
    xline(0, '--r');
end
sgtitle('Zscore for All Sites');
