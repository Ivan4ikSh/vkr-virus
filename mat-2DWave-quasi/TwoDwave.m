clear; clc;

%% Model parameters
R0 = 1.8;                  % basic reproduction ratio
Ub = 1e-3;                % mutation rate
a = 7;                     % immunity half-distance
N = 1e10;                  % total population size
% varx = 0.01;             % больше не используется (периодические границы)

%% Internal parameters
L = 50;                    % number of variants for each coordinate
Tmax = 1000;                % max time (units of recovery time)
stept = 0.5;               % time step
M = round(Tmax / stept);   % number of time points
SEED = 1;
print_interval = 100;      % print progress every 100 steps

rng(SEED);

%% Print parameters
fprintf('Grid size: %d x %d (periodic boundaries)\n', L, L);
fprintf('Time steps: %d\n', M);
fprintf('R0 = %.1f\n', R0);
fprintf('Ub = %g\n', Ub);
fprintf('a = %d\n', a);
fprintf('dt = %.1f\n', stept);
fprintf('init_infected = %.2e\n', 1e-2);
fprintf('N = %g\n', N);
fprintf('Seed = %d\n\n', SEED);

fprintf('%5s | %8s | %8s | %6s | %7s | %10s | %9s\n', 'Step', 'finf', 'frec', 'norm', 'wave_r', 'max_I', 'diversity');
fprintf('%5s-|-%8s-|-%8s-|-%6s-|-%7s-|-%10s-|-%9s\n', repmat('-',1,5), repmat('-',1,8), repmat('-',1,8), repmat('-',1,6), repmat('-',1,7), repmat('-',1,10), repmat('-',1,9));

%% Immunity matrix with periodic distance (only depends on column index j)
K_immunity = zeros(L, L);
for j1 = 1:L
    for j2 = 1:L
        d = min(abs(j1 - j2), L - abs(j1 - j2));
        K_immunity(j1, j2) = 1 - exp(-d / a);
    end
end

%% Initial conditions
I = zeros(L, L);                % Infected (fraction)
R = zeros(L, L);                % Recovered (fraction)
I(:, 16) = 1e-2 / L;            % infected line at column 16
R(:, 1:15) = (1 - 1e-2) / 15 / L; % recovered strip columns 1-15

Inew = zeros(L, L);
norm_ = zeros(1, M);
finf_ = zeros(1, M);
frec_ = zeros(1, M);
mean_X_ = zeros(1, M);
max_I_ = zeros(1, M);
diversity_ = zeros(1, M);

% Snapshot settings
n_snapshots = 4;
snapshot_indices = round(linspace(1, M, n_snapshots));
snapshots = cell(n_snapshots, 2);
snapshot_times = zeros(1, n_snapshots);
snapshots{1,1} = I;
snapshots{1,2} = R;
snapshot_times(1) = 0;
snap_idx = 2;

%% Main loop
for k = 1:M
    t = k * stept;
    
    % Diagnostics at beginning of step
    norm_(k) = sum(sum(I + R));
    finf_(k) = sum(sum(I)) / norm_(k);
    frec_(k) = sum(sum(R)) / norm_(k);   % temporary
    if sum(sum(I)) > 0
        % mean column index (j) of infected, weighted by I
        j_idx = 1:L;
        mean_X_(k) = sum(sum(I .* j_idx)) / sum(sum(I));
    else
        mean_X_(k) = NaN;
    end
    max_I_(k) = max(max(I));
    diversity_(k) = sum(sum(I > 0));
    
    % Print progress
    if mod(k, print_interval) == 0 || k == 1
        fprintf('%5d | %8.6f | %8.6f | %6.2f | %7.3f | %10.2e | %9d\n', ...
            k, finf_(k), frec_(k), norm_(k), mean_X_(k), max_I_(k), diversity_(k));
    end
    
    % Precompute sums over first index for efficiency
    Rsum = sum(R, 1);   % 1 x L, sum over i
    Isum = sum(I, 1);   % 1 x L
    
    % Compute Q and P vectors (dependent only on j)
    Qvec = Rsum * K_immunity;   % 1 x L, Q for each j
    Pvec = Isum * K_immunity;   % 1 x L, P for each j
    
    % Update all strains
    for i = 1:L
        for j = 1:L
            Q = Qvec(j);
            P = Pvec(j);
            
            Inew(i,j) = I(i,j) * (1 + stept * (R0 * Q - 1));
            R(i,j) = R(i,j) * (1 - stept * R0 * P) + stept * I(i,j);
            
            % ----- Random genetic drift (Poisson for low counts) -----
            exp_num = Inew(i,j) * N;
            if exp_num <= 10
                xm = round(6*exp_num);
                Inew(i,j) = sum(xm*rand(1,xm) < exp_num) / N;
            end
    
            % ----- Mutation with periodic boundaries -----
            % Neighbors with periodic indices
            i_up   = mod(i-2, L) + 1;
            i_down = mod(i,   L) + 1;
            j_left = mod(j-2, L) + 1;
            j_right= mod(j,   L) + 1;
            
            in = stept * Ub * (I(i_up, j) + I(i_down, j) + I(i, j_left) + I(i, j_right));
            out= stept * Ub * 4 * I(i, j);
            
            % Drift for incoming mutations
            exp_num_in = in * N;
            if exp_num_in <= 10
                xm = round(6*exp_num_in);
                in = sum(xm*rand(1,xm) < exp_num_in) / N;
            end
            
            Inew(i,j) = Inew(i,j) + in - out;
        end
    end
    
    % Cutoff: keep only strains with at least one individual
    Inew = Inew .* ( (Inew * N) >= 1 );
    I = Inew;
    
    % Store snapshot if needed
    if snap_idx <= n_snapshots && k == snapshot_indices(snap_idx)
        snapshots{snap_idx,1} = I;
        snapshots{snap_idx,2} = R;
        snapshot_times(snap_idx) = t;
        snap_idx = snap_idx + 1;
    end
end

% Finalise recovered fraction (cumulative)
for k = 1:M
    frec_(k) = sum(sum(R)) / norm_(k);
end

%% Create output folder
if ~exist('output', 'dir')
    mkdir('output');
end

%% Figure 1: Evolution snapshots
figure('Position', [100 100 1400 700]);
for t_idx = 1:n_snapshots
    subplot(2, n_snapshots, t_idx);
    imagesc(snapshots{t_idx,1});
    colorbar;
    title(sprintf('Infected t=%.2f', snapshot_times(t_idx)));
    xlabel('Y'); ylabel('X');
    axis xy equal tight;
    
    subplot(2, n_snapshots, t_idx + n_snapshots);
    imagesc(snapshots{t_idx,2});
    colorbar;
    title(sprintf('Recovered t=%.2f', snapshot_times(t_idx)));
    xlabel('Y'); ylabel('X');
    axis xy equal tight;
end
sgtitle(sprintf('R0=%.1f, a=%d, Ub=%g, N=%g (periodic)', R0, a, Ub, N));
saveas(gcf, fullfile('output', 'evolution_snapshots.png'));

%% Figure 2: Diagnostic plots
figure('Position', [100 100 1200 400]);
subplot(1,3,1);
plot((1:M)*stept, finf_, 'r-', 'LineWidth', 2);
%hold on;
%plot((1:M)*stept, frec_, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('Fraction'); legend('Infected','Recovered');
title('Population fractions');

subplot(1,3,2);
plot((1:M)*stept, norm_, 'k-', 'LineWidth', 2);
xlabel('Time'); ylabel('Total norm'); title('Conservation check (~1)');

subplot(1,3,3);
plot((1:M)*stept, mean_X_, 'g-', 'LineWidth', 2);
xlabel('Time'); ylabel('Mean column index'); title('Wave position (periodic)');

sgtitle(sprintf('R0=%.1f, a=%d, Ub=%g, N=%g (periodic)', R0, a, Ub, N));
saveas(gcf, fullfile('output', 'diagnostics.png'));