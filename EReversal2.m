function [Vm_traces, Vm_vect, delta_V_vect, E_rev] = EReversal2(filename, position_no, lin_reg_range, file_path)
% EReversal2 plots and fits a Delta_V vs. Vm curve for determination of the 
% reversal potential of an evoked response (i.e. iontophoresis, photo-uncaging, etc). 

close all;
k = 1;
%Load relevant files from path
abf = LoadAbf(filename);
UV = abf.data.Qswitch;
Vm = abf.data.Vm_1;

% plot the entire recording file
figure(1)
subplot(2,1,1)
plot(Vm)
ylabel('mV')
ylim([-100 -35])
title(position_no)
hold on
subplot(2,1,2)
ylabel('arb units')
xlabel('time (ms)')
plot(UV)
ylim([0 5])
xlim([0 length(Vm)])
hold on

% Identify stimulus indices 
UV_on = zeros(length(UV),1); % vector will indicate all 0/1 whether UV off/on
UV_start = zeros(length(UV),1); % vector will indicate 0/1 start of UV stimuli (trigger for Vm Responses)
UV_raw = find(UV > 1); % vector indicates raw Qswitch data

for i = 1:length(UV_raw);
    UV_on(UV_raw(i)) = 1; % vector wherein 0/1 indicates all times in which UV was off/on
end

UV_ind = find(UV_on > 0);

% This for loop will generate UV_start, which is binary 0/1 trigger for
% start of UV stimulus:
for i = 2:length(UV_ind) 
    if UV_ind(i) - UV_ind(i-1) > 11
        UV_start(UV_ind(i)) = 1;
    %else UV_start(i) = 0; % UV_start indicates only start times (trigger TTL) of UV stim
    end
end
UV_start(UV_ind(1)) = 1; % detect first stimulus, not generated in above for loop.

UV_ind = find(UV_start > 0 ) % TTL indices

% plot on figure 1:
figure(1)
subplot(2,1,1)
plot(UV_ind, -39, 'ro')
subplot(2,1,2)
plot(UV_ind, 4.5, 'ro')

% detect responses 
% plot single Vm responses separately
Vm_traces = zeros(20001,3); 


for i = 1:length(UV_ind) % **** CHANGE THESE INDICES TO INDICATE WHICH PULSES
    figure(2)
    title(position_no)
    supersubplot(2,4,4,k) %(2,round(length(UV_ind)/4),4,k) %(2, 4, 4, k) 
    t_start = UV_ind(i) - 5000;
    t_end = UV_ind(i) + 15000;
    plot(t_start:t_end, Vm(t_start: t_end));
    %ylim(y_axes);
    ylim([mean(Vm(t_start: t_end))-3 mean(Vm(t_start: t_end))+3]);
    length(Vm(t_start:t_end));
    Vm_traces(:,k) = Vm(t_start:t_end)';
    k = k+1;
end

% Calculate Vm (baseline) and deltaV for each response:
Vm_vect = zeros(k-1,1);
delta_V_vect = zeros(k-1,1);

for i = 1:k-1
    Vm_vect(i) = mean(Vm_traces(1:4000, i));
    delta_V_vect(i) = mean(Vm_traces(5200:6000, i))-Vm_vect(i);
end

figure
plot(Vm_vect(lin_reg_range), delta_V_vect(lin_reg_range), 'o', 'LineWidth', 3, =)
box off
line([-120 -30], [0 0])
xlabel('Vm (mV)')
ylabel('deltaV (mV)')
title(num2str([position_no]))
ylim([-3 3])

%Linear Regression!

x = Vm_vect(lin_reg_range);
y = delta_V_vect(lin_reg_range);
[R,p] = corrcoef(x,y);
C = cov(x,y);
m = R(1,2)*sqrt(C(2,2)/C(1,1));
b = mean(y) - m*mean(x);
E_rev = -b/m;

X = [-120 -30]; % chosen to match our X-axis of interest
Y = m * X + b;
hold on
plot(X,Y,'g--', 'LineWidth', 3);

Y_fit = m * x + b;
error = y-Y_fit;
squarederror = error.*error;% OR, error.^2 is el-by-el raising to 2nd power
meanerr = mean(squarederror);
text(-70, 2.5, 'MSE = ')
text(-50,2.5, num2str([meanerr]))
text(-70, 2, 'R = ')
text(-50, 2, num2str([R(2,1)]))
text(-70, 1.5, 'p = ')
text(-50, 1.5, num2str([p(2,1)]))
text(-70, 1, 'n = ')
text(-50, 1, num2str([length(Vm_vect)]))
text(-70, 0.5, 'x-intercept = Erev = ')
text(-40, 0.5, num2str([-b/m]))

save(file_path)

end

