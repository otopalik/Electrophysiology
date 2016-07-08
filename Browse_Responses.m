function [ output_args ] = Browse_Reponses(ErevMat)
%BROWSE_RESPONSES plots all evoked responses from a given position, as
%stored by the ErevMat in Vm_traces vector. It is intended for the purpose
%fo detecting bad responses for exclusion from analysis.

close all;
load(ErevMat);
[m,n] = size(Vm_traces);

k = 1;
for i = 1:n 
    figure(1)
    supersubplot(1,6,5,k) ;
    t_start = UV_ind(i) - 5000;
    t_end = UV_ind(i) + 15000;
    plot(t_start:t_end, Vm(t_start: t_end));
    title(k);
    %ylim(y_axes);
    ylim([mean(Vm(t_start: t_end))-3 mean(Vm(t_start: t_end))+3]);
    length(Vm(t_start:t_end));
    Vm_traces(:,k) = Vm(t_start:t_end)';
    k = k+1;
end



end

