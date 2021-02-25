delta_t = 0.01;
fin_t = 1.999;
twoPi = round(2*pi,13);

addpath /Users/domenicgermano/Workspace/ChasteDom/anim/matlab
% addpath /Users/germanod/workspace/Chaste/anim/matlab

addpath /Users/domenicgermano/Workspace/results/


bending_curvature_with_time = zeros(10,200);
for j=0:9
    bending_curvature_with_time_temp = zeros(1,200);
    i=0;
    for time=0:delta_t:fin_t
        i=i+1;
        File = join(['RandomToFlat_Rep2_',num2str(j),'/angle_string_',num2str(time,'%.3f'),'000.txt']);
        f = fopen(File, 'r');
        C = textscan(f, '%f', 'Delimiter', ',');
        fclose(f);
        vals = C{1};

    %     curvature_sum_t = sum((abs(vals(2:end) - twoPi)).^1.2);
        curvature_sum_t = sum((abs(round(vals(2:end) - twoPi,12))).^1.2);

        bending_curvature_with_time(j+1,i) = curvature_sum_t;
    end
end

%%

figure;
% semilogy(0:delta_t:fin_t,40*nobending_curvature_with_time,'linewidth',2)
hold on;
for j=1:10
    semilogy(0:delta_t:fin_t,40*bending_curvature_with_time(j,:),'color',[0.8 0.8 0.8],'linewidth',1.5)
end
semilogy(0:delta_t:fin_t,40*mean(bending_curvature_with_time),'k','linewidth',3)
hold off;
% axis([0 48 10^-6 10^4])
% xticks([0:6:48])
yticks([10^-6 10^-4 10^-2 10^0 10^2 10^4])
ax = gca;
set(ax, 'YScale', 'log')
ax.FontWeight='bold';
ax.FontSize=14;
xlabel('t (hrs)','Interpreter','latex')
ylabel('$U^{{Bending}}$','Interpreter','latex')
% legend('Without Bending Force','With Bending Force','location','southeast')