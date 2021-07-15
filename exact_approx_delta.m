syms Pm(p)
% q = 0.04
% Pm = 1/((p^2+q^2)*((1-p)^2+q^2))
% % rts=roots(chebfun(Pm,[0,1]))
% % len=length(rts)
% fplot(Pm,[0,1])
% 
% hold on 
% % q = 0.045
% % Pm = 1/((p^2+q^2)*((1-p)^2+q^2))
% % fplot(Pm,[0,1])
% 
% q = 0.05
% Pm = 1/((p^2+q^2)*((1-p)^2+q^2))
% fplot(Pm,[0,1])
% 
% q = 0.06
% Pm = 1/((p^2+q^2)*((1-p)^2+q^2))
% fplot(Pm,[0,1])
% hold off
% 
% ylim([0 inf])
% legend({'q = 0.04','q = 0.05','q = 0.06'},'Location','north')
% xlabel('p') 
% ylabel('P_{RIS}') 
% grid on


%syms q 
load('9apr_delta.mat')
iter = 5;
net_gain = zeros(iter,1);
f = zeros(iter,1);
percent = zeros(iter,1);
del = zeros(iter,1);
%q = 0.1
f = q_vals(1:5)
%f(1) = 1/25
q = f(1)
% add approx equation here
Pm = 1/((p^2+q^2)*((1-p)^2+q^2));
Pmin = 1/((1/4+q^2)*(1/4+q^2));
max_1 = (1 - sqrt(1-4*q^2))/2;
Pmax = 1/((max_1^2+q^2)*((1-max_1)^2+q^2));
r_start = 0.85;
r = r_start;
r_P = r*(Pmax - Pmin) + Pmin;
root_all = vpasolve(Pm - r_P, p, [0 0.5]);
if length(root_all) >= 2
    del1 = root_all(1);
    del2 = root_all(2);
elseif length(root_all) == 1
    del1 = 0;
    del2 = root_all(1);
end
del_P = Pmax;
del(1) = del2 - del1;
net_gain(1) = int(Pm - r_P, p, del1, del2);
%net_gain(1) = del(1)*(Pmax - r_P)/2;
percent(1)= r_start;

for i=2:iter
    %f(i) = 1/(25-i)
    q = f(i)
    Pm = 1/((p^2+q^2)*((1-p)^2+q^2));
    Pmin = 1/((1/4+q^2)*(1/4+q^2));
    max_1 = (1 - sqrt(1-4*q^2))/2;
    Pmax = 1/((max_1^2+q^2)*((1-max_1)^2+q^2));
    r = r_start;
    
    for k=1:r_start*1000
        
        r_P = r*(Pmax - Pmin) + Pmin;
        root_all = vpasolve(Pm - r_P, p, [0 0.5])
        if length(root_all) >= 2
            del1 = root_all(1);
            del2 = root_all(2);
        elseif length(root_all) == 1
            del1 = 0;
            del2 = root_all(1);
        end
        del_P = Pmax;
        del(i) = del2 - del1;
        net_gain(i) = int(Pm - r_P, p, del1, del2);
        %net_gain(i) = del(i)*(Pmax - r_P)/2
        if(net_gain(i)>=net_gain(1))
            percent(i)= r;
            break                                                          
        else
            r = r_start - (k/1000);
        end
    end
end
piecew_del = [0.0764782653484922;0.0847556750839333;0.0925526445048391;0.100862977810077;0.109187161251635]
piecew_per = [0.850000000000000;0.834000000000000;0.819000000000000;0.803000000000000;0.787000000000000]
plot(f,del)
hold on
plot(f,piecew_del)
plot(f,base_vals(1:5)./25)
hold off
% plot(f,percent*100)
% hold on
% plot(f,piecew_per*100)
% plot(f,perc_vals(1:5)*100)
% hold off
% ylim([50 100])
ylim([0 0.2])
%legend({'Analytical Model - Exact','Analytical Model - piecewise','SimRIS'},'Location','southwest')
legend({'Exact','Piecewise','SimRIS'},'Location','northwest','FontSize',13)
xlabel('q','FontSize',13) 
ylabel('Delta (\Delta)','FontSize',13) 
%ylabel('Percentage (x*)') 
grid on
%matlab2tikz('Fig6.tex');
%print -depsc delta_percent_all
print -depsc delta_all