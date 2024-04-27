clc; clear; close all

% theta = load('100.mat');
% theta = theta.mydata;

%n0 = 32;

color = {'#f94144','#f3722c','#f8961e','#f9c74f','#90be6d','#90be6d','#43aa8b','#4d908e','#577590'};
%%
close all
figure;

for i = 1:8

    n0 = 32;
    t_step = 600;
    fname = sprintf('%d.mat',t_step);
    theta = load(fname);

    theta = theta.mydata;
    comm_range = (i-1)*n0 + 1 : (i)*n0;


    data = theta(t_step,comm_range) .* 1i;

    data = exp(data);

    subplot(2,4,i);
    scatter(real(data),imag(data),'MarkerFaceColor',string(color(i)),'MarkerEdgeColor','none')
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    xticks([-1 0 1])
    yticks([-1 0 1])

end

%%

mu = [2 1];
Sigma = [9 3; 3 9];

x1 = -13:0.2:13;
x2 = -13:0.2:13;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

y = mvnpdf(X,mu,Sigma);
y = reshape(y,length(x2),length(x1));

surf(x1,x2,y)
% axis([-3 3 -3 3 0 0.4])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density')

%% Chimera index for different beta

beta = linspace(0,0.8,30);
load('chi.mat'); X = mydata; X = X(1:30);

plot(beta,X,'ko','MarkerSize',7,'MarkerFaceColor','k')
set(gca,'TickLabelInterpreter','latex','FontSize',17)
set(gca,'LineWidth',0.8,'Box','on')
xlabel('$\beta$','Interpreter','latex')
ylabel('$\chi$','Interpreter','latex')



