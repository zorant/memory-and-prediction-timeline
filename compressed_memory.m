% Toy example of scale-invariant, log-compressed memory representation
%
% The implementation follows from equations in Shankar & Howard, Neural
% Computaiton 2012 and Shankar & Howard Journal of Machine Learning Research,
% 2013. Equations in the text refer to the 2012 paper.
%
% Zoran Tiganj (zoran.tiganj@gmail.com) on May 22 2015  


%% Initialize parameters and allocate space

buff_len = 50; % number of output cells ("time cells") - T
k = 10; % order of the spatial derivative between intermediate and output layer
    %(higher value results in more narrow firing fields of time cells) 
N = buff_len+2*k; % number of cells in the intermediate layer (leaky integrators) - t
    % k edges will be ignored because the spatial derivative will be
    % incomplete
Taustar_min = 1; % peak time of the first time cell
Taustar_max = 10; % peak time of the last time cell

% Create power-law growing Taustarlist (peak times of time cells) and
% corresponding s 
alpha = (Taustar_max/Taustar_min)^(1/buff_len)-1;
Taustarlist = Taustar_min*(1+alpha).^(-(k):(buff_len+(k)-1));
s = k./Taustarlist; % eq 2.3


% Create DerivMatrix that will be used to computed a spatial derivative of t to
% compute T
DerivMatrix = zeros(N,N);
 for i = 2:N-1
   DerivMatrix(i,i-1) = -(s(i+1)-s(i))/(s(i)-s(i-1))/(s(i+1)-s(i-1));
   DerivMatrix(i,i) = ((s(i+1)-s(i))/(s(i)- s(i-1))-(s(i)-s(i-1))/(s(i+1)-s(i)))/(s(i+1)-s(i-1));
   DerivMatrix(i,i+1) = (s(i)-s(i-1))/(s(i+1)-s(i))/(s(i+1)-s(i-1));
 end
  
% Create time vector (tau) with length tau_len and resolution dtau
tau_len = 10000; dtau = 0.001; tau = 0:dtau:tau_len*dtau-dtau;

% Create input signal f of length N; in this example it contains two "delta" pulses 
f = zeros(1,tau_len); f(4000) = 1; f(8000) = 1;

% Alocate space for t and T
t = zeros(N,tau_len);
T = zeros(N,tau_len);


%% Run the main loop 
% Storing t and T for every tau_index is done only for plotting purpose,
% both t and T are time local (t is computed using only its single previous value
% - we're taking first order time-derivative and T is computed as a linear
% combination of time-local values of t)
tic
for tau_index = 2:tau_len
    t(:,tau_index) = t(:,tau_index-1)+((-s'.*t(:,tau_index-1)+f(tau_index))*dtau); % eq 2.1
    t_diff = DerivMatrix^k*t(:,tau_index); % perform k-th order derivative 
        % of t with respect to s
    T(:,tau_index) = (-1)^k*s.^(k+1)'.*t_diff/factorial(k)'; % eq 2.3
end
toc



%% Plot the results

figure

subplot(3,1,1)
plot(f,'LineWidth',2);
hold on
scatter(10000-Taustarlist(k+1:end-k)*1000,T(k+1:end-k,end)/max(T(k+1:end-k,end)),'Marker','.');
legend('Input signal as a function of time',...
    'Normalized memory representation as a function of \tau* at time = 10000',...
    'Location', 'northwest')
set(gca,'YTickLabel',{})
set(gca,'YTick',[])

% Notice that the reconstruction peaks at k/(k+1) of the input peak

for i=1:5:buff_len % plot every 5th to make the figure look clean
    subplot(3,1,2)
    hold on, plot(t(i+k,:)) 
    title('Intermediate layer - leaky integrators (only 10 shown for clarity, out of 50 )')

    subplot(3,1,3)
    hold on, plot(T(i+k,:)) % T only makes sense between k+1 and end-k, 
        % since we're taking k-th spatial derivative, edges can not be
        % properly computed
    title('Output layer - time cells (only 10 shown for clarity, out of 50 )')
end

xlabel('Time')

set(gcf,'color','w')

% safe the figure
export_fig('compressed_memory.png')
%saveas(gcf,'compressed_memory.fig')
