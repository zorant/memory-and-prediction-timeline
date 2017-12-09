% Temporal associations and predition added to the toy example of 
% scale-invariant, log-compressed memory representation.
%
% Some parts of this implementation follow from equations in Shankar & 
% Howard, Neural Computaiton 2012 and Shankar & Howard Journal of Machine 
% Learning Research, 2013. Equations in the text reffer to the 2012 paper.
%M
% Zoran Tiganj (zoran.tiganj@gmail.com), first version created on February 
% 23 2016  

clear all

%% Initialize parameters and alocate space

params.video_on = 0; % make a video of the input signal and the prediction
params.figure_on = 1; % draw the results after each time step

buff_len = 50; % number of output cells ("time cells") - T
k = 4; % order of the spatial derivative between intermediate and output layer
    %(higher value results in more narrow firing fields of time cells) 
N = buff_len+2*k; % number of cells in the intermediate layer (leaky integrators) - t
    % k edges will be ignored becase the spatial derivative will be
    % incomplete
Taustar_min = 0.005; % peak time of the first time cell
Taustar_max = 0.1; % peak time of the last time cell

% Create power-law growing Taustarlist (peak times of time cells) and
% corresponding s 
alpha = (Taustar_max/Taustar_min)^(1/buff_len)-1;
Taustarlist = Taustar_min*(1+alpha).^(-(k):(buff_len+(k)-1));
s = k./Taustarlist; % eq 2.3


% Create DerivMatrix that will be used to computed a spatial derivateve of t to
% compute T
DerivMatrix = zeros(N,N);
 for i = 2:N-1
   DerivMatrix(i,i-1) = -(s(i+1)-s(i))/(s(i)-s(i-1))/(s(i+1)-s(i-1));
   DerivMatrix(i,i) = ((s(i+1)-s(i))/(s(i)- s(i-1))-(s(i)-s(i-1))/(s(i+1)-s(i)))/(s(i+1)-s(i-1));
   DerivMatrix(i,i+1) = (s(i)-s(i-1))/(s(i+1)-s(i))/(s(i+1)-s(i-1));
 end
  
% Create time vector (tau) with length tau_len and resolution dtau
tau_len = 300; dtau = 0.001; tau = 0:dtau:tau_len*dtau-dtau;

% Create input signals composed of four different stimuli that activate 
% sequantially spaced by 10 time steps
%f(1,:) = rand(1,tau_len) > 0.98; % random input
f(1,:) = zeros(1,tau_len); 
f(1,[10,90,150,250]) = 1; % input at fixed time steps
f(2,:) = circshift(f(1,:)',10)';
f(3,:) = circshift(f(2,:)',10)';
f(4,:) = circshift(f(3,:)',10)';

% To demonstrate properties of S, comment the input definition above and
% use the one below where second input predicts third and fourth input
% equaly often making them similar in S (value in S(3,3), S(3,4), S(4,3)
% and S(4,4) will be similar.
% f(1,:) = rand(1,tau_len) > 0.98;
% f(2,:) = circshift(f(1,:)',10)';
% tmp = circshift(f(2,:)',10)';
% i = find(tmp==1);
% i = i(randperm(length(i)));
% f(3,i(1:round(length(i)/2))) = 1;
% f(4,i(round(length(i)/2)+1:end)) = 1;

I=size(f,1); % number of unique input stimuli

% Allocate space 
t=zeros(I,N,length(tau)); % leaky integrators
T=zeros(I,N,length(tau)); % time cells
M=zeros(I,I,buff_len); % 3 tensor that associates present (f) with 
    %compressed memory (T)
p=zeros(size(f)); % prediction
S=zeros(size(f,1),size(f,1)); % matrix that assosiates what is predicted 
    %(p) and what reapply happened (f) - 'semantic similarity'

% define variables for video encoding
if params.video_on
    v = VideoWriter('scale_free_prediction.mp4','MPEG-4');
    v.Quality = 100;
    open(v);
end
dur = 100;

figure

%% Run the main loop 
for tau_index = 2:tau_len
    for i=1:I % construct seprate memory representation for each stimulus
        t(i,:,tau_index)=t(i,:,tau_index-1)+(-s.*squeeze(t(i,:,tau_index-1))+f(i,tau_index))*dtau; % eq 2.1
        t_diff=DerivMatrix^k*t(i,:,tau_index)'; % perform k-th order derivative 
            % of t with respect to s
        T(i,:,tau_index)=(-1)^k*s.^(k+1)'.*t_diff/factorial(k)'; % eq 2.3
    end
    for n=1:buff_len % 
        M(:,:,n)=M(:,:,n)+f(:,tau_index)*T(:,n+k,tau_index)'; % eq 2.11 update 
            %assocaition between f and T at every lag n
        S=S+p(:,tau_index-1)*f(:,tau_index)'; % update association between
            % p and f
        p(:,tau_index)=p(:,tau_index)+M(:,:,n)*T(:,n+k,tau_index); % eq 2.12 probe
            % M by T to compute the ovelap between current temporal context
            % (T) to the average temporal context (M)
    end
    
    % p can be normalized to represent probability of a stimulus happening
%     if max(p(:,tau_index)) > 0
%         p(:,tau_index) = p(:,tau_index)./sum(p(:,tau_index)); % normalize p
%     end

    % plot the result at current time step
    plot(f','*','MarkerSize',14)
    ax = gca;
    ax.ColorOrderIndex = 1;
    hold on, plot(p')
    hold off
    set(gca, 'xlim', [tau_index-dur,tau_index])
    set(gca,'xtick',[tau_index-dur:10:tau_index])
    set(gca,'fontsize', 14)
    box off
    set(gca,'YColor', 'w')
    set(gca,'ylim',[0.001,1.1])
    xlabel('Time step')
    %legend('Input 1', 'Input 2', 'Input 3', 'Input 4', 'Prediction 1',...
    %     'Prediction 2', 'Prediction 3',...
    %      'Prediction 4','Orientation','horizontal','Location','bestoutside')
    set(gcf,'color','w')
    if params.figure_on
        drawnow;
    end
    if params.video_on
        % write to the video file
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end

if params.video_on
    % close the video file
    close(v);
end

% to plot the entire time series and the prediction
% figure,plot(f')
% ax = gca;
% ax.ColorOrderIndex = 1;
% hold on, plot(p')


figure
for kk=1:4
    %subplot(2,2,kk),surf(squeeze(M(kk,:,:)))
    subplot(2,2,kk),image(squeeze(M(kk,:,:)*128))
    title(sprintf('Input stimulus %i',kk))
    xlabel('\tau^*')
    ylabel('Stimulus in the past')
    set(gca,'xtick',[1,13,24,36,47])
    set(gca,'xticklabel',{-[5,10,20,40,80]})
end
colormap jet
suptitle('Temrapol associations stored in the 3-tensor M') 
set(gcf,'color','w')

% save the figure
export_fig 'M_heatmap.png' -r1200
saveas(gcf,'M_heatmap.fig')
