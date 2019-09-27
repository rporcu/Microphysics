clear all
close all 
clc

% required user input
uio = 'uio_vel_p_0.dat'; % file name to process 
To = 5.0d0; % truncate time from startup 
Tb = 5.0d0; % bin averaging window size

%import data exa
x = importdata(uio, ' ', 1); %assumed to be in t, x1(t), x2(t), ... 

%size data
t = x.data(:,1);
x = x.data(:,3);
N = size(x,1);
I = size(x,2);

%truncate startup period
for ni = 1:N
  if (t(ni) > To)
    no = ni;
    break;
  end
end

%time average data into non-overlapping bins
K = round((t(N) - To)/Tb)
X = zeros(K,I);
S = zeros(K,I);

kk    = 1;
Ts    = To + Tb;
bigT  = 0.0d0;
for ni = no:N-1
  dt = (t(ni+1) - t(ni-1))/2.0;
  %if x(ni) .ne. 0  <-- if holes
  bigT = bigT + dt;
  X(kk,:) = X(kk,:) + x(ni,:)*dt;
  S(kk,:) = S(kk,:) + x(ni,:).^2*dt;
  %check if next time is in next non-overlapping bin 
  if (t(ni+1) >= Ts || ni == N-1)
    X(kk,:) = X(kk,:)./bigT;
    S(kk,:) = S(kk,:)./bigT;
    S(kk,:) = sqrt(S(kk,:) - X(kk,:).^2);
    if (kk == K)
      break;
    else
      kk    = kk + 1;
      Ts    = Ts + Tb;
      bigT  = 0.0d0;
    end
  end
end

% calc NOBM time averaging error 
Fs = tinv(1-0.05/2,K-1);
muhatX = mean(X);
sighatX = std(X).*Fs/sqrt(double(K));
muhatS = mean(S);
sighatS = std(S).*Fs/sqrt(double(K));

% print
fprintf(1,'\n\nx^(i)\t ta mean\t   95%% CI \t ta std\t   95%% CI\n');
for ii = 1:I
  fprintf(1,'%12.6e\t%12.6e\t%12.6e\t%12.6e\n',muhatX(ii),sighatX(ii),muhatS(ii),sighatS(ii));
end

