clear all
close all 
clc
c = 299792458; 


% required user input
uio = 'uio_vel_p_'; % file name to process 
To = 5.0d0; % truncate time from startup 
Tb = 5.0d0; % bin averaging window size

%import data exa


% import data
for ii = 1:5
    iox = importdata(strcat(uio,num2str(ii+3),'.dat'), ' ', 1); %assumed to be in t, x1(t), x2(t), ... 
    n(:,ii) = iox.data(:,2);
    x(:,ii) = iox.data(:,3);  %assumed to be informat (t np up vp wp) 
end

%size data
t = iox.data(:,1);
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
bigT  = zeros(1,I);
for ni = no:N-1
  dt = (t(ni+1) - t(ni-1))/2.0;
  %check if no particles in Ith bin
  for ii = 1:I
    if (n(ni,ii) > 0)
      bigT(ii) = bigT(ii) + dt;
      X(kk,ii) = X(kk,ii) + x(ni,ii)*dt;
      S(kk,ii) = S(kk,ii) + x(ni,ii).^2*dt;
    end
  end
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
      bigT  = zeros(1,I);
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

