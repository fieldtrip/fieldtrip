function test_old_filtering

% MEM 1gb
% WALLTIME 00:10:00


% This script is for testing ideas about how to correct
% for edge artifacts for different filter types.
% The plot shows the original signal in blue,
% the filtered signal without correction in red,
% and the filtered signal with correction in green.

% (C) 2010, S. Klanke

T = 2;
Fs = 500;
Fc = 10;
N = 4; % filter order

t = linspace(0, T, T*Fs);  % time axis for 2 seconds of data at 500 Hz
x = 0.2*t + 1 + 0.25*sin(2*pi*Fc*t) + 0.1*randn(size(t));

figure(1);
titles = {'Forward', 'Backward', 'FiltFilt'};

for k=1:4
  Wn = Fc/(0.5*Fs);
  
  switch k
    case 1
      [B,A] = butter(N, Wn, 'low');
      label = 'Lowpass';
    case 2
      [B,A] = butter(N, Wn, 'high');
      label = 'Highpass';
    case 3
      [B,A] = butter(N, [0.5*Wn 2*Wn]);
      label = 'Bandpass';
    case 4
      [B,A] = butter(N, [0.5*Wn, 2*Wn], 'stop');
      label = 'Notch/Stop';
  end
  
  % Calculate DC gain of the filter, so we now how much 
  % of the constant part of a signal the filter lets through.
  dcGain = polyval(B,1) / polyval(A,1);
  
  % forward filtering (_c = with correction)
  % the correction as shown here only works for 1-D signals,
  % but can easily repmat'ed for arbitrary matrices
  xf   = filter(B,A,x);
  xf_c = filter(B,A,x - x(1)) + dcGain*x(1);
  
  % backward filtering (_c = with correction)
  xb   = fliplr(filter(B,A,fliplr(x)));
  xb_c = fliplr(filter(B,A,fliplr(x) - x(end))) + dcGain*x(end);

  % two-pass filtering doesn't need correction
  x2 = filtfilt(B,A,x);
  
  subplot(4,3, (k-1)*3 + 1);
  plot(t,x,t,xf, 'r', t,xf_c, 'g');
  ylabel(label);
  
  subplot(4,3, (k-1)*3 + 2);
  plot(t,x,t,xb, 'r', t,xb_c, 'g');
  
  subplot(4,3, (k-1)*3 + 3);
  plot(t,x,t,x2, 'r');
end

for i=1:3
  subplot(4,3,i);
  title(titles{i});
end

