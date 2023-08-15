function gen_gfsk( varargin )
% Sim parameters
close all

SAMP    = 64;                           % Number of samples in the output data
ITS     = 100;                          % Number of iterations per SNR point
SNRW    = [40 30 20 10];                % Set of SNRs of signals 
FOFS    = [-200 -100 0 100 200]*1e3;    % Set of Frequency Offsets between TX and RX

% For SNRW and FOFS I would use all of them for training and when testing
% test for one impairement at a time and generate two plots
%  1. Detection Rate of each vs SNR with 0 offset
%  2. Detection Rate of each vs Frequency Offset at 10 dB SNR
% With 0 offset I would also plot detection performance of each modualtion
% vs SAMP 64, 32, 16, 8 

MI      = 0.5;     % Transmit modulation index.
NSYMS   = 32;      % Number of symbols in the simulation
OSR     = 4;       % Oversample rate for normal radio part of sim
RXBW    = 1.2;     % 6dB RX symbol bandwidth filter in MHz.  Usually fixed.
BT1     = 0.5;     % BT Gauss-filter-1 for wanted communications
Fs      = 1e6;     % Symbol rate
PE      = 0;       % Plot enable

fid1 = fopen('data_iq.csv','w');
fid2 = fopen('data_ph.csv','w');

spSequence = [...
    0 0 0 0 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0
    1 1 0 0 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0
    2 0 1 0 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0
    3 1 1 0 0 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1
    4 0 0 1 0 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1
    5 1 0 1 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0
    6 0 1 1 0 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1
    7 1 1 1 0 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1
    8 0 0 0 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1
    9 1 0 0 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1
    10 0 1 0 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1
    11 1 1 0 1 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0
    12 0 0 1 1 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0
    13 1 0 1 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1
    14 0 1 1 1 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0
    15 1 1 1 1 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0];

spSequence = spSequence(:, 6:end);  % strip off index, keep only spreading seq


Fso   = Fs*OSR;
Brx   = sinc( (-OSR*8/2:OSR*8/2)*RXBW/OSR )/OSR.*hanning(OSR*8+1).'*RXBW; % LP filter 6dB bw
bg1   = gaussfir( BT1, 3, OSR );  % GFSK filter for transmit

for sn = 1:length(SNRW) % For each SNR
    
    snr = SNRW(sn);
    
    for it = 1:ITS % For each iteration to build statistics
        
        for frq = 1:length(FOFS)
                    nfrq = FOFS(frq)/Fso;
            
        % Generate BLE Signal - copied from standard
        syms = randi([0 1], 1, NSYMS); % Construct random symbol sequence from wanted transmitter
        wf   = [zeros(1, 1.5*OSR) reshape( repmat((syms*2-1), OSR, 1), 1, [] ) zeros(1, 0.5*OSR)]*(Fs/4); % FSK
        fwf  = filter(bg1, 1, [wf zeros(1,(length(bg1)-1)/2)]);  % Convert FSK into GFSK and
        fwf  = fwf((length(bg1)-1)/2+1:end);                     % remove filter delay
        ph   = 2*pi*cumsum(fwf)*MI/0.5/Fso;
        wsc  = exp(1i*ph);                                         % Wanted signal (clean) GFSK modulated, with specific MI
        

        fwsc = wsc.*exp(1j*2*pi*(0:length(wsc)-1)*nfrq);
        
        wmp  = fwsc + 10^(-(snr-10*log10(OSR))/20)*[1 1j]/sqrt(2)*randn(2, length(fwsc)); % Add wanted receiver noise
        lps  = filter(Brx, 1, [wmp zeros(1,(length(Brx)-1)/2)]); % LP receiver channel filter.
        lps  = lps(OSR+(length(Brx)-1)/2+1:end);                % Calibrate out filter delay.
        dw   = imag( lps(OSR+1:end) .* conj(lps(1:end-OSR)) ); % Descriminator signal.
        off = 6+round(rand * length(dw)/4);
        
        fprintf(fid2,'%d, %d, %f, BLE, ',it, snr, nfrq);
        fprintf(fid2,'%f, ',dw(off:off+SAMP-1));
        fprintf(fid2,'\n');
        fprintf(fid1,'%d, %d, %f, BLE, ',it, snr, nfrq);
        fprintf(fid1,'%f, ',real(lps(off:off+SAMP-1)));
        fprintf(fid1,'\n');
        fprintf(fid1,'%d, %d, %f, BLE, ',it, snr, nfrq);
        fprintf(fid1,'%f, ',imag(lps(off:off+SAMP-1)));
        fprintf(fid1,'\n');

        if PE
        figure(1)
        subplot(311)
        plot(dw(off:off+SAMP-1))
        hold on
        subplot(312)
        plot(real(lps(off:off+SAMP-1)))        
        hold on
        subplot(313)
        plot(imag(lps(off:off+SAMP-1)))
        hold on
        end
        
        % Generate OQPSK Signal - copied from standard 

        %[K, L] = size( spSequence );        %
        %spreadingGain = L/log2(K);

        %Rb = 250e3;                 % information bit rate
        %Rc = Rb*spreadingGain;      % chip rate 
        numByte = 1;
        packet = randi( [0,1], 1, 8*(numByte) );	% 0/1 binary payload 
        symbolIdx = 2.^(0:3)* reshape(packet, 4, length(packet)/4 ) +1;  % symbol index \in [0:15]+1
        chip = spSequence( symbolIdx', :)';     % TX column-wise 
        chip = chip(:)';                        % flatten, TX sequentially
        chipI = chip(1:2:end)*2-1;	% In-phase   spreading sequence 1/0 -> 1/-1
        chipQ = chip(2:2:end)*2-1;	% quadrature spreading sequence 1/0 -> 1/-1
        %Fs = OSR*Rc/2;              % sampling freq of DBB:  chip rate in I/Q is 1M
        pulseFil = sin( (0:OSR-1)/OSR *pi )';    
        chipIfil = pulseFil*chipI;      chipIfil = chipIfil(:)';
        chipQfil = pulseFil*chipQ;      chipQfil = chipQfil(:)';
        sigCL = [ chipIfil zeros(1, OSR/2)] + 1i*[zeros(1, OSR/2)  chipQfil];
        
        sigFO = sigCL.*exp(1j*2*pi*(0:length(sigCL)-1)*nfrq);
        
        sigBB  = sigFO + 10^(-(snr-10*log10(OSR))/20)*[1 1j]/sqrt(2)*randn(2, length(sigFO)); % Add wanted receiver noise
        sigPH = diff(unwrap(angle(sigBB)))/pi*OSR;
        
        fprintf(fid2,'%d, %d, %f, Zigbee, ',it, snr, nfrq);
        fprintf(fid2,'%f, ',sigPH(off:off+SAMP-1));
        fprintf(fid2,'\n');
        fprintf(fid1,'%d, %d, %f, Zigbee, ',it, snr, nfrq);
        fprintf(fid1,'%f, ',real(sigBB(off:off+SAMP-1)));
        fprintf(fid1,'\n');
        fprintf(fid1,'%d, %d, %f, Zigbee, ',it, snr, nfrq);
        fprintf(fid1,'%f, ',imag(sigBB(off:off+SAMP-1)));
        fprintf(fid1,'\n');
        
        if PE
        figure(2)
        subplot(311)
        plot(sigPH(off:off+SAMP-1))
        hold on
        subplot(312)
        plot(real(sigBB(off:off+SAMP-1)))        
        hold on
        subplot(313)
        plot(imag(sigBB(off:off+SAMP-1)))
        hold on
        end
        
        
        % Generate Noise Only
        wmp  = 10^(-(10-10*log10(OSR))/20)*[1 1j]/sqrt(2)*randn(2, length(fwsc)); % Add wanted receiver noise
        nps  = filter(Brx, 1, [wmp zeros(1,(length(Brx)-1)/2)]); % LP receiver channel filter.
        nps  = nps(OSR+(length(Brx)-1)/2+1:end);                % Calibrate out filter delay.
        %dwn   = imag( nps(OSR+1:end) .* conj(nps(1:end-OSR)) ); % Descriminator signal.
        
        dwn = diff(unwrap(angle(nps)))/pi*OSR;
        
        off = 6+round(rand * length(dwn)/4);
        
        
        fprintf(fid2,'%d, %d, %f, Noise, ',it, snr, nfrq);
        fprintf(fid2,'%f, ',dwn(off:off+SAMP-1));
        fprintf(fid2,'\n');
        fprintf(fid1,'%d, %d, %f, Noise, ',it, snr, nfrq);
        fprintf(fid1,'%f, ',real(nps(off:off+SAMP-1)));
        fprintf(fid1,'\n');
        fprintf(fid1,'%d, %d, %f, Noise, ',it, snr, nfrq);
        fprintf(fid1,'%f, ',imag(nps(off:off+SAMP-1)));
        fprintf(fid1,'\n');

        if PE
        figure(3)
        subplot(311)
        plot(dwn(off:off+SAMP-1))
        hold on
        subplot(312)
        plot(real(nps(off:off+SAMP-1)))        
        hold on
        subplot(313)
        plot(imag(nps(off:off+SAMP-1)))
        hold on
        end
        
        
    end
    end
end

fclose('all');

end

