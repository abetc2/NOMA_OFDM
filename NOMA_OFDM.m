%% %%% Description
%%%%%% In this code one communication user and multiple sensing users send
%%%%%% their signals in the same frequency band. The algorithm used in this
%%%%%% code allows to separate signals at the receiver. For details you can
%%%%%% the paper "CSI-based Sensing with NOMA of Multiple Sensing Users for
%%%%%% ISAC". Also you can reach the paper from the following link:
%%%%%% https://ieeexplore.ieee.org/document/10464774
%% %%% Abbreviations
%%%%%% SU used for sensing user, CU used for communication user
%% %%% Analyzed Conditions
%%%%%% Only ISAC-ICE is indicated in this code for the condition 1,2,3,4
%%%%%% SUs. For CO-ISAC you can add just simple othogonal condition without
%%%%%% using the proposed algorithm.

clear all; close all; clc;

ToBeAnalyzedSUNo = 4; % Number of analyzed SUs.

for SUNo = 1:ToBeAnalyzedSUNo  %running the code for 4 number of SU condition
    
    %% Defining parameters

    num_sensing_user = SUNo;          % number of sensing users
    nFFT       = 256;                 % fft size
    nSym       = 1*10^3;              % number of symbols
    CP         = 16;                  % Cyclic Prefix length
    M          = 2;                   % 2 for BPSK and 4 for QPSK
    EbN0dB     = 0:5:30;              % bit to noise ratios for testing conditions
    nBitPerSym = nFFT*log2(M);        % data bit number per symbol
    nTap       = 8;                   % number of channel taps for CU and SUs
    % defining sensing pilots for all sensing users
    sF = zeros(nSym, nFFT, num_sensing_user);
    for i = 1:num_sensing_user
        sF(: , i:2*num_sensing_user:nFFT , i) = 2;
        sF(: , i+num_sensing_user:2*num_sensing_user:nFFT , i) = -2;
    end
    % defining the size of vectors that are used in the for loop
    xF       = zeros(nSym,nFFT);
    mse_sens = zeros(numel(EbN0dB),1);
    nErr     = zeros(numel(EbN0dB),1);
    st1      = zeros(nSym , nFFT , num_sensing_user);
    st       = zeros(nSym , nFFT+CP , num_sensing_user);
    hts      = zeros(nSym, nTap, num_sensing_user);
    hFs      = zeros(nSym, nFFT, num_sensing_user);
    sht      = zeros(nSym,nTap+nFFT+CP-1,num_sensing_user);
    hF_sens  = zeros(nSym,nFFT);
    hF_ofdm  = zeros(nSym,nFFT);
    sF_est   = zeros(nSym, floor(nFFT/num_sensing_user), num_sensing_user);
    s_t      = zeros(nSym, nFFT, num_sensing_user);
    hFs_est  = zeros(nSym, nFFT, num_sensing_user);
    dd_yF    = zeros(nSym,nFFT);
    sF_est2  = zeros(nSym, nFFT, num_sensing_user);

    %% Begining of the Algorithm

    for ii = 1:numel(EbN0dB)

        %% CU Transmitter

        ipBit = 1.*(rand(1,log2(M)*nBitPerSym*nSym) > 0.5); % random 1's and 0's for data
        ipMod = qammod(ipBit.',M,'InputType','bit','UnitAveragePower',true); % Modulation
        ipMod = reshape(ipMod,nBitPerSym,nSym).'; % grouping into multiple symbols
        xF = ipMod;  % CU signal preparing in freq
        xt1 = (nFFT/sqrt(nFFT))*ifft(xF.').'; % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1
        xt = [xt1(:,[nFFT-CP+1:nFFT]) xt1]; % Appending cylic prefix

        %% SUs Transmitters

        for i=1:num_sensing_user   % sensing pilots in time for each SU
            st1(:,:,i) = (nFFT/sqrt(nFFT))*ifft(sF(:,:,i) .').';
        end
        for i=1:num_sensing_user   % Appending cylic prefix for sensing for each SU
            st(:,:,i) = [st1(:,[nFFT-CP+1:nFFT],i) st1(:,:,i)];
        end

        %% multipath channel

        % communication user channel
        ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + 1i*randn(nSym,nTap));  % channel taps in time
        hF = (fft(ht,nFFT,2)); % computing and storing the frequency response of the channel, for use at recevier (assuming CU channel is known)
        xht=zeros(nSym,nTap+nFFT+CP-1); % the CU signal after channel
        for jj = 1:nSym
            xht(jj,:) = conv(ht(jj,:),xt(jj,:)); % convoluting channel and CU signal
        end
        xt = xht;
        % multipath channel for sensing
        for i=1:num_sensing_user
            hts(:,:,i) = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + 1i*randn(nSym,nTap)); % channel taps in time
            hFs(:,:,i) = fft(hts(:,:,i),nFFT,2); % computing and storing the frequency response of the channel, for computing SU channel estimation error
            for jj = 1:nSym
                sht(jj,:,i) = conv(hts(jj,:,i),st(jj,:,i)); % convoluting channel and SU signal
            end
        end
        xt = reshape(xt.',1,nSym*((nFFT+CP)+nTap-1)); % Concatenating multiple symbols to form a long vector
        sht2 = zeros(nSym,((nFFT+CP)+nTap-1));
        for i=1:num_sensing_user  % summing all sensing users
            sht2 = sht2 + sht(:,:,i);
        end
        sht3 = reshape(sht2.',1,nSym*((nFFT+CP)+nTap-1)); % Concatenating multiple symbols to form a long vector
        % Gaussian noise
        snr = EbN0dB(ii);
        sigma = sqrt(1/(2*(10^(snr/10)))); % Eb/N0
        nt = (sigma).*[randn(1,numel(xt)) + 1i*randn(1,numel(xt))]; % noise in time domain
        % Adding all signals,
        yt = xt + sht3 + nt;

        %% Receiver

        yt = reshape(yt.',((nFFT+CP)+nTap-1),nSym).'; % formatting the received vector into symbols
        yt = yt(:,[CP+1:(nFFT+CP)]); % removing cyclic prefix
        yF = (sqrt(nFFT)/nFFT)*(fft(yt.')).'; % converting to frequency domain
        yFR = yF;
        sFR = yFR;
        L=2;
        % Beginning of the proposed algorithm
        while L>0
            yF_ofdm = yFR;
            for i=1:num_sensing_user
                sF_est(:,:,i) = sFR(:,i:num_sensing_user:(nFFT-mod(nFFT,num_sensing_user)))./sF(:,i:num_sensing_user:(nFFT-mod(nFFT,num_sensing_user)),i); % LS estimation
                sF_est2(:,:,i) = spline(i:num_sensing_user:(nFFT-mod(nFFT,num_sensing_user)), sF_est(:,:,i), 1:nFFT);
                s_t(:,:,i) = (nFFT/sqrt(nFFT))*ifft(sF_est2(:,:,i).').'; % DFT based estimation
                hFs_est(:,:,i) = (sqrt(nFFT)/nFFT)*(fft(s_t(:,1:CP,i).', nFFT)).'; % DFT based estimation
                yF_ofdm = yF_ofdm - hFs_est(:,:,i).*sF(:,:,i); % extracting signal of OFDM
            end
            yF_eq = yF_ofdm./hF; % equalization and obtaining modulated data and pilots
            yMod = yF_eq;
            ipBitHat = qamdemod(yMod.',M,'OutputType','bit','UnitAveragePower',true).'; % to obtain data do demodulation
            dd_yF = qammod(ipBitHat.',M,'InputType','bit','UnitAveragePower',true).'; % modulation of data
            yF_itr = yFR - dd_yF.*hF; % deleting OFDM signal
            %       sFR = yF_itr;
            L=L-1;
        end %end of the proposed algorithm
        ipBitR = ipBitHat;
        ipBitR = reshape(ipBitR.',1,nSym*(1)*nFFT);
        hF_ofdm = hF;
        mse_sens(ii) = mean(mean(abs(hFs_est(:,:,1)-hFs(:,:,1)).^2)); % calculating mse error of sensing channel
        nErr(ii) = size(find(ipBitR - ipBit),2); % counting the errors
    end
    mse_sens2(SUNo,:) = mse_sens;
    simBer(SUNo,:) = nErr/(nSym*nBitPerSym);
end

figure();
semilogy(EbN0dB,mse_sens2(1,:),'o-','LineWidth',2);
hold on; grid on;
semilogy(EbN0dB,mse_sens2(2,:),'o-','LineWidth',2);
semilogy(EbN0dB,mse_sens2(3,:),'o-','LineWidth',2);
semilogy(EbN0dB,mse_sens2(4,:),'o-','LineWidth',2);
legend("2 Iteration, 1 SU","2 Iteration, 2 SU","2 Iteration, 3 SU","2 Iteration, 4 SU");
title("Mean Square Error for Channel (SU)");
ylabel("MSE"), xlabel("SNR (dB)");

figure();
semilogy(EbN0dB,simBer(1,:),'o-','LineWidth',2);
hold on; grid on;
semilogy(EbN0dB,simBer(2,:),'o-','LineWidth',2);
semilogy(EbN0dB,simBer(3,:),'o-','LineWidth',2);
semilogy(EbN0dB,simBer(4,:),'o-','LineWidth',2);
legend("2 Iteration, 1 SU","2 Iteration, 2 SU","2 Iteration, 3 SU","2 Iteration, 4 SU");
title("Bit Error Rate for CU");
ylabel("BER"), xlabel("SNR (dB)");

figure();
semilogy(EbN0dB,simBer(1,:),'o-','LineWidth',2);
hold on; grid on;
semilogy(EbN0dB,simBer(2,:),'o-','LineWidth',2);
semilogy(EbN0dB,simBer(3,:),'o-','LineWidth',2);
semilogy(EbN0dB,simBer(4,:),'o-','LineWidth',2);
legend("2 Iteration, 1 SU","2 Iteration, 2 SU","2 Iteration, 3 SU","2 Iteration, 4 SU");
title("Bit Error Rate for CU");
ylabel("BER"), xlabel("SNR (dB)");
