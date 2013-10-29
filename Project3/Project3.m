%{
EE-2347 Project 3
Contributers: Brian McRee
10/13/2013

Published open-source on Git: LegallyKF5RCL/MatlabProject3
%}

clc;
clear all;
close all;

%MessageData = '5';
%save('Message.txt', MessageData);
FileP = fopen('Message.txt');
Message = fscanf(FileP, '%s');
disp(Message);

BinaryCharacters = dec2bin(Message);    %convert to binary
StringLen = length(BinaryCharacters);   %find how many characters we have
CharLen = 7;                            %length of char will always be 7
TotalLen = StringLen * CharLen;

%make a single array for the message
BinaryMessage = zeros([1,(TotalLen)]);
BinaryMessage = char(BinaryMessage);

for i = 1:StringLen
    for j = 1:CharLen
        BinaryMessage(1, i*j) = BinaryCharacters(i,j);
    end 
end

SampleFreq = 50000000;      %sampling frequency (2.5x max expected freq)
BitLength = .001;           %time in seconds alloted for each bit
Fs = SampleFreq*BitLength;  %ratio of samples to 1 second

X = linspace(0,BitLength,Fs);
Y = zeros([TotalLen,Fs]);
DecodedMessage = zeros([TotalLen,Fs]);
DecodedBinary = char(zeros([1,(TotalLen)]));

WaveSynthTime = tic();
for q = 1:TotalLen
    if BinaryMessage(q) == '0'
        for k = 1:Fs
            Y(q,k) = sin(20000000*2*pi*X(k));
        end
    else
        for k = 1:Fs
            Y(q,k) = sin(10000000*2*pi*X(k));
        end
    end
end
WaveSynthTimeElapsed = toc(WaveSynthTime);
disp('Time encoding:');
disp(WaveSynthTimeElapsed);

%begin decoding
DecodeTime = tic();
for w = 1:TotalLen
    DecodedMessage(w, 1:Fs) = fft(Y(w, 1:Fs));
    DecodedMessage(w, 1:Fs) = fftshift(DecodedMessage(w, 1:Fs));
    Magnitude = abs(DecodedMessage(w,1:Fs)).^2;
    [MaximumIndex, Index] = max(Magnitude(Fs/2+2:Fs));
    Index = Index / BitLength;
    if Index == 20000000
        DecodedBinary(w) = '0';
    else
        DecodedBinary(w) = '1';
    end
        
end
DecodeTimeElapsed = toc(DecodeTime);
disp('Time decoding:');
disp(DecodeTimeElapsed);

DecodedBinary = int8(DecodedBinary);

for ii = 1:StringLen
    for jj = 1:CharLen
         if DecodedBinary(ii*jj) == 48
            A(ii,jj) = 0;
         else
            A(ii,jj) = 1;
         end
    end
    %A(ii, 1:CharLen) = bin2dec(char(A(ii, 1:CharLen)));
end


    
%{
PROJECT NOTES

Reqs
    -minimum characters for text message is 100
    -'0' encoded as 20Mhz; '1' encoded as 10 MHz
    -length of time alloted for each bit is 1 ms

functions
    -fft()
        -fast fourior transform
    -bin2dec()
        -interprets a binary string and makes it into a decimal
    -dec2bin()
        -converts decimal number into binary
            -makes an array of '1' or '0' characters 
            -entering a 1xn array of numbers results in a nxm matrix of
                binary strings
                -where m is the minimum binary string length to represent
                    all of the binary numbers
    -xcorr()
        -cross correlation
        -is the autocorrelation sequence for the input vector
            -what?
                -this needs further research
    -goertzel()
        -discrete fourier transform of the input data using a second order
            goertzel algorithm
        -if data is a matrix, each column separately
    -int8()
        -convert thing into an 8bit number
    -char()
        -convert thing into a character array
    -input()
        -prompt user for input
    -disp()
        -display text or array
%}

%{
FUNCTION TESTING
EVERYTHING BELOW SHOULD BE COMMENTED OUT

%works
disp('TEST PRINT GO');      

%works
SingleNumber = 15;
SingleNumberConverted2binary = dec2bin(SingleNumber);
disp(SingleNumberConverted2binary);
disp(SingleNumberConverted2binary(2));

%works
ArrayOfNumbers = zeros([1,17]);
for n = 1:17
    ArrayOfNumbers(n) = n - 1;
end
ArrayOfNumbersConverted2binary = dec2bin(ArrayOfNumbers);
disp(ArrayOfNumbersConverted2binary);
disp(ArrayOfNumbersConverted2binary(17,1));
disp(ArrayOfNumbersConverted2binary(17,2));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%PROOF OF CONCEPT FOR FFT~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
tic();

%this is how you make a sin wave
Fs = 50000000;
BitLength = .001;
ActualFs = Fs*BitLength;

%Time = 1/Fs;
X = linspace(0,BitLength,ActualFs);
Y = zeros([1,ActualFs]);
for i = 1:ActualFs
   Y(i) = sin(20000000*2*pi*X(i));
end
disp('WE HAVE MADE THE WAVE');

Len = length(Y);

%this is how you make a sin wave
DiscFFT = fft(Y);
disp('WE HAVE FINISHED THE FOURIER TRANSFORM');

DiscFFT = fftshift(DiscFFT);
disp('WE HAVE SHIFTED THE THING');

Magnitude = abs(DiscFFT(1:Len)).^2;
disp('WE HAVE FOUNDED THE MAGNITUDE');

[MaximumIndex, Index] = max(Magnitude(ActualFs/2+2:ActualFs));
Index = Index / BitLength;
disp(Index);

f = linspace(-ActualFs/2, ActualFs/2, ActualFs);
%plot(X,DiscFFT);
plot(f,Magnitude);

toc();
disp('Finished...')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%}


