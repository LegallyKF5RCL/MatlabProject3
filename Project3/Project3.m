%{
EE-2347 Project 3
Contributers: Brian McRee, Jacob Sanchez
10/13/2013

Published open-source on Git: LegallyKF5RCL/MatlabProject3
%}

clc;
clear all;
close all;

TotalTime = tic();

FileP = fopen('Message.txt');
Message = fscanf(FileP, '%s');
disp(Message);

BinaryCharacters = dec2bin(Message);    %convert to binary
StringLen = length(BinaryCharacters);   %find how many characters we have
CharLen = 7;                            %length of char will always be 7
TotalLen = StringLen * CharLen;

%make a single array for the message
BinaryMessage = int8(zeros([StringLen,CharLen]));

for i = 1:StringLen
    for j = 1:CharLen
            BinaryMessage(i, j) = int8(BinaryCharacters(i,j));
    end
end

%define variables and constants
SampleFreq = 50000000;              %sampling frequency (2.5x max expected freq)
BitLength = .001;                   %time in seconds alloted for each bit
Fs = SampleFreq*BitLength;          %ratio of samples to 1 second

%make matricies...
AxisX = linspace(-Fs/2, Fs/2, Fs);  %variable for plotting the frequency domain
X = linspace(0,BitLength,Fs);       %create time domain
Y = zeros([TotalLen,Fs]);                           %create a zero matrix for the signals
DecodedMessage = zeros([TotalLen,Fs]);              %create a zero matrix for decoding the message
DecodedBinary = int8(zeros([StringLen,CharLen]));   %create a zero matrix for the bits after the message is decoded

WaveSynthTime = tic();              %start encoding clock


for q = 1:StringLen                     %for each row
    for h = 1:CharLen                   %for each value in the row
        if BinaryMessage(q,h) == 48     %if it is equal to 48
           for k = 1:Fs                 
               %write the value as a 20MHz sine wave
               Y(((q - 1) * CharLen + h),k) = sin(20000000*2*pi*X(k));
           end
        else
           for k = 1:Fs
               %write the value as a 10MHz sine wave
               Y(((q - 1) * CharLen + h),k) = sin(10000000*2*pi*X(k));
           end
        end
    end
end

%output time
WaveSynthTimeElapsed = toc(WaveSynthTime);      %stop encoding clock
disp('Time encoding:');                         %print time taken to encode
disp(WaveSynthTimeElapsed);                     %print time taken to encode

%begin decoding
DecodeTime = tic();         %start decode clock
for w = 1:StringLen         
    for d = 1:CharLen
        DecodedMessage(((w - 1) * CharLen + d), 1:Fs) = fft(Y(((w - 1) * CharLen + d), 1:Fs));                      %get the frequency domain
        DecodedMessage(((w - 1) * CharLen + d), 1:Fs) = fftshift(DecodedMessage(((w - 1) * CharLen + d), 1:Fs));    %shift it to where 0Hz is at center
        Magnitude = abs(DecodedMessage(((w - 1) * CharLen + d),1:Fs)).^2;                     %find the magnitudes of the freq components
        [MaximumIndex, Index] = max(Magnitude(Fs/2+2:Fs));              %find the index value where there is a maximum magnitude
        %NOTE: this accounts for positive frequencyies and offsets by two
        %for the zero frequency and the "sum" which apparently is part of fft()    
        if Index == 20000                   %if the freq is 20k (with respect to 1ms)
            DecodedBinary(w, d) = 48;       %its a 0
        else                                %otherwise
            DecodedBinary(w, d) = 49;       %its a 1
        end
    end
end

DecodeTimeElapsed = toc(DecodeTime);        %stop decode clock
disp('Time decoding:');                     %le print
disp(DecodeTimeElapsed);                    %le print the time

for ii = 1:StringLen                        %for each row
    for jj = 1:CharLen                      %for each bit in the row
         if DecodedBinary(ii,jj) == 48      %if that bit is a decimal 48
            A(ii,jj) = 0;                   %set the bit as a 0
         else                               %otherwise
            A(ii,jj) = 1;                   %set the bit as a 1
         end
    end
end

for kk = 1:StringLen                        %loop through everything
    for asdf = 1:CharLen
        B(kk,asdf) = dec2bin(A(kk,asdf));   %convert the decimals to a binary sequence
    end
end
B = bin2dec(B);     %convert the binary sequence into a decimal number that represents the message character

disp('~~~~~');
disp(Message);
disp(char(B'));
disp('~~~~~');


EndTime = toc(TotalTime);
disp('Total Time:');                     %le print
disp(EndTime);                    %le print the time
  

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


