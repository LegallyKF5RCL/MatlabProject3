%{
EE-2347 Project 3 DECODER
Contributers: Brian McRee, Jacob Sanchez
10/13/2013

Published open-source on Git: LegallyKF5RCL/MatlabProject3
%}

%clc;
clear all;
close all;

%load the data
LoadingTime = tic();            %start loading timer
Y = load('EncodedMessage.txt'); %load the file zzzz
LoadingEnd = toc(LoadingTime);  %end loading timer
disp('Loading Time: ');         %print
disp(LoadingEnd);               %display time elapsed for loading

%declare variables and constants
StringLen = size(Y)/7;                  %find how many characters we have
StringLen = StringLen(1,1);             %this index holds the actual value
CharLen = 7;                            %length of char will always be 7
TotalLen = StringLen * CharLen;         %total length in the total amount of elements
SampleFreq = 50000000;              %sampling frequency (2.5x max expected freq)
BitLength = .001;                   %time in seconds alloted for each bit
Fs = SampleFreq*BitLength;          %ratio of samples to 1 second

%make matricies...
AxisX = linspace(-Fs/2, Fs/2, Fs);              %variable for plotting the frequency domain
X = linspace(0,BitLength,Fs);                   %create time domain
%Y = zeros([TotalLen,Fs]);                          %create a zero matrix for the signals
DecodedMessage = zeros([TotalLen,Fs]);              %create a zero matrix for decoding the message
DecodedBinary = int8(zeros([StringLen,CharLen]));   %create a zero matrix for the bits after the message is decoded


%begin decoding
DecodeTime = tic();         %start decode clock
for w = 1:StringLen         
    for d = 1:CharLen
        DecodedMessage(((w - 1) * CharLen + d), 1:Fs) = fft(Y(((w - 1) * CharLen + d), 1:Fs));                      %get the frequency domain
        DecodedMessage(((w - 1) * CharLen + d), 1:Fs) = fftshift(DecodedMessage(((w - 1) * CharLen + d), 1:Fs));    %shift it to where 0Hz is at center
        Magnitude = abs(DecodedMessage(((w - 1) * CharLen + d),1:Fs)).^2;                     %find the magnitudes of the freq components
        [MaximumIndex, Index] = max(Magnitude(Fs/2+2:Fs));              %find the index value where the maximum magnitude is
        %NOTE: this accounts for positive frequencies and offsets by two
        %for the zero frequency and the "sum" which apparently is part of fft()    
        if Index == 20000                   %if the freq is 20k (with respect to 1ms)
            DecodedBinary(w, d) = 48;       %its a 0
        else                                %otherwise
            DecodedBinary(w, d) = 49;       %its a 1
        end
    end
end

DecodeTimeElapsed = toc(DecodeTime);        %stop decode clock
disp('Time decoding:');                     %print
disp(DecodeTimeElapsed);                    %print the time

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
disp(char(B'));
disp('~~~~~');

disp('Finished...');