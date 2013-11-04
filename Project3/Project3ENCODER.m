%{
EE-2347 Project 3 ENCODER
Contributers: Brian McRee, Jacob Sanchez
10/13/2013

Published open-source on Git: LegallyKF5RCL/MatlabProject3
%}

clc;
clear all;
close all;

Message = input('Input Message (must not be all numbers): ', 's');    %prompt user for input   

%initialize
BinaryCharacters = dec2bin(Message);    %convert to binary
StringLen = length(BinaryCharacters);   %find how many characters we have
CharLen = 7;                            %length of char will always be 7
TotalLen = StringLen * CharLen;         %total amount of 

%make a single array for the message
BinaryMessage = int8(zeros([StringLen,CharLen]));   

%rearrange things
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

save('EncodedMessage.txt', 'Y', '-ascii');      

disp('Finished...');