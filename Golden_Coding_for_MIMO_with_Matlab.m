function output=goldenml();
clear
SNRbit=0:2:20;
M=4; %PSK
Mx=M;
L=4;
SNRsyb=log2(M)*(10.^(SNRbit/10))/2; % average symbol snr per transmit antenna
datalength=10^4;
q=(1+sqrt(5))/2;
qe=(1-sqrt(5))/2;
alpha=1+sqrt(-1)-sqrt(-1)*q;
alphae=1+sqrt(-1)-sqrt(-1)*qe;
symbRange=symbRangef;
for a=1:length(SNRbit)
errorx=0;
errorbit=0;
insym1=randint(1,datalength,[0 M-1]); % symbols1
insym2=randint(1,datalength,[0 M-1]); % symbols2
insym3=randint(1,datalength,[0 M-1]); % symbols1
insym4=randint(1,datalength,[0 M-1]); % symbols2

trsym1=dmodce(insym1,1,1,'psk',M); %modulated symbols1
trsym2=dmodce(insym2,1,1,'psk',M); %modulated symbols2
trsym3=dmodce(insym3,1,1,'psk',M); %modulated symbols1
trsym4=dmodce(insym4,1,1,'psk',M); %modulated symbols2

noise =sqrt(0.5/SNRsyb(a))*(randn(L,datalength) + i*randn(L,datalength));
channel = sqrt(0.5)*(randn(L,datalength)+i*randn(L,datalength));

r1=channel(1,:).*a11 +channel(2,:).*a21+ noise(1,:);
r2=channel(1,:).*a12+channel(2,:).*a22+noise(2,:);
r3=channel(3,:).*a11 +channel(4,:).*a21+ noise(3,:);
r4=channel(3,:).*a12+channel(4,:).*a22+noise(4,:);

for x=1:256
    trsym1a=symbRange(x,1);
    trsym2a=symbRange(x,2);
    trsym3a=symbRange(x,3);
    trsym4a=symbRange(x,4);
    sum1(x,:)=abs(r1-channel(1,:).*sqrt(0.2).*alpha.*(trsym1a+trsym2a*q)-channel(2,:).*sqrt(0.2)*sqrt(-1)*alphae*(trsym3a+trsym4a*qe)).^2 + abs(r3-channel(3,:).*sqrt(0.2).*alpha.*(trsym1a+trsym2a*q)-channel(4,:).*sqrt(0.2).*sqrt(-1)*alphae.*(trsym3a+trsym4a*qe)).^2+abs(r2-channel(1,:).*sqrt(0.2).*alpha.*(trsym3a+trsym4a*q)-channel(2,:).*sqrt(0.2).*alphae.*(trsym1a+trsym2a*qe)).^2 +abs(r4-channel(3,:).*sqrt(0.2).*alpha.*(trsym3a+trsym4a*q)-channel(4,:).*sqrt(0.2).*alphae.*(trsym1a+trsym2a*qe)).^2;
end

[D1 F1]=min(sum1,[],1); %D min value, F: min value row

mindistind1a=ddemodce([symbRange(F1,1) ],1,1,'psk',Mx) ;
mindistind2a=ddemodce([symbRange(F1,2) ],1,1,'psk',Mx);
mindistind3a=ddemodce([symbRange(F1,3) ],1,1,'psk',Mx);
mindistind4a=ddemodce([symbRange(F1,4) ],1,1,'psk',Mx);

error1= sum(transpose(mindistind1a(1,1))~=insym1);
error2= sum(transpose(mindistind2a(1,1))~=insym2);
error3= sum(transpose(mindistind3a(1,1))~=insym3);
error4= sum(transpose(mindistind4a(1,1))~=insym4);

errorx=errorx+error1+error2+error3+error4;

gray_code_mindist = [];
gray_code_insym = [];

for c=1:1:datalength
    gray_code_mindist1(1,2*c-1:1:2*c)= binary_to_gray_code(convert_decimal_to_binary(mindistind1a(c)));
    gray_code_insym1(1,2*c-1:1:2*c) = binary_to_gray_code(convert_decimal_to_binary(insym1(1,c)));
    gray_code_mindist2(1,2*c-1:1:2*c)= binary_to_gray_code(convert_decimal_to_binary(mindistind2a(c)));
    gray_code_insym2(1,2*c-1:1:2*c) = binary_to_gray_code(convert_decimal_to_binary(insym2(1,c)));
    gray_code_mindist3(1,2*c-1:1:2*c)= binary_to_gray_code(convert_decimal_to_binary(mindistind3a(c)));
    gray_code_insym3(1,2*c-1:1:2*c) = binary_to_gray_code(convert_decimal_to_binary(insym3(1,c)));
    gray_code_mindist4(1,2*c-1:1:2*c)= binary_to_gray_code(convert_decimal_to_binary(mindistind4a(c)));
    gray_code_insym4(1,2*c-1:1:2*c) = binary_to_gray_code(convert_decimal_to_binary(insym4(1,c)));
end

error3a= sum(gray_code_mindist1~= gray_code_insym1);
error4a=sum(gray_code_mindist2~= gray_code_insym2);
error5a= sum(gray_code_mindist3~= gray_code_insym3);
error6a=sum(gray_code_mindist4~= gray_code_insym4);

errorbit=error3a+error5a+error4a+error6a;

SER2(a)=(errorx)/(datalength*4);
BER(a)=(errorbit)/(4*datalength*log2(M));
end

semilogy(SNRbit,BER,'k*-');
xlabel('SNR,dB');
ylabel('SER');
title('BPSK, 2X2 ');

function symbRangef=symbRangef()
b=1;
c=-1;
d=sqrt(-1);
e=-sqrt(-1);
ax=[ b c d e];
for kx=1:1:4
    gx=ax(:,kx);
    for cx=1:1:4
         bx=ax(:,cx);
          for dx=1:1:4
         ex=ax(:,dx);
            for mx=1:1:4
                nx=ax(:,mx);
   
symbRangef(64*(kx-1)+16*(cx-1)+4*(dx-1)+mx,1:1:4)=   [gx bx ex nx];
            end
          end
    end
end 
