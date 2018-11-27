datalength=total_data_length;
symbRange=dmodce([0:M-1],1,1,'psk',M);
tt=number_of_iterations; ;
SNR=28;
snr=(10^(SNR/10));   % calculate current SNRsd
E=snr;
for a=1:length(G)
 toterror1=0;
 toterror2=0;
 toterror3=0;
 toterror4=0;
 totError1=0;
g=10^(G(a)/10);
Grd=(1+sqrt(g^(-1)))^2;
Gsr=Grd*g;
Gsd=1;
AN1=2*Gsr*Kt;
AN2=Ks*(1+2*Gsr*Kt*snr)./(Grd*(1-Ks)*snr);
AD=(1+2*Gsr*Kt*snr)./(2*Grd*(1-Kt)*(1-Ks)*snr);
A0=2*Kt;

tic
for trial=1:1:tt 

insym1=randint(1,datalength,[0 M-1]); % symbols1
insym2=randint(1,datalength,[0 M-1]); % symbols2

trsym1=dmodce(insym1,1,1,'psk',M); %modulated symbols1
trsym2=dmodce(insym2,1,1,'psk',M);  %modulated symbols2


noise =sqrt(0.5)*(randn(4,datalength) + sqrt(-1)*randn(4,datalength));
channel = sqrt(0.5)*(randn(3,datalength)+sqrt(-1)*randn(3,datalength));
channelsr = sqrt(0.5)*(randn(1,1)+sqrt(-1)*randn(1,1));%
channelrd=sqrt(0.5)*(randn(1,1)+sqrt(-1)*randn(1,1));
channelsd=  sqrt(0.5)*(randn(1,1)+sqrt(-1)*randn(1,1));%SD
channel(1,:)=channelsr;
channel(2,:)=channelrd;
channel(3,:)=channelsd;

A1= AN1./(AD+abs(channelrd).^2);
A2= AN2./(AD+abs(channelrd).^2);

 h0=sqrt(E)*sqrt(A1).*channel(1,:).*channel(2,:);
 h1=sqrt(E)*sqrt(A2).*channel(3,:);
 h3=sqrt(E)*sqrt(A0)*channel(3,:);


r1=h3.*trsym1+ noise(1,:);
r2=h0.*trsym1 +h1.*trsym2+ noise(2,:);
r3=h3.*trsym2+ noise(3,:);
r4=h0.*conj(-trsym2)+h1.*conj(trsym1)+noise(4,:);


d1=conj(h0).*r2+ h1.*conj(r4); 
d2=conj(h1).*r2-h0.*conj(r4);

for x=1:M
sum1(x,:)=abs(d1-symbRange(x)).^2 + abs(r1-h3*symbRange(x)).^2;
sum2(x,:)=abs(d2-symbRange(x)).^2+ abs(r3-h3*symbRange(x)).^2;   
end

[D1 F1]=min(sum1,[],1);  %D min value, F: min value ya sahip olan row
[D2 F2]=min(sum2,[],1);

mindistind1=F1-1;
mindistind2=F2-1;

error1= sum(mindistind1(1:1:datalength)~=insym1(1:1:datalength));
error2= sum(mindistind2(1:1:datalength)~=insym2(1:1:datalength));

toterror1=toterror1+error1;
toterror2=toterror2+error2;

gray_code_mindist = [];
gray_code_insym = [];

for c=1:1:datalength
gray_code_mindist1(1,K*c-(K-1):1:K*c)= binary_to_gray_code(convert_decimal_to_binary(mindistind1(c),K));
gray_code_insym1(1,K*c-(K-1):1:K*c) = binary_to_gray_code(convert_decimal_to_binary(insym1(c),K));
gray_code_mindist2(1,K*c-(K-1):1:K*c)= binary_to_gray_code(convert_decimal_to_binary(mindistind2(c),K));
gray_code_insym2(1,K*c-(K-1):1:K*c) = binary_to_gray_code(convert_decimal_to_binary(insym2(c),K));
end

errorbit= sum(gray_code_mindist1~= gray_code_insym1)+sum(gray_code_mindist2~= gray_code_insym2);
totError1= totError1 +errorbit;
end
SER2(a)=(toterror1+toterror2)/((datalength)*2*tt);
BER(a)=totError1/(K*(datalength)*2*tt);
toc
end
