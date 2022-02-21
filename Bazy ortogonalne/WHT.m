clear;
close;

%% Prep
Ref = imread("green_frog.jpg");
Pois = imread("green_frog_poiss.jpg");
[M, N, L] = size(Pois);

Err(1,1) = "WHT";
Err(1,2) = "immse to Ref";
Err(1,3) = "immse to Pois";
Err(end+1,1) = "Pois";
Err(end,2) = immse(Pois,Ref);
Err(end,3) = immse(Pois,Pois);
%% forward Anscombe transform
Pois = double(Pois);
Anscombe = 2 * real((Pois + 3/8).^(1/2));
Anscombe = uint8(Anscombe);
Pois = uint8(Pois);


%% forward WHT
T1 = 0.1;
T2 = 1;
Q1 = 1;
Q2 = 64;

%wht requires square matrix size of n = 2^x
x = max(M,N);
x=nextpow2(x);
x=2^x;
WHT_Ref = zeros(x,x);
Anscombe_wht = zeros(x,x);
[R, C, K] = size(Anscombe);

for r = 1:R
    for c = 1:C
        for k = 1:K
            Anscombe_wht(r,c,k) = Anscombe(r,c,k);
        end
    end
end

WHT_Ref(:,:,1) = fwht(double(Anscombe_wht(:,:,1)));
WHT_Ref(:,:,2) = fwht(double(Anscombe_wht(:,:,2)));
WHT_Ref(:,:,3) = fwht(double(Anscombe_wht(:,:,3)));


WHT_T1 = WHT_Ref;
WHT_T1(abs(WHT_T1)<=T1) = 0;
figure(2),subplot(3,2,1), imshow(log(abs(WHT_T1)),'DisplayRange', [0 255]);
% figure(2),subplot(3,2,3), imshow(WHT_T1,'DisplayRange', [0 255]);
title("WHT\_T1");

WHT_Q1 = WHT_Ref;
WHT_Q1 = (Q1.*WHT_Q1 + 1/2)./Q1;

figure(2),subplot(3,2,2), imshow(log(abs(WHT_Q1)),'DisplayRange', [0 255]);
% figure(2),subplot(3,2,4), imshow(WHT_Q1,'DisplayRange', [0 255]);
title("WHT\_Q1");

WHT_T2 = WHT_Ref;
WHT_T2(abs(WHT_T2)<=T2) = 0;
figure(2),subplot(3,2,3), imshow(log(abs(WHT_T2)),'DisplayRange', [0 255]);
% figure(2),subplot(3,2,5), imshow(WHT_T2,'DisplayRange', [0 255]);
title("WHT\_T2");

WHT_Q2 = WHT_Ref;
WHT_Q2 = (Q2.*WHT_Q2 + 1/2)./Q2;

figure(2),subplot(3,2,4), imshow(log(abs(WHT_Q2)),'DisplayRange', [0 255]);
% figure(2),subplot(3,2,6), imshow(WHT_Q1,'DisplayRange', [0 255]);
title("WHT\_Q2");
%% inverse WHT and Anscombe
I_T1 = zeros(x,x,L);
for l=1:L
    I_T1(:,:,l) = ifwht(WHT_T1(:,:,l));
end
% I_T1 = 1/4 * I_T1.^2 - 1/8;
I_T1 = 1/4 * I_T1.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_T1 - 11/8 * ...
    1/I_T1.^2 + 5/8 * real((3/2)^(1/2))*1/I_T1.^3;
I_T1 = I_T1(1:R,1:C,:);
I_T1 = uint8(I_T1);
Err(end+1,1) = "I_T1";
Err(end,2) = immse(I_T1,Ref);
Err(end,3) = immse(I_T1,Pois);


I_Q1 = zeros(x,x,L);
for l=1:L
    I_Q1(:,:,l) = ifwht(WHT_Q1(:,:,l));
end
% I_Q1 = 1/4 * I_Q1.^2 - 1/8;
I_Q1 = 1/4 * I_Q1.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_Q1 - 11/8 * ...
    1/I_Q1.^2 + 5/8 * real((3/2)^(1/2))*1/I_Q1.^3;
I_Q1 = uint8(I_Q1);

I_Q1 = I_Q1(1:R,1:C,:);
Err(end+1,1) = "I_Q1";
Err(end,2) = immse(I_Q1,Ref);
Err(end,3) = immse(I_Q1,Pois);


I_T2 = zeros(x,x,L);
for l=1:L
    I_T2(:,:,l) = ifwht(WHT_T2(:,:,l));
end
% I_T2 = 1/4 * I_T2.^2 - 1/8;
I_T2 = 1/4 * I_T2.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_T2 - 11/8 * ...
    1/I_T2.^2 + 5/8 * real((3/2)^(1/2))*1/I_T2.^3;
I_T2 = uint8(I_T2);

I_T2 = I_T2(1:R,1:C,:);
Err(end+1,1) = "I_T2";
Err(end,2) = immse(I_T2,Ref);
Err(end,3) = immse(I_T2,Pois);


I_Q2 = zeros(x,x,L);
for l=1:L
    I_Q2(:,:,l) = ifwht(WHT_Q2(:,:,l));
end
% I_Q2 = 1/4 * I_Q2.^2 - 1/8;
I_Q2 = 1/4 * I_Q2.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_Q2 - 11/8 * ...
    1/I_Q2.^2 + 5/8 * real((3/2)^(1/2))*1/I_Q2.^3;
I_Q2 = uint8(I_Q2);

I_Q2 = I_Q2(1:R,1:C,:);
Err(end+1,1) = "I_Q2";
Err(end,2) = immse(I_Q2,Ref);
Err(end,3) = immse(I_Q2,Pois);
%% show images

% differences
figure(3), subplot(3,2,1), imshow(uint8(Ref-Pois)*10);
title("Ref-Pois");
figure(3), subplot(3,2,3), imshow(uint8(Ref-I_T1)*10);
title("Diff\_T1");
figure(3), subplot(3,2,4), imshow(uint8(Ref-I_Q1)*10);
title("Diff\_Q1");
figure(3), subplot(3,2,5), imshow(uint8(Ref-I_T2)*10);
title("Diff\_T2");
figure(3), subplot(3,2,6), imshow(uint8(Ref-I_Q2)*10);
title("Diff\_Q2");
% display images
figure(1)
subplot(3,2,1), imshow(Ref),title("Ref");
subplot(3,2,2), imshow(Pois),title("Pois");
subplot(3,2,3), imshow(I_T1),title("I\_T1");
figure(1), subplot(3,2,4), imshow(I_Q1);
title("I\_Q1");
figure(1), subplot(3,2,5), imshow(I_T2);
title("I\_T2");
figure(1), subplot(3,2,6), imshow(I_Q2);
title("I\_Q2");


% %differences
% figure(3), subplot(2,2,1), imshow(uint8(Ref-I_T1)*10);
% title("Diff\_T1");
% figure(3), subplot(2,2,2), imshow(uint8(Ref-I_Q1)*10);
% title("Diff\_Q1");
% figure(3), subplot(2,2,3), imshow(uint8(Ref-I_T2)*10);
% title("Diff\_T2");
% figure(3), subplot(2,2,4), imshow(uint8(Ref-I_Q2)*10);
% % images
% figure(1), subplot(2,2,1), imshow(I_T1);
% title("I\_T1");
% 
% subplot(2,2,2), imshow(I_Q1);
% title("I\_Q1");
% 
% figure(1), subplot(2,2,3), imshow(I_T2);
% title("I\_T2");
% 
% subplot(2,2,4), imshow(I_Q2);
% title("I\_Q2");
