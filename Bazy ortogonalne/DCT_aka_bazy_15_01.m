clear;
close;

%% Prep
Ref = imread("green_frog.jpg");
Pois = imread("green_frog_poiss.jpg");
[M, N, L] = size(Pois);

Err(1,1) = "DCT";
Err(1,2) = "immse to Ref";
Err(1,3) = "immse to Pois";
Err(end+1,1) = "Pois";
Err(end,2) = immse(Pois,Ref);
Err(end,3) = immse(Pois,Pois);

%% Anscombe transform
Pois = double(Pois);
Anscombe = 2 * real((Pois + 3/8).^(1/2));
% Anscombe = uint8(Anscombe);
Pois = uint8(Pois);

%% DCT
T1 = 1;
T2 = 8;
Q1 = 1;
Q2 = 8;
DCTRef = zeros(M,N,L);
% transform
for l = 1:L
    DCTRef(:,:,l) = dct2(Anscombe(:,:,l));
end
figure(2),subplot(3,2,1), imshow(log(abs(DCTRef)),'DisplayRange', [0 255]);
title("Ref");
% operations on coefficients
DCT_T1 = DCTRef;
DCT_T1(abs(DCT_T1)<=T1) = 0;
figure(2),subplot(3,2,3), imshow(log(abs(DCT_T1)),'DisplayRange', [0 255]);
title("T1");

DCT_Q1 = DCTRef;
DCT_Q1 = (Q1.*DCT_Q1 + 1/2)./Q1;
figure(2),subplot(3,2,4), imshow(log(abs(DCT_Q1)),'DisplayRange', [0 255]);
title("Q1");

DCT_T2 = DCTRef;
DCT_T2(abs(DCT_T2)<=T2) = 0;
figure(2),subplot(3,2,5), imshow(log(abs(DCT_T2)),'DisplayRange', [0 255]);
title("T2");

DCT_Q2 = DCTRef;
DCT_Q2 = (Q2.*DCT_Q2 + 1/2)./Q2;
figure(2),subplot(3,2,6), imshow(log(abs(DCT_Q2)),'DisplayRange', [0 255]);
title("Q2");

% inverse transforms
I_T1 = zeros(M,N,L);
for l=1:L
    I_T1(:,:,l) = idct2(DCT_T1(:,:,l));
end
% I_T1 = 1/4 * I_T1.^2 - 1/8;
I_T1 = 1/4 * I_T1.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_T1 - 11/8 * ...
    1/I_T1.^2 + 5/8 * real((3/2)^(1/2))*1/I_T1.^3;
I_T1 = uint8(I_T1);
Err(end+1,1) = "I_T1";
Err(end,2) = immse(I_T1,Ref);
Err(end,3) = immse(I_T1,Pois);

I_Q1 = zeros(M,N,L);
for l=1:L
    I_Q1(:,:,l) = idct2(DCT_Q1(:,:,l));
end
% I_Q1 = 1/4 * I_Q1.^2 - 1/8;
I_Q1 = 1/4 * I_Q1.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_Q1 - 11/8 * ...
    1/I_Q1.^2 + 5/8 * real((3/2)^(1/2))*1/I_Q1.^3;
I_Q1 = uint8(I_Q1);
Err(end+1,1) = "I_Q1";
Err(end,2) = immse(I_Q1,Ref);
Err(end,3) = immse(I_Q1,Pois);

I_T2 = zeros(M,N,L);
for l=1:L
    I_T2(:,:,l) = idct2(DCT_T2(:,:,l));
end
% I_T2 = 1/4 * I_T2.^2 - 1/8;
I_T2 = 1/4 * I_T2.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_T2 - 11/8 * ...
    1/I_T2.^2 + 5/8 * real((3/2)^(1/2))*1/I_T2.^3;
I_T2 = uint8(I_T2);
Err(end+1,1) = "I_T2";
Err(end,2) = immse(I_T2,Ref);
Err(end,3) = immse(I_T2,Pois);

I_Q2 = zeros(M,N,L);
for l=1:L
    I_Q2(:,:,l) = idct2(DCT_Q2(:,:,l));
end
% I_Q2 = 1/4 * I_Q2.^2 - 1/8;
I_Q2 = 1/4 * I_Q2.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/I_Q2 - 11/8 * ...
    1/I_Q2.^2 + 5/8 * real((3/2)^(1/2))*1/I_Q2.^3;
I_Q2 = uint8(I_Q2);
Err(end+1,1) = "I_Q2";
Err(end,2) = immse(I_Q2,Ref);
Err(end,3) = immse(I_Q2,Pois);

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

