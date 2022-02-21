clear;
close;

%% Prep
Ref = imread("green_frog.jpg");
Pois = imread("green_frog_poiss.jpg");

% wname = "bior1.1";
% wname = "bior2.2";
wname = "bior4.4";
Err(1,1) = wname;
Err(1,2) = "immse to Ref";
Err(1,3) = "immse to Pois";
Err(end,3) = "Non zeros";
Err(end+1,1) = "Pois";
Err(end,2) = immse(Pois,Ref);
Err(end,3) = immse(Pois,Pois);

[I_T1, T1_non]= falkoweT(Pois, 3, wname);
Err(end+1,1) = "I_T1";
Err(end,2) = immse(I_T1,Ref);
Err(end,3) = immse(I_T1,Pois);
Err(end,4) = T1_non;

[I_T2, T2_non] = falkoweT(Pois, 50, wname);
Err(end+1,1) = "I_T2";
Err(end,2) = immse(I_T2,Ref);
Err(end,3) = immse(I_T2,Pois);
Err(end,4) = T2_non;

[I_Q1, Q1_non] = falkoweQ(Pois, 2^4, wname);
Err(end+1,1) = "I_Q1";
Err(end,2) = immse(I_Q1,Ref);
Err(end,3) = immse(I_Q1,Pois);
Err(end,4) = Q1_non;

[I_Q2, Q2_non] = falkoweQ(Pois, 2^5, wname);
Err(end+1,1) = "I_Q2";
Err(end,2) = immse(I_Q2,Ref);
Err(end,3) = immse(I_Q2,Pois);
Err(end,4) = Q2_non;


figure(1);
subplot(3,2,1), imshow(Ref),title("Ref");
subplot(3,2,2), imshow(Pois),title("Pois");
subplot(3,2,3), imshow(I_T1),title("I\_T1");
subplot(3,2,4), imshow(I_Q1),title("I\_Q1");
subplot(3,2,5), imshow(I_T2),title("I\_T2");
subplot(3,2,6), imshow(I_Q2),title("I\_Q2");


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

function [I_T, non] = falkoweT(Img, T, wname)
    Img_ansc = Anscombe(Img);
    [cA, cH, cV, cD] = dwt2(Img_ansc, wname);
    cA_T = cA;
    cA_T(abs(cA_T)<=T) = 0;
    cH_T = cH;
    cH_T(abs(cH_T)<=T) = 0;
    cV_T = cV;
    cV_T(abs(cV_T)<=T) = 0;
    cD_T = cD;
    cD_T(abs(cD_T)<=T) = 0;
    
    Anon = nnz(cA_T);
    Hnon = nnz(cH_T);
    Vnon = nnz(cV_T);
    Dnon = nnz(cD_T);
    non = Anon+Hnon+Vnon+Dnon;

    I_T = idwt2(cA_T, cH_T, cV_T, cD_T, wname);
    I_T = iAnscombe(I_T);
end

function [I_Q, non] = falkoweQ(Img, Q, wname)
    Img_ansc = Anscombe(Img);
    [cA, cH, cV, cD] = dwt2(Img_ansc, wname);
    cA_Q = (Q .* cA + 1/2)./Q;
    cH_Q = (Q .* cH + 1/2)./Q;
    cV_Q = (Q .* cV + 1/2)./Q;
    cD_Q = (Q .* cD + 1/2)./Q;
    
    Anon = nnz(cA_Q);
    Hnon = nnz(cH_Q);
    Vnon = nnz(cV_Q);
    Dnon = nnz(cD_Q);
    non = Anon+Hnon+Vnon+Dnon;

    I_Q = idwt2(cA_Q, cH_Q, cV_Q, cD_Q, wname);
    I_Q = iAnscombe(I_Q);
end

function Ansc = Anscombe(Img)
    Img = double(Img);
    Ansc = 2 * real((Img + 3/8).^(1/2));
end

function iAnsc = iAnscombe(Img)
    Img = 1/4 * Img.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/Img - 11/8 * ...
        1/Img.^2 + 5/8 * real((3/2)^(1/2))*1/Img.^3;
    iAnsc= uint8(Img);
end
