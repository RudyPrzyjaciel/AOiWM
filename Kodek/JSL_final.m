clear;
close all;

Img = imread("Tiger.jpg");
% Img = imread("Wiosna-winniczka.jpg");
[M, N, L] = size(Img);

Q = Qtable();
% IMG -> JSL

% RGB -> YCbCr
I_ycbcr = rgb2ycbcr(Img);

% Forward Anscombe
I_anscombe = Anscombe(I_ycbcr);

% Chrominance Donwsampling
I_chrom = ChromDown(I_anscombe);

% Discrete Cosine Transform & Quantisation
I_dct = zeros(size(I_chrom));
for l = 1:3
    I_dct(:,:,l) = blockproc(I_chrom(:,:,l), [8 8], @(Img) dct2(Img.data));
    I_dct(:,:,l) = blockproc(I_dct(:,:,l), [8 8], @(Img) round(Img.data ./ Q));
end

HuffmanCell = cell(0,2);
for offsetM = 1:M/8
    for offsetN = 1:N/8
        block(1:8,1:8,1:3) = I_dct((offsetM-1)*8+(1:8),(offsetN-1)*8+(1:8),1:3);
        
        % ZigZag
        I_zig = zeros(1,64,3);
        for l = 1:L
            I_zig(:,:,l) = zigzag(block(:,:,l));
        end
        
        % Huffman Encoding
        for l = 1:L
            [I_huff, dict] = HuffmanEnc(I_zig(:,:,l));
            HuffmanCell{end+1, 1} = I_huff;
            HuffmanCell{end, 2} = dict;
        end
    end
end

save("tiger.jsl","HuffmanCell","-mat");
% save("winniczek.jsl","HuffmanCell","-mat");

%% JSL -> IMG

load("tiger.jsl","HuffmanCell","-mat");
% load("winniczek.jsl","HuffmanCell","-mat");
index = 0;
for cell = 1:3:size(HuffmanCell,1)
    % Huffman Decoding
    for l = 1:L
        code = HuffmanCell{cell + (l-1),1};
        dict = HuffmanCell{cell + (l-1),2};
       DecodedHuff(:,:,l) = huffmandeco(double(code), dict);
    end
    % Inverse Zigzag
    I_Izig = zeros(8,8,3);
        for l = 1:L
            I_Izig(:,:,l) = izigzag(DecodedHuff(:,:,l),8,8);   
        end
        % Pasting block onto final image net
        %512x1024
        I_re_dct(1+8*(floor(index/128)):8 + 8*(floor(index/128)), ...
                1 + 8*(mod(index,128)):8 + 8*(mod(index,128)),:) ...
                = I_Izig(:,:,:);
        %1024x2048
%         I_re_dct(1+8*(floor(index/256)):8 + 8*(floor(index/256)), ...
%                 1 + 8*(mod(index,256)):8 + 8*(mod(index,256)),:) ...
%                 = I_Izig(:,:,:);
        index = index + 1;
end

% Inverse Quantisation & Discrete Cosine Transform 
I_idct = zeros(size(I_re_dct));
for l = 1:3
    I_idct(:,:,l) = blockproc(I_re_dct(:,:,l), [8 8], @(Img) Img.data .* Q);
    I_idct(:,:,l) = blockproc(I_idct(:,:,l), [8 8], @(Img) idct2(Img.data));
end

% Inverse Anscombe Transform
I_ians = iAnscombe(I_idct);
% Inverse Colour Conversion
I_rgb = ycbcr2rgb(I_ians);

figure(2);
montage({(Img - I_rgb)});
figure(1);
montage({Img, I_rgb});


%% fun defs
% Huffman Encoding without RLE
function [output, dictionary] = HuffmanEnc(Array) 
    symbols = unique(Array);
    sym = hist(Array,symbols); %#ok<HIST> 
    if size(symbols,2) == 1 
        dictionary = {symbols,symbols};
        output = logical(Array);
    else
        
        counts = [symbols;sym];
        counts(2,:) = normalize(counts(2,:),'norm',1);
        [dictionary, ~] = huffmandict(counts(1,:),counts(2,:));
        output = logical(huffmanenco(Array,dictionary));
    end
    if size(dictionary,2)==0
        if symbols == 0 
        dictionary = [0, 0];
        output = zeros(1,64);
        output = logical(output);
        else
        dictionary = [1, 1];
        output = ones(1,64);
        output = logical(output);
        end
    end
end

% ZigZag
function output = zigzag(in)

    % Alexey S. Sokolov a.k.a. nICKEL, Moscow, Russia
    % June 2007
    % alex.nickel@gmail.com
    % https://www.mathworks.com/matlabcentral/fileexchange/15317-zigzag-scan
    % initializing the variables
    %----------------------------------
    h = 1;
    v = 1;
    vmin = 1;
    hmin = 1;
    vmax = size(in, 1);
    hmax = size(in, 2);
    i = 1;
    output = zeros(1, vmax * hmax);
    %----------------------------------
    while ((v <= vmax) && (h <= hmax))
        
        if (mod(h + v, 2) == 0)                 % going up
            if (v == vmin)       
                output(i) = in(v, h);        % if we got to the first line
                if (h == hmax)
	          v = v + 1;
	        else
                  h = h + 1;
                end
                i = i + 1;
            elseif ((h == hmax) && (v < vmax))   % if we got to the last column
                output(i) = in(v, h);
                v = v + 1;
                i = i + 1;
            elseif ((v > vmin) && (h < hmax))    % all other cases
                output(i) = in(v, h);
                v = v - 1;
                h = h + 1;
                i = i + 1;
            end
            
        else                                    % going down
           if ((v == vmax) && (h <= hmax))       % if we got to the last line
                output(i) = in(v, h);
                h = h + 1;
                i = i + 1;
            
           elseif (h == hmin)                   % if we got to the first column
                output(i) = in(v, h);
                if (v == vmax)
	          h = h + 1;
	        else
                  v = v + 1;
                end
                i = i + 1;
           elseif ((v < vmax) && (h > hmin))     % all other cases
                output(i) = in(v, h);
                v = v + 1;
                h = h - 1;
                i = i + 1;
           end
        end
        if ((v == vmax) && (h == hmax))          % bottom right element
            output(i) = in(v, h);
            break
        end
    end
end

% Inverse ZigZag
function output = izigzag(in, vmax, hmax)
% Alexey S. Sokolov a.k.a. nICKEL, Moscow, Russia
% June 2007
% alex.nickel@gmail.com
% https://www.mathworks.com/matlabcentral/fileexchange/15317-zigzag-scan
    % initializing the variables
    %----------------------------------
    h = 1;
    v = 1;
    vmin = 1;
    hmin = 1;
    output = zeros(vmax, hmax);
    i = 1;
    %----------------------------------
    while ((v <= vmax) && (h <= hmax))
        if (mod(h + v, 2) == 0)                % going up
            if (v == vmin)
                output(v, h) = in(i);
                if (h == hmax)
	          v = v + 1;
	        else
                  h = h + 1;
                end
                i = i + 1;
            elseif ((h == hmax) && (v < vmax))
                output(v, h) = in(i);
                i;
                v = v + 1;
                i = i + 1;
            elseif ((v > vmin) && (h < hmax))
                output(v, h) = in(i);
                v = v - 1;
                h = h + 1;
                i = i + 1;
            end
            
        else                                   % going down
           if ((v == vmax) && (h <= hmax))
                output(v, h) = in(i);
                h = h + 1;
                i = i + 1;
            
           elseif (h == hmin)
                output(v, h) = in(i);
                if (v == vmax)
	          h = h + 1;
	        else
                  v = v + 1;
                end
                i = i + 1;
           elseif ((v < vmax) && (h > hmin))
                output(v, h) = in(i);
                v = v + 1;
                h = h - 1;
                i = i + 1;
           end
        end
        if ((v == vmax) && (h == hmax))
            output(v, h) = in(i);
            break
        end
    end
end

% Chrominance Downsampling
function output = ChromDown(Img)
    [M, N, ~] = size(Img);
    output = zeros(size(Img));
    newChrom = zeros(M/2, N/2, 2);
    box = zeros(2,2);
    for i = 1:2 % 1: Cb, 2: Cr
        for m = 1:M/2
            for n = 1:N/2
                for mb = 1:2
                    for nb = 1:2
                        box(mb,nb) = Img(2*(m-1)+mb, 2*(n-1)+nb,i+1);
                    end
                end
                newChrom(m,n,i) = mean(box, "all");
            end
        end
    end
    newChrom = imresize(newChrom,2);

    output(:,:,1) = Img(:,:,1);
    output(:,:,2) = newChrom(:,:,1);
    output(:,:,3) = newChrom(:,:,2);
end

% Anscombe transform
function Ansc = Anscombe(Img)
    Img = double(Img);
    Ansc = 2 * real((Img + 3/8).^(1/2));
end

% inverse Anscombe transform
function iAnsc = iAnscombe(Img)
    Img = 1/4 * Img.^2 - 1/8 + 1/4 * real((3/2)^(1/2)) * 1/Img - 11/8 * ...
        1/Img.^2 + 5/8 * real((3/2)^(1/2))*1/Img.^3;
    iAnsc= uint8(Img);
end

% Quantisation Table
function Q = Qtable()
%     Q = [1 1 1 1 1 1 1 1;
%         1 1 1 1 1 1 1 1;
%         1 1 1 1 1 1 1 1;
%         1 1 1 1 1 1 1 1;
%         1 1 1 1 1 1 1 1;
%         1 1 1 1 1 1 1 1;
%         1 1 1 1 1 1 1 1;
%         1 1 1 1 1 1 1 1];
    Q = [4, 3, 4, 4, 4, 6, 11, 15;
        3, 3, 3, 4, 5, 8, 14, 19;
        3, 4, 4, 5, 8, 12, 16, 20;
        4, 5, 6, 7, 12, 14, 18, 20;
        6, 6, 9, 11, 14, 17, 21, 23;
        9, 12, 12, 18, 23, 22, 25, 21;
        11, 13, 15, 17, 21, 23, 25, 21;
        13, 12, 12, 13, 16, 19, 21, 21];
%     Q = [16 11 10 16 24 40 51 61;
%          12 12 14 19 26 58 60 55;
%          14 13 16 24 40 57 69 56;
%          14 17 22 29 51 87 80 62;
%          18 22 37 56 68 109 103 77;
%          24 35 55 64 81 104 113 92;
%          49 64 78 87 103 121 120 101;
%          72 92 95 98 112 100 103 99];
end