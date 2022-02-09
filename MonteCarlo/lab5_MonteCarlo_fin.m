clear;
% close;

I = uint8(imread("tesla_front.jpg"));
I = imcomplement(I);
[M, N, L] = size(I);

treshold = 67;


GS = GrayScale(I);
Poiss = Poisson(GS, -7);

Filtered = Median(Poiss);
Filtered = Convolution(Filtered);

for i = 1:3
    Filtered = Median(Filtered);
    Filtered = Convolution(Filtered);
end
FilteredMC = BWTresh(Filtered, treshold);

frontSurface = MonteCarlo(FilteredMC, 10000, 100);

tic
    points = nnz(FilteredMC);
    frontSurface = 100 * (points/(M*N));
t = toc;

figure(1);
subplot(2,3,1), imshow(Poiss);
subplot(2,3,2), imshow(Filtered);
subplot(2,3,3), imshow(FilteredMC);

function GS = GrayScale(Img)
[M, N, L] = size(Img);
GS = zeros(M,N);

    for m = 1:M
        for n = 1:N
            GS(m,n,:) = sum(Img(m,n,:)/L);
        end
    end
GS = uint8(GS);
end

function Ps = Poisson(Img, i)
    [M, N, L] = size(Img);
    P = zeros(M,N); 
    Lambda = 2.^(i);
    
    for m = 1:M
        for n = 1:N
            P(m,n) = poissrnd(Lambda * (sum(Img(m,n,:))/L))/(Lambda);
        end
    end
    Ps = uint8(P);
end

function Convo = Convolution(Img)
troj = @(x) (abs(x) < 1) .* (1 - abs(x));

%setting triangle filter for 1D
n = 7; 
xs = linspace(-1,1,n);
X = troj(xs);

%setting above filter into 2D
Z = X.' * X;

%normalizing filter
Z(:) = Z./(sum(Z))/(n-2);
Z(isnan(Z))=0;

Convo = imfilter(Img, Z);
end

function Med = Median(Img)

[M, N, L] = size(Img);
%size of filter window
s = 5;
window = zeros(s,s);

edgex = floor(s/2);
edgey = floor(s/2);

Med = zeros(M,N);

for x = s:(M-s)
   for y = s:(N-s)
       i = 1;
       j = 1;
       for fx = 1:s
           for fy = 1:s
               window(i, j) = Img(x + fx - edgex,y + fy - edgey);
               i = i + 1;
           end
           j = j + 1;
           i = 1;
       end
       med_vec = reshape(window', s.*s, 1);
       med_vec = sort(med_vec);
       Med(x,y) = med_vec((s.*s + 1)/2);
   end
end
Med = uint8(Med);
end

function blackWhite = BWTresh(Img, treshold)
[M,N,L] = size(Img);
blackWhite = zeros(M,N);
for m = 1:M
    for n = 1:N
        if (sum(Img(m,n,:)/L) < treshold)
            blackWhite(m,n) = 0;
        else
            blackWhite(m,n) = 255;
        end
    end
end
blackWhite = uint8(blackWhite);
end

function monte = MonteCarlo(Img, nPoints, SurfaceOfImg)
    tic

    [M,N,L] = size(Img);
    sz = [N,M];
    [Ir,Ic] = ind2sub(sz,randperm(prod(sz),nPoints));
    I = [Ir',Ic']';

    monteCount = 0;
    for n = 1:nPoints
        if Img(I(2,n), I(1,n)) ~= 0
            monteCount = monteCount + 1;
        end
    end

    monte = SurfaceOfImg * (monteCount/nPoints);
    
    t = toc;

    figure(1);
    hold on;
    subplot(2,3,4), imshow(Img);
    hold on;
    plot(I(1,:), I(2,:), '.r', "MarkerSize", 2);
    title("Surface: " + monte + "%  N = 100" + " t = " + t);
end