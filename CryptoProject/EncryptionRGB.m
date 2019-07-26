%Closing and clearing all open screens
close all;
clear all;

%Reading and opening the image
image = imread('house.tiff');

%figure
%imshow(image)
%figure
%imhist(image)

%Initialize the size of the image
[u,v,w] = size(image);

%RGB is seperated into 3 matrices
img_red = image(:,:,1); % Red channel
img_green = image(:,:,2); % Green channel
img_blue = image(:,:,3); % Blue channel

%Initiating the matrices required for Arnold Cat Map
shuffled_red = zeros(u,v, 'uint8');
shuffled_green = zeros(u,v, 'uint8');
shuffled_blue = zeros(u,v, 'uint8');
twobytwo = zeros(2,2);
threebyone = zeros(3,1);
shuffled = zeros(size(image), 'uint8');

%---------------------ARNOLD CAT MAP BEGINS----------------------

%Defining matrix A of Arnold Cat Map 3D
a = input('Enter the value of p in arnold cat map: ');
b = input('Enter the value of q in arnold cat map: ');
c = input('Enter the value of r in arnold cat map: ');
d = input('Enter the value of s in arnold cat map: ');
matrixA = [1 a b; c ((a*c)+1) (b*c); d (a*b*c*d) ((b*d)+1)];

%Elements of MatrixA mod 256
matrixA = mod(matrixA, 256);

matrixA = matrixA*matrixA;
matrixB = inv(matrixA);
matrixB = abs(matrixB);
%C = matrixA*matrixB;
%C
%whos matrixA matrixB

for k = 1:w
    for i = 1:u
        for j = 1:v
            threebyone = mod(matrixA * [i; j; k], 256);
            shuffled(threebyone(1,1) + 1, threebyone(2,1) + 1, mod(threebyone(3,1),w) + 1) = image (i,j,k);
        end
    end
end


%threebyone
%figure
%imshow(shuffled)
%pp = primpoly(16,'all');
%pp
%g_f = pp;
%g_f = gf(pp,16);
%gf
%gf = mod(gf,256);
%----------------------ARNOLD CAT MAP ENDS-----------------------%
%----------------------CHAOTIC MAP INITIALIZE-------------------------------%
R1 = shuffled(:,:,1); % Red channel
G1 = shuffled(:,:,2); % Green channel
B1 = shuffled(:,:,3); % Blue channel

%lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
chaotic = lorentz(u*v,0,35,28,3,1,1,1,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v
    chaotic(i) = floor(mod( (abs( chaotic(i) - floor(chaotic(i)) ) * 1e14), 256 ));
end
chaotic = uint8(chaotic);
%whos R1 chaotic
%----------------------CHAOTIC MAP INITIALIZING DONE-------------------%
%-----------------BITXOR WITH CHAOTIC AND R,G,B LAYERS SEPARATELY------%
R1 = de2bi(R1);
G1 = de2bi(G1);
B1 = de2bi(B1);
chaotic = de2bi(chaotic);

em = bitxor(R1, chaotic);
%to make it into a 1D matrix from 2D of [65536,8]
reshape(em,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = em(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

R2 = ei;

em = bitxor(G1, chaotic);
%to make it into a 1D matrix from 2D of [65536,8]
reshape(em,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = em(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

G2 = ei;
em = bitxor(B1, chaotic);
%to make it into a 1D matrix from 2D of [65536,8]
reshape(em,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = em(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

B2 = ei;
%-----------------BITXOR WITH CHAOTIC AND R,G,B LAYERS SEPARATELY DONE------%
%--------------COLUMN INTERLEAVING----------------------%
col = reshape([R2(:) G2(:) B2(:)]',3*size(R2,1), [])';
R2 = col(:,1:256);
G2 = col(:,257:512);
B2 = col(:,513:768);
%--------------COLUMN INTERLEAVING DONE----------------%
%--------------ROW INTERLEAVING---------------------%
row = reshape([R2(:) G2(:) B2(:)]',3*size(R2,1), []);
R2 = row(1:256,:);
G2 = row(257:512,:);
B2 = row(513:768,:);
%---------------ROW INTERLEAVING DONE---------------%

%----------------COMBINING THE LAYERS TO FORM A COLOUR IMAGE-------------%
%whos R2 G2 B2
I(:,:,1) = R2;
I(:,:,2) = G2;
I(:,:,3) = B2;

%figure
%imshow(I);
%figure
%imhist(I);

%---------------Image Encryption done until chaotic----------------%
%---------------Key Image Encryption Begin-------------------------%

%Reading and opening the image
image = imread('girl.jpg');

%figure
%imshow(image)
%figure
%imhist(image)

%Initialize the size of the image
[u,v,w] = size(image);

%RGB is seperated into 3 matrices
img_red = image(:,:,1); % Red channel
img_green = image(:,:,2); % Green channel
img_blue = image(:,:,3); % Blue channel

%Initiating the matrices required for Arnold Cat Map
shuffled_red = zeros(u,v, 'uint8');
shuffled_green = zeros(u,v, 'uint8');
shuffled_blue = zeros(u,v, 'uint8');
twobytwo = zeros(2,2);
threebyone = zeros(3,1);
shuffled = zeros(size(image), 'uint8');

%---------------------ARNOLD CAT MAP BEGINS----------------------

%Defining matrix A of Arnold Cat Map 3D
%a = input('Enter the value of p in arnold cat map: ');
%b = input('Enter the value of q in arnold cat map: ');
%c = input('Enter the value of r in arnold cat map: ');
%d = input('Enter the value of s in arnold cat map: ');
matrixA = [1 a b; c ((a*c)+1) (b*c); d (a*b*c*d) ((b*d)+1)];

%Elements of MatrixA mod 256
matrixA = mod(matrixA, 256);

matrixA = matrixA*matrixA;
matrixB = inv(matrixA);
matrixB = abs(matrixB);
%C = matrixA*matrixB;
%C
%whos matrixA matrixB

for k = 1:w
    for i = 1:u
        for j = 1:v
            threebyone = mod(matrixA * [i; j; k], 256);
            shuffled(threebyone(1,1) + 1, threebyone(2,1) + 1, mod(threebyone(3,1),w) + 1) = image (i,j,k);
        end
    end
end


%threebyone
%figure
%imshow(shuffled)
%pp = primpoly(16,'all');
%pp
%g_f = pp;
%g_f = gf(pp,16);
%gf
%gf = mod(gf,256);
%----------------------ARNOLD CAT MAP ENDS-----------------------%
%----------------------CHAOTIC MAP INITIALIZE-------------------------------%
R1 = shuffled(:,:,1); % Red channel
G1 = shuffled(:,:,2); % Green channel
B1 = shuffled(:,:,3); % Blue channel

%lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
chaotic = lorentz(u*v,0,35,28,3,1,1,1,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v
    chaotic(i) = floor(mod( (abs( chaotic(i) - floor(chaotic(i)) ) * 1e14), 256 ));
end
chaotic = uint8(chaotic);
%whos R1 chaotic
%----------------------CHAOTIC MAP INITIALIZING DONE-------------------%
%-----------------BITXOR WITH CHAOTIC AND R,G,B LAYERS SEPARATELY------%
R1 = de2bi(R1);
G1 = de2bi(G1);
B1 = de2bi(B1);
chaotic = de2bi(chaotic);

em = bitxor(R1, chaotic);
%to make it into a 1D matrix from 2D of [65536,8]
reshape(em,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = em(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

R2 = ei;

em = bitxor(G1, chaotic);
%to make it into a 1D matrix from 2D of [65536,8]
reshape(em,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = em(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

G2 = ei;
em = bitxor(B1, chaotic);
%to make it into a 1D matrix from 2D of [65536,8]
reshape(em,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = em(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

B2 = ei;
%-----------------BITXOR WITH CHAOTIC AND R,G,B LAYERS SEPARATELY DONE------%
%--------------COLUMN INTERLEAVING----------------------%
col = reshape([R2(:) G2(:) B2(:)]',3*size(R2,1), [])';
R2 = col(:,1:256);
G2 = col(:,257:512);
B2 = col(:,513:768);
%--------------COLUMN INTERLEAVING DONE----------------%
%--------------ROW INTERLEAVING---------------------%
row = reshape([R2(:) G2(:) B2(:)]',3*size(R2,1), []);
R2 = row(1:256,:);
G2 = row(257:512,:);
B2 = row(513:768,:);
%---------------ROW INTERLEAVING DONE---------------%

%----------------COMBINING THE LAYERS TO FORM A COLOUR IMAGE-------------%
%whos R2 G2 B2
K(:,:,1) = R2;
K(:,:,2) = G2;
K(:,:,3) = B2;

%figure
%imshow(I);
%figure
%imhist(I);

%---------------Key Image Encryption done until chaotic----------------%
%---------------Bitxor between I & K---------------------------------------%
R1 = I(:,:,1); % Red channel
G1 = I(:,:,2); % Green channel
B1 = I(:,:,3); % Blue channel

R2 = K(:,:,1); % Red channel
G2 = K(:,:,2); % Green channel
B2 = K(:,:,3); % Blue channel

R1 = de2bi(R1);
G1 = de2bi(G1);
B1 = de2bi(B1);

R2 = de2bi(R2);
G2 = de2bi(G2);
B2 = de2bi(B2);

R_ei = bitxor(R1,R2);
G_ei = bitxor(G1,G2);
B_ei = bitxor(B1,B2);

%to make it into a 1D matrix from 2D of [65536,8]
reshape(R_ei,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = R_ei(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

R_f = ei;

%to make it into a 1D matrix from 2D of [65536,8]
reshape(G_ei,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = G_ei(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

G_f = ei;

%to make it into a 1D matrix from 2D of [65536,8]
reshape(B_ei,[65536*8, 1]);
%converting the 65536*8 binary matrix into a*b matrix after intialisation
ei = zeros(u, v, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = B_ei(index:index+7);
  substring = num2str(substring);
  ei(count) = bin2dec(substring);
  count = count+1;
end

B_f = ei;

encryptedImage(:,:,1) = R_f;
encryptedImage(:,:,2) = G_f;
encryptedImage(:,:,3) = B_f;

%----------------------------------------------------------------------%
%whos K I
figure
imshow(encryptedImage)
figure
imhist(encryptedImage)
