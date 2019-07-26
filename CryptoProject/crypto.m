%closing and clearing all open screens
close all;
clear all;

%reading and opening the image
image = imread('D:\Shreyas\CryptoProject\lena.jpg');
figure
imshow(image)
figure
imhist(image)
%taking size of the image
[a,b] = size(image);

%initiating shuffled image matrix
shuffled = zeros(a,b, 'uint8');
twobytwo = zeros(2,2);

%number of iterations of catmap
num_iter = input('Enter the number of arnold cat map iterations: ');

%arnold cat map
%Defining matrix A = [1 p; p pq+1]
p = input('Enter the value of p in arnold cat map: ');
q = input('Enter the value of q in arnold cat map: ');
matrixA = [1 p; p (p*q+1)];

%initiate A^n as I
%eye(2,2) is identity matrix of 2X2 dimensions
matrixApowerN = eye(2,2);

% We are taking num_iter as input and calculating A^num_iter
for j = 1:num_iter
    matrixApowerN = matrixApowerN * matrixA;
end

%We take modulus 256 on each element of matrix matrixApowerN
for k = 1:2
    for l = 1:2
        matrixApowerN(k,l) = mod (matrixApowerN(k,l), 256);
    end
end

%We have now obtained the matrix M, that is, matrixApowerN
%For each pixel, [x' y'] = M[x y] mod256
%Calculating shuffled image pixels using above relation minus the mod256

for m = 1:a
    for n =1:b
        twobytwo = mod(matrixApowerN * [m; n], 256);
        shuffled(twobytwo(1) + 1, twobytwo(2) + 1) = image (m, n);
    end
end

%lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
chaotic = lorentz(a*b,0,35,28,3,1,1,1,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:a*b
    chaotic(i) = floor(mod( (abs( chaotic(i) - floor(chaotic(i)) ) * 1e14), 256 ));
end

%making shuffled matrix as 1d matrix and naming as shuffled2
shuffled2 = shuffled(:);
shuffled2 = uint8(shuffled2);

%converting decimal to binary
shuffled2 = de2bi(shuffled2);
chaotic = de2bi(chaotic);

%encrypted matrix = shuffled XOR key sequence
encrypted_matrix = xor(shuffled2, chaotic);

%to make it into a 1D matrix from 2D of [65536,8]
reshape(encrypted_matrix,[65536*8, 1]);

%converting the 65536*8 binary matrix into a*b matrix after intialisation
encrypted_image = zeros(a, b, 'uint8');
count = 1;
for index = 1 : 8 : 256*256*8
  substring = encrypted_matrix(index:index+7);
  substring = num2str(substring);
  encrypted_image(count) = bin2dec(substring);
  count = count+1;
end
whos encrypted_image
%Converting matrix into grayscale
encrypted_image = mat2gray(encrypted_image);
figure
imshow(encrypted_image, [])
figure
imhist(encrypted_image)

