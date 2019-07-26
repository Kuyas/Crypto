%Closing and clearing all open screens
close all;
clear all;

%Reading and opening the image
image = imread('house256.tiff');

figure
imshow(image)

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

%---------------------ARNOLD CAT MAP BEGINS----------------------%

%Defining matrix A of Arnold Cat Map 3D
%a = input('Enter the value of p in arnold cat map: ');
%b = input('Enter the value of q in arnold cat map: ');
%c = input('Enter the value of r in arnold cat map: ');
%d = input('Enter the value of s in arnold cat map: ');
a=5;
b=3;
c=7;
d=2;
matrixA = [1 a b; c ((a*c)+1) (b*c); d (a*b*c*d) ((b*d)+1)];

%Elements of MatrixA mod 256
matrixA = mod(matrixA, 256);

for k = 1:w
    for i = 1:u
        for j = 1:v
            threebyone = mod(matrixA * [i; j; k], 256);
            shuffled(threebyone(1,1) + 1, threebyone(2,1) + 1, mod(threebyone(3,1),w) + 1) = image (i,j,k);
        end
    end
end

%Arnold cat map has been performed on the Image that has been provided
%The resulting shuffled image is stored in 'shuffled'

%----------------------ARNOLD CAT MAP ENDS-----------------------%


%----------------------CHAOTIC MAP-------------------------------%
R1 = shuffled(:,:,1); % Red channel
G1 = shuffled(:,:,2); % Green channel
B1 = shuffled(:,:,3); % Blue channel

%Lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
%Runge-kutta step size is initiated as 0.001.
chaotic = lorentz(u*v*w,0,35,28,3,1,1,1,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v*w
    chaotic(i) = floor(mod((abs( chaotic(i) - floor(chaotic(i))) * 1e14), 256));
end
chaotic = uint8(chaotic);

%Chaotic map of size u*v has been generated in a 1D matrix

%---------------------CHAOTIC MAP ENDS-----------------------%

image = imread('girl2.tiff');

%Initialize the size of the image
[u,v,w] = size(image);
%RGB is seperated into 3 matrices
img_red = image(:,:,1); % Red channel
img_green = image(:,:,2); % Green channel
img_blue = image(:,:,3); % Blue channel

%---------------------SHUFFLING BEGINS-----------------------
image1d = image(:);
[a,b] = sort(chaotic);

for k = 1:u*v*w
    shuffled123(b(k)) = image1d(k);
end


for k = 1:u*v
    shuffled_1(k) = shuffled123(k);
end

for k = u*v + 1:u*v*2
    shuffled_2(k - (u*v)) = shuffled123(k);
end

for k = u*v*2 + 1:u*v*3
    shuffled_3(k - (u*v*2)) = shuffled123(k);
end

shuffled_1 = reshape(shuffled_1, [u,v]);
shuffled_2 = reshape(shuffled_2, [u,v]);
shuffled_3 = reshape(shuffled_3, [u,v]);

shuffledim2 = zeros(u,v,w);
shuffledim2 = uint8(shuffledim2);

shuffledim2(:,:,1) = shuffled_1;
shuffledim2(:,:,2) = shuffled_2;
shuffledim2(:,:,3) = shuffled_3;

%------------------SHUFFLING DONE----------------------

shuffled = de2bi(shuffled);
shuffledim2 = de2bi(shuffledim2);

final = bitxor(shuffled,shuffledim2);

reshape(final,[u*v*w*8,1]);
encrypted_image = zeros(u,v,w);
encrypted_image = uint8(encrypted_image);
count = 1;
for index = 1 : 8 : u*v*w*8
  substring = final(index:index+7);
  substring = num2str(substring);
  encrypted_image(count) = bin2dec(substring);
  count = count+1;
end

%figure
%imshow(encrypted_image)
ei1 = encrypted_image;


%-------------------------------------------------------------------
%---------------------------------------------------------------------%


%Reading and opening the image
image = imread('jelly.tiff');

%figure
%imshow(image)

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

%---------------------ARNOLD CAT MAP BEGINS----------------------%

%Defining matrix A of Arnold Cat Map 3D
%a = input('Enter the value of p in arnold cat map: ');
%b = input('Enter the value of q in arnold cat map: ');
%c = input('Enter the value of r in arnold cat map: ');
%d = input('Enter the value of s in arnold cat map: ');
a=5;
b=3;
c=7;
d=2;
matrixA = [1 a b; c ((a*c)+1) (b*c); d (a*b*c*d) ((b*d)+1)];

%Elements of MatrixA mod 256
matrixA = mod(matrixA, 256);

for k = 1:w
    for i = 1:u
        for j = 1:v
            threebyone = mod(matrixA * [i; j; k], 256);
            shuffled(threebyone(1,1) + 1, threebyone(2,1) + 1, mod(threebyone(3,1),w) + 1) = image (i,j,k);
        end
    end
end

%Arnold cat map has been performed on the Image that has been provided
%The resulting shuffled image is stored in 'shuffled'

%----------------------ARNOLD CAT MAP ENDS-----------------------%


%----------------------CHAOTIC MAP-------------------------------%
R1 = shuffled(:,:,1); % Red channel
G1 = shuffled(:,:,2); % Green channel
B1 = shuffled(:,:,3); % Blue channel

%Lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
%Runge-kutta step size is initiated as 0.001.
chaotic = lorentz(u*v*w,0,35,28,3,1,1,1,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v*w
    chaotic(i) = floor(mod((abs( chaotic(i) - floor(chaotic(i))) * 1e14), 256));
end
chaotic = uint8(chaotic);

%Chaotic map of size u*v has been generated in a 1D matrix

%---------------------CHAOTIC MAP ENDS-----------------------%

image = imread('girl2.tiff');

%Initialize the size of the image
[u,v,w] = size(image);
%RGB is seperated into 3 matrices
img_red = image(:,:,1); % Red channel
img_green = image(:,:,2); % Green channel
img_blue = image(:,:,3); % Blue channel

%---------------------SHUFFLING BEGINS-----------------------
image1d = image(:);
[a,b] = sort(chaotic);

for k = 1:u*v*w
    shuffled123(b(k)) = image1d(k);
end


for k = 1:u*v
    shuffled_1(k) = shuffled123(k);
end

for k = u*v + 1:u*v*2
    shuffled_2(k - (u*v)) = shuffled123(k);
end

for k = u*v*2 + 1:u*v*3
    shuffled_3(k - (u*v*2)) = shuffled123(k);
end

shuffled_1 = reshape(shuffled_1, [u,v]);
shuffled_2 = reshape(shuffled_2, [u,v]);
shuffled_3 = reshape(shuffled_3, [u,v]);

shuffledim2 = zeros(u,v,w);
shuffledim2 = uint8(shuffledim2);

shuffledim2(:,:,1) = shuffled_1;
shuffledim2(:,:,2) = shuffled_2;
shuffledim2(:,:,3) = shuffled_3;

%------------------SHUFFLING DONE----------------------

shuffled = de2bi(shuffled);
shuffledim2 = de2bi(shuffledim2);

final = bitxor(shuffled,shuffledim2);

reshape(final,[u*v*w*8,1]);
encrypted_image = zeros(u,v,w);
encrypted_image = uint8(encrypted_image);
count = 1;
for index = 1 : 8 : u*v*w*8
  substring = final(index:index+7);
  substring = num2str(substring);
  encrypted_image(count) = bin2dec(substring);
  count = count+1;
end

ei2 = encrypted_image;

results = NPCR_and_UACI( ei1, ei2, 1, 255 )

%---------------------DECRYPTION-------------------------%
%----ei1 = encrypted image 1; ei2 = encrypted image 2----%
%--- Step 1 XOR ( Self invertible )
%--- Step 2 Arnold Inverse sort reverse ---------------------%

%ei1 = ei1(:);
%shuffledim2 = shuffledim2(:);
ei1 = de2bi(ei1);
%shuffledim2 = de2bi(shuffledim2);

decryptArnold = bitxor(ei1,shuffledim2);

reshape(decryptArnold,[u*v*w*8,1]);
decrypted_image = zeros(u,v,w);
decrypted_image = uint8(decrypted_image);
count = 1;
for index = 1 : 8 : u*v*w*8
  substring = decryptArnold(index:index+7);
  substring = num2str(substring);
  decrypted_image(count) = bin2dec(substring);
  count = count+1;
end




%{
twobytwo = zeros(2,2);
threebyone = zeros(3,1);
shuffled = zeros(size(image), 'uint8');

decryptArnold1 = decrypted_image;
a=5;
b=3;
c=7;
d=2;
matrixA = [1 a b; c ((a*c)+1) (b*c); d (a*b*c*d) ((b*d)+1)];
matrixB = inv(matrixA);
matrixB = ceil(matrixB);
matrixB = mod(matrixB, 256);

image = decryptArnold1;
[u,v,w] = size(image);
for k = 1:w
    for i = 1:u
        for j = 1:v
            threebyone = mod(matrixB * [i; j; k], 256);
            shuffled(threebyone(1,1) + 1, threebyone(2,1) + 1, mod(threebyone(3,1),w) + 1) = image (i,j,k);
        end
    end
end
%}
[u,v,w] = size(decrypted_image)
iteration = 1;
% Initialize image.
oldScrambledImage = decrypted_image;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
N = u;
while iteration <= 3 * N
	% Scramble the image based on the old image.
    for k = 1:w
    	for row = 1 : v % y
        	for col = 1 : u % x
            	c = mod((2 * col) + row, N) + 1; % x coordinate
                r = mod(col + row, N) + 1; % y coordinate
                % Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
                currentScrambledImage(row, col, :) = oldScrambledImage(r, c, :);
            end
    	end
    end
    
	% Display the current image.
	caption = sprintf('Iteration #%d', iteration);
	fprintf('%s\n', caption);
	subplot(1, 2, 2);
	imshow(currentScrambledImage);
	axis on;
	title(caption, 'FontSize', fontSize);
	drawnow;
	
	% Insert a delay if desired.
% 	pause(0.1);
	
	% Save the image, if desired.
	filename = sprintf('Arnold Cat Iteration %d.png', iteration);
% 	imwrite(currentScrambledImage, filename);
	fprintf('Saved image file %s to disk.\n', filename);
	
%	if immse(currentScrambledImage, grayImage) == 0
%		caption = sprintf('Back to Original after %d Iterations.', iteration);
%		fprintf('%s\n', caption);
%		title(caption, 'FontSize', fontSize);
%		break;
%	end
	
	% Make the current image the prior/old one so we'll operate on that the next iteration.
	oldScrambledImage = currentScrambledImage;
	% Update the iteration counter.
	iteration = iteration+1;
end

figure
imshow(image)

