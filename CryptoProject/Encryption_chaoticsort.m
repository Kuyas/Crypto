%Closing and clearing all open screens
close all;
clear all;

image = imread('house256.tiff');

%Initialize the size of the image
[u,v,w] = size(image);
%RGB is seperated into 3 matrices
img_red = image(:,:,1); % Red channel
img_green = image(:,:,2); % Green channel
img_blue = image(:,:,3); % Blue channel

%----------------------CHAOTIC MAP-------------------------------%

%Lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
%Runge-kutta step size is initiated as 0.001.
chaotic = lorentz(u*v*w,0,35,28,3,22,32,18,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v*w
    chaotic(i) = floor(mod((abs( chaotic(i) - floor(chaotic(i))) * 1e14), 256));
end
chaotic = uint8(chaotic);

%Chaotic map of size u*v has been generated in a 1D matrix
%---------------------CHAOTIC MAP ENDS-----------------------%
%---------------------SHUFFLING BEGINS-----------------------
[a,b] = sort(chaotic);
image1d = image(:);

for k = 1:u*v*w
    shuffled(b(k)) = image1d(k);
end


for k = 1:u*v
    shuffled_1(k) = shuffled(k);
end

for k = u*v + 1:u*v*2
    shuffled_2(k - (u*v)) = shuffled(k);
end

for k = u*v*2 + 1:u*v*3
    shuffled_3(k - (u*v*2)) = shuffled(k);
end

shuffled_1 = reshape(shuffled_1, [u,v]);
shuffled_2 = reshape(shuffled_2, [u,v]);
shuffled_3 = reshape(shuffled_3, [u,v]);

shuffledim1 = zeros(u,v,w);
shuffledim1 = uint8(shuffledim1);

shuffledim1(:,:,1) = shuffled_1;
shuffledim1(:,:,2) = shuffled_2;
shuffledim1(:,:,3) = shuffled_3;

%------------------SHUFFLING DONE----------------------

image = imread('jelly.tiff');

%Initialize the size of the image
[u,v,w] = size(image);
%RGB is seperated into 3 matrices
img_red = image(:,:,1); % Red channel
img_green = image(:,:,2); % Green channel
img_blue = image(:,:,3); % Blue channel

%----------------------CHAOTIC MAP-------------------------------%

%Lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
%Runge-kutta step size is initiated as 0.001.
chaotic = lorentz(u*v*w,0,35,28,3,22,32,18,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v*w
    chaotic(i) = floor(mod((abs( chaotic(i) - floor(chaotic(i))) * 1e14), 256));
end
chaotic = uint8(chaotic);

%Chaotic map of size u*v has been generated in a 1D matrix
%---------------------CHAOTIC MAP ENDS-----------------------%
%---------------------SHUFFLING BEGINS-----------------------
[a,b] = sort(chaotic);
image1d = image(:);

for k = 1:u*v*w
    shuffled(b(k)) = image1d(k);
end


for k = 1:u*v
    shuffled_1(k) = shuffled(k);
end

for k = u*v + 1:u*v*2
    shuffled_2(k - (u*v)) = shuffled(k);
end

for k = u*v*2 + 1:u*v*3
    shuffled_3(k - (u*v*2)) = shuffled(k);
end

shuffled_1 = reshape(shuffled_1, [u,v]);
shuffled_2 = reshape(shuffled_2, [u,v]);
shuffled_3 = reshape(shuffled_3, [u,v]);

shuffledim1 = zeros(u,v,w);
shuffledim1 = uint8(shuffledim1);

shuffledim3(:,:,1) = shuffled_1;
shuffledim3(:,:,2) = shuffled_2;
shuffledim3(:,:,3) = shuffled_3;

%------------------SHUFFLING DONE----------------------



image = imread('girl2.tiff');

%Initialize the size of the image
[u,v,w] = size(image);
%RGB is seperated into 3 matrices
img_red = image(:,:,1); % Red channel
img_green = image(:,:,2); % Green channel
img_blue = image(:,:,3); % Blue channel

%Lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
%Runge-kutta step size is initiated as 0.001.
chaotic = lorentz(u*v*w,0,35,28,3,2,2,8,0.001);

%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v*w
    chaotic(i) = floor(mod((abs( chaotic(i) - floor(chaotic(i))) * 1e14), 256));
end
chaotic = uint8(chaotic);

%Chaotic map of size u*v has been generated in a 1D matrix
[a,b] = sort(chaotic);

%---------------------SHUFFLING BEGINS-----------------------
image1d = image(:);

for k = 1:u*v*w
    shuffled(b(k)) = image1d(k);
end


for k = 1:u*v
    shuffled_1(k) = shuffled(k);
end

for k = u*v + 1:u*v*2
    shuffled_2(k - (u*v)) = shuffled(k);
end

for k = u*v*2 + 1:u*v*3
    shuffled_3(k - (u*v*2)) = shuffled(k);
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


shuffledim1 = de2bi(shuffledim1);
shuffledim2 = de2bi(shuffledim2);
shuffledim3 = de2bi(shuffledim3);

final = bitxor(shuffledim1,shuffledim2);

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
ei1 = encrypted_image;

figure
imshow(ei1)

final = bitxor(shuffledim3,shuffledim2);

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
