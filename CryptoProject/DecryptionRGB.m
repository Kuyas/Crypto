%------------Encrypted Image Read------------%
%Closing and clearing all open screens
close all;
clear all;
eg = imread("house.tiff");
key = imread("key.tiff");

%--------------------------------------------------------------------%
%------------------------XOR WITH KEY IMAGE START--------------------%
%--------------------------------------------------------------------%

eg = de2bi(chaotic);
key = de2bi(key);

OIG = bitxor(eg,key);


%--------------------------------------------------------------------%
%------------------------XOR WITH KEY IMAGE END----------------------%
%--------------------------------------------------------------------%


%--------------------------------------------------------------------%
%------------------------OIG XOR CHAOTIC START-----------------------%
%--------------------------------------------------------------------%
[u,v,w] = size(eg);
%RGB is seperated into 3 matrices
R1 = eg(:,:,1); % Red channel
G1 = eg(:,:,2); % Green channel
B1 = eg(:,:,3); % Blue channel

%----------chaotic map initialize--------------%
%lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
chaotic = lorentz(u*v,0,35,28,3,1,1,1,0.001);
%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v
    chaotic(i) = floor(mod( (abs( chaotic(i) - floor(chaotic(i)) ) * 1e14), 256 ));
end
chaotic = uint8(chaotic);
%-----------------BITXOR WITH CHAOTIC AND R,G,B LAYERS SEPARATELY------%
R1 = de2bi(R1);
G1 = de2bi(G1);
B1 = de2bi(B1);
chaotic = de2bi(chaotic);


z = 0;
while z<3 
    if z == 0
        oig = R1;
    elseif z==1
        oig = G1;
    elseif z == 2
        oig = B1;
    end
        
em = bitxor(oig, chaotic);
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
if z == 0
        R2 = ei;
    elseif z==1
        G2 = ei;
    elseif z == 2
        B2 = ei;
end
z++;
end

shuffle(:,:,1) = R2;
shuffle(:,:,2) = G2;
shuffle(:,:,3) = B2;

%--------------------------------------------------------------------%
%------------------------OIG XOR CHAOTIC END-------------------------%
%--------------------------------------------------------------------%

%--------------------------------------------------------------------%
%------------------------------UNSORT--------------------------------%
%--------------------------------------------------------------------%

%----------chaotic map initialize--------------%
%lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
chaotic = lorentz(u*v*w,0,35,28,3,1,1,1,0.001);
%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:u*v*w
    chaotic(i) = floor(mod( (abs( chaotic(i) - floor(chaotic(i)) ) * 1e14), 256 ));
end
chaotic = uint8(chaotic);
c_sort = reshape(chaotic);
[c d] = sort(chaotic);

%--------------------------------------------------------------------%
%--------------------------------UNSORT------------------------------%
%--------------------------------------------------------------------%



%-------------END-----------------------------$