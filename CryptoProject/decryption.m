%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%DECRYPTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of iterations of catmap
num_iter = 5
%arnold cat map
%Defining matrix A = [1 p; p pq+1]
p = 1
q = 1
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
%%%%%%%%%%%%%
%ENCRYPTION P1 = q%%%
imq = zeros(256,256);
[a,b] = size(imq);

%initiating shuffled image matrix
shuffledq = zeros(a,b);
twobytwoq = zeros(2,2);

for m = 1:a
    for n =1:b
        twobytwoq = mod(matrixApowerN * [m; n], 256);
        shuffledq(twobytwoq(1) + 1, twobytwoq(2) + 1) = imq (m, n);
    end
end

%lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
chaotic = lorentz(a*b,0,35,28,3,1,1,1,0.001);
for i = 1:a*b
    chaotic(i) = floor(mod( (abs( chaotic(i) - floor(chaotic(i)) ) * 1e14), 256 ));
end
%making shuffled matrix as 1d matrix and naming as p1q
plq = shuffledq(:);
plq = uint8(plq);

%converting decimal to binary
plq = de2bi(plq);
chaotic = de2bi(chaotic);


ep1 = xor(plq, chaotic);
reshape(ep1,[65536*8, 1]);
c1 = zeros(a,b);

count = 1;
for index = 1 : 8 : 256*256*8
  substring = ep1(index:index+7);
  substring = num2str(substring);
  c1(count) = bin2dec(substring);
  count = count+1;
end
k = bitxor(c1,imq);

disp(k);

tf = isequal(k,chaotic);
disp(chaotic);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%ENCRYPTION P2 = p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imp = zeros(256,256);
imp(241,1) = 1;
imp(241,6) = 2;
[a,b] = size(imp);
s_p = zeros(a,b);
t_p = zeros(2,2);

for m = 1:a
    for n =1:b
        t_p = mod(matrixApowerN * [m; n], 256);
        s_p(t_p(1) + 1, t_p(2) + 1) = imp (m, n);
    end
end
%disp(s_p);

%lorentz chaotic mapping solutions lorentz(n,level,s,r,b,x0,y0,z0,h)
chaotic = lorentz(a*b,0,35,28,3,1,1,1,0.001);
%Obtaining final sequence K with ( |x - floor(x)| X 10^14 )mod256
for i = 1:a*b
    chaotic(i) = floor(mod( (abs( chaotic(i) - floor(chaotic(i)) ) * 1e14), 256 ));
end
%making shuffled matrix as 1d matrix and naming as p2
p2 = s_p(:);
p2 = uint8(p2);
pn = zeros(65536,1);

pn(:,1) = p2(:,1);

%chaotic = uint8(chaotic);
%disp(pn);
pn = uint8(pn);
%converting decimal to binary
pn = de2bi(pn);
chaotic = de2bi(chaotic);
pn = pn(:,1);
ep2 = xor(pn, chaotic);

%to make it into a 1D matrix from 2D of [65536,8]
reshape(ep2,[65536*8, 1]);

%converting the 65536*8 binary matrix into a*b matrix after intialisation
c2 = zeros(a, b);

count = 1;
for index = 1 : 8 : 256*256*8
  substring = ep2(index:index+7);
  substring = num2str(substring);
  c2(count) = bin2dec(substring);
  count = count+1;
end

s2 = bitxor(c2,k);

%disp(s2);


