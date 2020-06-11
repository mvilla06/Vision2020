clc
close all
clear all
I = imread("image_gray.jpg");
figure(1);
imshow(I);
#Promedio
K = [1, 1, 1; 1, 1, 1; 1, 1, 1];
K = K/9;

#Gaussiano

K = [1, 2, 1; 2, 4, 2; 1, 2, 1];
K = K/16;

#Laplaciano
K=[-1, -1, -1; -1, 8, -1; -1, -1, -1];

#Laplaciano de Gaussiano
K = [0, 0, 1, 0, 0;
      0, 1, 2, 1, 0;
      1, 2, -16, 2, 1;
      0, 1, 2, 1, 0;
      0, 0, 1, 0, 0];
      
#Derivativos
#X
K=[1, 0, -1; 2, 0, -2; 1, 0, -1];
#Y
K = [1, 2, 1; 0, 0, 0; -1, -2, -1];

#Enfatizador
K=[0, -1, 0; -1, 5, -1; 0, -1, 0];


R = conv2(I, K, "same");
R = round(R);
R = uint8( R);
figure(2);
imshow(R);