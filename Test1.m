% Test program for compute_manifolds program 
% Initial 
close all
clear all 
clc

% Test Examples 

A1=[0 0 1 0 0 1 0 0;
    0 0 0 0 0 0 0 1;
    1 1 0 0 0 0 1 0;
    1 1 1 0 0 0 0 0;
    1 1 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0; 
    1 1 1 0 0 0 0 0;
    0 0 0 0 1 0 0 0];

A2=[0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 1 0;
    0 0 1 1 0 0 1 0;
    0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0];

A = {A1,A2};
% Test Invasive case and A as a cell 
Coupling = 1;
% A1 
tic
[manifolds1 partitions1]= compute_manifolds(A,Coupling); 
toc
tic
[manifolds2 partitions2]= compute_manifolds_original(A,Coupling);
toc
%unique(manifolds1-manifolds2)

%{
% Test Invasive cases 
Coupling = 1;
% A1 
tic
[manifolds1 partitions1]= compute_manifolds(A1,Coupling); 
toc
tic
[manifolds2 partitions2]= compute_manifolds_original(A1,Coupling);
toc
unique(manifolds1-manifolds2)
% A2 
tic
[manifolds1 partitions1]= compute_manifolds(A1,Coupling); 
toc
tic
[manifolds2 partitions2]= compute_manifolds_original(A1,Coupling);
toc
unique(manifolds1-manifolds2)

%  Test Non-invasive cases 
Coupling = 2; 

% A1 
tic
[manifolds1 partitions1]= compute_manifolds(A1,Coupling); 
toc
tic
[manifolds2 partitions2]= compute_manifolds_original(A1,Coupling);
toc
unique(manifolds1-manifolds2)
% A2 
tic
[manifolds1 partitions1]= compute_manifolds(A2,Coupling); 
toc
tic
[manifolds2 partitions2]= compute_manifolds_original(A2,Coupling);
toc
unique(manifolds1-manifolds2)


% Test Real Adjacency 

Ar = [0          0         1.9968 0         0         1.9976 0         0        ;
         0          0         0         0         0         0         1.0254 1.9892;
         1.9989  1.9922 0         0.9933 0         0         1.9868  0       ;
         1.9844  1.9894 2.0016 0         0         0         1.0013  0       ;
         2.0192  1.9823 0.9972 0.9922 0         0          0.9856 0       ;
         0           0        0         0         0         0          2.0130 0       ;
         1.9935  1.9895 1.9885 1.0055 0         0          0         0       ;
         0           0         0         0        2.0091 0          0         0       ;];
Coupling =1; 
%tol = 0.01;
%tol = 0.02;
%tol = 0.022;
tol = 0.025;
tic
[manifolds1 partitions1]= compute_manifolds(Ar,Coupling,tol); 
toc
tic
[manifolds2 partitions2]= compute_manifolds_original(Ar,Coupling,tol);
toc
unique(manifolds1-manifolds2)
%}