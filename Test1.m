clc
clear all

% Test Example from the Article 

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
Coupling = 1;
%G = A1;
%Grow = G*ones(size(G,2),1);
%    [Grows,IA,IC] = unique(Grow)
At = A1;
At = A2;
tic
[manifolds1 partitions1]= compute_manifolds_new(At,Coupling); 
toc
tic
[manifolds2 partitions2]= compute_manifolds(At,Coupling);
toc
Coupling = 2; 
%manifolds = compute_manifolds_mod(A,Coupling);