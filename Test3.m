clear all 
clc

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


Coupling = 1;

A = {A1,A2};

G1=A1;
G2=A2;
N=size(A1,1);
r=size(A,2);
%%%%%%%%%%%%%%%%%%%%%
Gs=[];
for ii=1:r
    if norm(mod(A{1,ii},1))==0
    As=A{1,ii}*ones(N,1);
    else 
        if nargin~=3
            disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
            return
        else
     epsilon=varargin{1};
     tolu=epsilon/max(abs(As(:)));
      [Au,IA,As] = unique(As,tolu); % As indicates which sums of row are equal w.r.t. the given tolerance
        end
    end
    Gs=[Gs As];
end

IC=[];
for ii=1:N
    index=[];
    for jj=1:r
        temp=(Gs(ii,jj)==Gs(:,jj));
        index=[index temp];
    end
    IC=[IC prod(index,2)];
end
IC=unique(IC','rows');
IC=IC';
IC=IC*[1:size(IC,2)]';
    