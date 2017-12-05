function [manifolds,partitions]=compute_manifolds(G,couplingtype,varargin);

% manifolds: all possible partial synchronization manifolds
% partitions: all possible partitions of the systems

% use: 
% 1) adjacency matrix consisting of integer elments:
%
% manifolds=compute_manifolds(G,couplingtype);
% G: adjacency matrix 
% coupling type: 1 (invasive coupling) or 2 (non-invasive coupling)
%
% 
% 2) real adjacency matrix not consisting of integer elments:
%
% manifolds=compute_manifolds(G,couplingtype,tolerance);
% G: adjacency matrix 
% coupling type: 1 (invasive coupling) or 2 (non-invasive coupling)
% tolerance: tolerance in comparison of row-sums (finite precision
% arithmetic)
% 
%
%
% 3) adjacency matrix that can be decomposed as
%   sum_{i=1}^m G_i r_i   with {r1,... r_m} rationally independent
%    and G_i consisting of positive integer elments, i=1,...,m
%
% manifolds=compute_manifolds(G,couplingtype);
% G: structured array of length m; G{i}=G_i 
% coupling type: 1 (invasive coupling) or 2 (non-invasive coupling)
% tolerance: tolerance in comparison of row-sums
% 
%
% Based on the the input type of G, the function automatically
% decides between case 3) and cases 1)-2)






manifolds=[];

%
%
if iscell(G)==1, 
    N=length(G{1});
else N=length(G);
end

% part 1: create all combinations

A=[0;1];

for k=1:N-2;
   B=[];
   for l=1:size(A,1);
       for m=0:max(A(l,:))+1, 
       B=[B; [A(l,:) m]  ];
       end
   end
   A=B;
end

%disp('list of all partitions')

A=[zeros(size(A,1),1) A];
partitions=A;
A=A(1:end-1,:);

%disp('click to continue')
%pause

% part 2: check manifolds - conditions with row sums





if ~iscell(G),
N=length(G);






if norm(mod(G,1))==0,
disp('integer adjacency matrix')





manifolds=[];

if couplingtype==1;
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsum(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end

elseif couplingtype==2;
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsum(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end



else disp('invalid coupling type number')
return
end
    










else % non-integer row sum

if nargin~=3
disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
return
else
epsilon=varargin{1};



manifolds=[];

if couplingtype==1;
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsumnint(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end

elseif couplingtype==2;
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsumnint(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end



else disp('invalid coupling type number')
end
    


end
end








    
 














%
% iscell loop
else

mm=length(G);


manifolds=[];

if couplingtype==1;
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    %%
    for kk=1:mm,
          Gs{kk}=EE*G{kk}*EE';
    end
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
		 %%	            
		  for kk=1:mm,	
		        block= Gs{kk}(positionsh,positionsv); 
                blocksum=rowsum(block);
                 if blocksum==0,
                 isman=0;
                end         
            end
            %%
            %%
            if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end

elseif couplingtype==2;
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrices
    %%
    for kk=1:mm,
         Gs{kk}=EE*G{kk}*EE';
    end
    %%
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
		 %%
            for kk=1:mm,
                 block= Gs{kk}(positionsh,positionsv); 
                 blocksum=rowsum(block);
                 if blocksum==0,
                    isman=0;
                 end
            end
            %%         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end



else disp('invalid coupling type number')
return
end



























%
% iscell loop
%
end




















   
    
% inline functions    
    
function y=rowsum(bl);
    [nn]=size(bl);
    bl2=bl*ones(nn(2),1);
    bl2=abs(bl2-bl2(1)*ones(nn(1),1));
   if max(bl2)==0,
       y=1;
   else
       y=0;
   end 
end
    

function y=rowsumnint(bl);
    [nn]=size(bl);
    bl2=bl*ones(nn(2),1);
    bl2=abs(bl2-bl2(1)*ones(nn(1),1));
   if max(bl2)<epsilon,
       y=1;
   else
       y=0;
   end 
end





% end main function
end