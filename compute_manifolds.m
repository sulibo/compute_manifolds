function [PS,manifolds,partitions]=compute_manifolds(G,couplingtype,F,varargin)
% Version 3.0
% manifolds: all possible partial synchronization manifolds
% partitions: all possible partitions of the systems
% PS: array of structures with - manifolds 
%                                           - R: reordering matrix 
%                                           - Ared/Lred: synchronization error dynamics 
%                                              related matrices (Ared for invasive coupling, Lred
%                                              for non-invasive coupling)
%                                           - Asyn: adjancency matrix of the reduced
%                                              network 
% 

% use: 
% 1) adjacency matrix consisting of integer elments:
%
% manifolds=compute_manifolds(G,couplingtype,F);
% G: adjacency matrix 
% coupling type: 1 (invasive coupling) or 2 (non-invasive coupling)
% F: indicator of dyanmics types
% 
% 2) real adjacency matrix not consisting of integer elments:
%
% manifolds=compute_manifolds(G,couplingtype,tolerance,F);
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
%
% 4) F indicates the type of dynamics in each node



%% Initialization
manifolds=[];
partitions=[];

% Check system dynamics differences 
 if length(unique(F))==length(F)
    disp('All system dynamics are different. No partial synchronization manifolds exist.')
    PS=[];
    return
 end
% Dimensions
if iscell(G)==1 
    N=length(G{1});
else
    N=length(G);
end

%% part 1: create all combinations
% G is a cell 
if couplingtype ==1
    if nargin ==4
        tolerance = varargin{1};
    [A, partitions]=generate_partition_invasive(G,F, tolerance);
    else 
    [A, partitions]=generate_partition_invasive(G,F);
    end
else
    [A, partitions]=generate_partition_universal(G,F);
end

%disp('click to continue')
%pause

%% part 2: check manifolds - conditions with row sums
R={}; %cell for storage of reorder matrix
if ~iscell(G)
N=length(G);

if norm(mod(G,1))==0
disp('integer adjacency matrix')

manifolds=[];

if couplingtype==1
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs)
        if rijs(m)~=rijs(m-1)
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1)
        for tel2=1:(length(indices)-1)
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsum(block);
            if blocksum==0
               isman=0;
            end         
           if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1
        manifolds=[manifolds;rij]; 
        R{end+1}=EE';
    end
end

elseif couplingtype==2
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs)
        if rijs(m)~=rijs(m-1)
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1)
        for tel2=1:(length(indices)-1)
          if tel1~=tel2
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsum(block);
            if blocksum==0
               isman=0;
            end         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1
        manifolds=[manifolds;rij]; 
        R{end+1}=EE';
    end
end



else disp('invalid coupling type number')
return
end
    



else % non-integer row sum

if nargin~=4
disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
return
else
epsilon=varargin{1};



manifolds=[];

if couplingtype==1
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs)
        if rijs(m)~=rijs(m-1)
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1)
        for tel2=1:(length(indices)-1)
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsumnint(block);
            if blocksum==0
               isman=0;
            end         
           if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1
        manifolds=[manifolds;rij]; 
        R{end+1}=EE';
    end
end

elseif couplingtype==2
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs)
        if rijs(m)~=rijs(m-1)
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1)
        for tel2=1:(length(indices)-1)
          if tel1~=tel2
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsumnint(block);
            if blocksum==0
               isman=0;
            end         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1
        manifolds=[manifolds;rij]; 
        R{end+1}=EE';
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

if couplingtype==1
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    %%
    for kk=1:mm
          Gs{kk}=EE*G{kk}*EE';
    end
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs)
        if rijs(m)~=rijs(m-1)
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1)
        for tel2=1:(length(indices)-1)
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
		 %%	            
		  for kk=1:mm
		        block= Gs{kk}(positionsh,positionsv); 
                blocksum=rowsum(block);
                 if blocksum==0
                 isman=0;
                end         
            end
            %%
            %%
            if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1
        manifolds=[manifolds;rij]; 
        R{end+1}=EE';
    end
end

elseif couplingtype==2
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrices
    %%
    for kk=1:mm
         Gs{kk}=EE*G{kk}*EE';
    end
    %%
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs)
        if rijs(m)~=rijs(m-1)
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1)
        for tel2=1:(length(indices)-1)
          if tel1~=tel2
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
		 %%
            for kk=1:mm
                 block= Gs{kk}(positionsh,positionsv); 
                 blocksum=rowsum(block);
                 if blocksum==0
                    isman=0;
                 end
            end
            %%         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1
        manifolds=[manifolds;rij]; 
        R{end+1}=EE';
    end
end



else disp('invalid coupling type number')
return
end



%
% iscell loop
%
end

%% part 3 dynamics seperation 
if isempty(manifolds)
    PS=[];
    return
end
if ~iscell(G)
   A=G;
   if couplingtype==1 % Invasive coupling 
   PS(size(manifolds,1))=struct('manifold',[],'R',[],'Ared',[],'Asyn',[]);

    for kk=1:size(manifolds,1)
        manifold=manifolds(kk,:);
        [Ared, Asyn] = dynamic_seperation_invasive(A,manifold,R{kk});
        PS(kk)=struct('manifold',manifold,'R',R{kk},'Ared',Ared,'Asyn',Asyn);
    end
   end
   if couplingtype==2 % Non-invasive coulpling 
   PS(size(manifolds,1))=struct('manifold',[],'R',[],'Lred',[],'Asyn',[]);
    for kk=1:size(manifolds,1)
        manifold=manifolds(kk,:);
        [Lred, Asyn] = dynamic_seperation_non_invasive(A,manifold,R{kk});
        PS(kk)=struct('manifold',manifold,'R',R{kk},'Lred',Lred,'Asyn',Asyn);
    end
   end
else  % G is a cell
   if couplingtype==1 % Invasive coupling 
   PS(size(manifolds,1))=struct('manifold',[],'R',[],'Ared',[],'Asyn',[]);
    for kk=1:size(manifolds,1)
        manifold=manifolds(kk,:);
        Ared={};
        Asyn={};
          for ll=1:size(G,2)
              A=G{1,ll};
              [Ared{1,end+1}, Asyn{1,end+1}] = dynamic_seperation_invasive(A,manifold,R{kk});
          end
        PS(kk).manifold=manifold;
        PS(kk).R=R{kk};
        PS(kk).Ared=Ared;
        PS(kk).Asyn=Asyn;
    end
   end
   if couplingtype==2 % Non-invasive coulpling 
   PS(size(manifolds,1))=struct('manifold',[],'R',[],'Ared',[],'Asyn',[]);
       for kk=1:size(manifolds,1)
        manifold=manifolds(kk,:);
        Lred={};
        Asyn={};
        for ll=1:size(G,2)
              A=G{1,ll};
              [Lred{1,end+1}, Asyn{1,end+1}] = dynamic_seperation_non_invasive(A,manifold,R{kk});
        end
        PS(kk).manifold=manifold;
        PS(kk).R=R{kk};
        PS(kk).Lred=Lred;
        PS(kk).Asyn=Asyn;
       end
   end
   
end


   
    
%% inline functions    
    
function y=rowsum(bl)
    [nn]=size(bl);
    bl2=bl*ones(nn(2),1);
    bl2=abs(bl2-bl2(1)*ones(nn(1),1));
   if max(bl2)==0
       y=1;
   else
       y=0;
   end 
end
    
function y=rowsumnint(bl)
    [nn]=size(bl);
    bl2=bl*ones(nn(2),1);
    bl2=abs(bl2-bl2(1)*ones(nn(1),1));
   if max(bl2)<epsilon
       y=1;
   else
       y=0;
   end 
end

function [A,partitions]=generate_partition_universal(G,F)
% Initialization 
    A=[];
    partitions=[];
% Dimension
    if iscell(G)==1 
    N=length(G{1});
    else N=length(G);
    end

    A=[0];
    for k=2:N
        B=[];
        for l=1:size(A,1)
            index=(F(k)==F(1:k-1));     %  check the system k have the same dynamics with systems 1,2,...,k-1 
            Al=A(l,:)';                          %  the values of system k can be chosen from based on the row sums
            P=unique(Al(index));         %  remove duplicated entries
            P=[P;max(Al)+1];              %  add the situation where the k system is in a new group 
            for m=1:length(P)
                B=[B;A(l,:) P(m)];
            end
         end
         A=B;
    end
    partitions=A;
    A=A(1:end-1,:);
end

function [A,partitions]=generate_partition_invasive(G, F, varargin)
    A=[];
    partitions=[];
    % Dimension
    if iscell(G)==1 
    N=length(G{1});
    else N=length(G);
    end
    % G is a matrix
    if ~iscell(G)
    Grow = G*ones(size(G,2),1);
        if norm(mod(G,1))==0
        %disp('integer adjacency matrix')
        [Grows,IA,IC] = unique(Grow);
        else % non-integer row sum
            if nargin~=3
            disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
            return
            else
        tolerance=varargin{1};
        tol=tolerance/max(abs(Grow(:)));
        [Grows,IA,IC] = uniquetol(Grow,tol);
            end
        end
        
        if length(IA)==N
            disp('No partial synchronization manifolds exist')
            return
        end
        Grow=[IC F];
         Skip=[];
        for ii=1:N
            if any(ii==Skip)
                continue;
            end
            index=[];
            for jj=1:size(Grow,2)
                temp=(Grow(ii,jj)==Grow(:,jj));
                index=[index temp];
            end
            index_com=prod(index,2); % find common indicator of constant row sum of all matrices 
            Skip=[Skip;find(index_com)]; % store which nodes can be skipped
            IC=[IC index_com];
        end       
        IC=IC*[1:size(IC,2)]';      
            if size(IC,2)==N
            disp('No partial synchronization manifolds exist')
            return
            end
        A=[0];
        for k=2:N
            B=[];
            for l=1:size(A,1)
                index=(IC(k)==IC(1:k-1));  %  check the system k have the same row sum with systems 1,2,...,k-1 
                Al=A(l,:)';                          %  the values of system k can be chosen from based on the row sums
                P=unique(Al(index));         %  remove duplicated entries
                P=[P;max(Al)+1];              %  add the situation where the k system is in a new group 
                for m=1:length(P)
                    B=[B;A(l,:) P(m)];
                end
            end
            A=B;
        end
        partitions=A;
        A=A(1:end-1,:);
  % G is a cell 
    else 
        r=size(G,2);  % G={G1, G2, ... Gr} 
        Grow=[]; % matrix to store sum of row in G1, G2, ..., Gr
        for ii=1:r
            if norm(mod(G{1,ii},1))==0
                Gp=G{1,ii}*ones(N,1);
            else 
                if nargin~=3
                    disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
                    return
                else
                    Gp=G{1,ii}*ones(N,1);
                    tolerance=varargin{1};
                    tol=tolerance/max(abs(Gp(:)));
                    [Gu,IA,Gp] = uniquetol(Gp,tol) % Gp indicates which sums of row are equal w.r.t. the given tolerance
                    if length(IA)==N
                        disp('No partial synchronization manifolds exist')
                        return
                    end
                end
            end
            Grow=[Grow Gp];
        end
        Grow=[Grow F];  % add dynamics indicator
        IC=[];
        Skip=[];
        for ii=1:N
            if any(ii==Skip)
                continue;
            end
            index=[];
            for jj=1:(r+1)
                temp=(Grow(ii,jj)==Grow(:,jj));
                index=[index temp];
            end
            index_com=prod(index,2); % find common indicator of constant row sum of all matrices 
            Skip=[Skip;find(index_com)]; % store which nodes can be skipped
            IC=[IC index_com];
        end
        
            IC=IC*[1:size(IC,2)]';
            if size(IC,2)==N
            disp('No partial synchronization manifolds exist')
            return
            end

        
        A=[0];
        for k=2:N
            B=[];
            for l=1:size(A,1)
                index=(IC(k)==IC(1:k-1));
                Al=A(l,:)';
                P=unique(Al(index));
                P=[P;max(Al)+1];
                for m=1:length(P)
                    B=[B;A(l,:) P(m)];
                end
            end
            A=B;
        end
        partitions=A;
        A=A(1:end-1,:);
    end 
end 
 function [Ared,Asyn]=dynamic_seperation_invasive(A,manifold,R)
    % Initialization 
    Ared=[];
    Asyn=[];
    % A: adjacency matrix 
    % manifold: one manifold presented as a row vector (manifolds{i,:})
    % R: reordering matrix associated with one manifolds, R{i}
        val=unique(manifold);
        kappa=sum(bsxfun(@eq,val,manifold(:))); % number of nodes in each cluster
        K=length(unique(manifold));  % number of clusters 
        Are=[];
        Are=R'*A*R;
        T1=zeros(N-K,N);
        T2=zeros(N-K,N);
        indexr=1;
        indexc=1;
        Rs=[];
        indices=[];
        for mm=1:length(kappa)
            indices=[indices indexc];
            T1([indexr:indexr+kappa(mm)-2],[indexc+1:indexc+kappa(mm)-1])=eye(kappa(mm)-1);
            T2([indexr:indexr+kappa(mm)-2],indexc)=ones(kappa(mm)-1,1);
            indexr=indexr+kappa(mm)-1;
            indexc=indexc+kappa(mm);
        end
        Rij=[];
        for ii=1:K
            for jj=1:K
                Rij(ii,jj)=sum(Are(indices(ii),[indices(jj):indices(jj)+kappa(jj)-1]));
            end
        end
        Rs=sum(Rij,2);
        Ared=T1*Are*T1'-T2*Are*T1';
        Asyn=Rij;        
    end
 function [Lred,Asyn]=dynamic_seperation_non_invasive(A,manifold,R)
     % Initialization 
     Lred=[];
     Asyn=[];
        % A: adjacency matrix 
        % manifold: one manifold presented as a row vector (manifolds{i,:})
        % R: reordering matrix associated with one manifolds, R{i}
        val=unique(manifold);
        kappa=sum(bsxfun(@eq,val,manifold(:))); % number of nodes in each cluster
        K=length(unique(manifold));  % number of clusters 
        Are=[];
        Are=R'*A*R;
        T1=zeros(N-K,N);
        T2=zeros(N-K,N);
        indexr=1;
        indexc=1;
        indices=[];
        for mm=1:length(kappa)
            indices=[indices indexc];
            T1([indexr:indexr+kappa(mm)-2],[indexc+1:indexc+kappa(mm)-1])=eye(kappa(mm)-1);
            T2([indexr:indexr+kappa(mm)-2],indexc)=ones(kappa(mm)-1,1);
            indexr=indexr+kappa(mm)-1;
            indexc=indexc+kappa(mm);
        end
        Rij=[];
        for ii=1:K
            for jj=1:K
                Rij(ii,jj)=sum(Are(indices(ii),[indices(jj):indices(jj)+kappa(jj)-1]));
                if ii==jj
                    Rij(ii,jj)=0;
                end
            end
        end
        Ared=T1*Are*T1'-T2*Are*T1';
        Rs=sum(Are,2);
        Rs(indices)=[];
        Dred=diag(Rs);
        Lred=Dred-Ared;
        Asyn=Rij;        
 end
% end main function
end
