A1=[0 0 1 0 0 1 0 0;
    0 0 0 0 0 0 0 1;
    1 1 0 0 0 0 1 0;
    1 1 1 0 0 0 0 0;
    1 1 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0; 
    1 1 1 0 0 0 0 0;
    0 0 0 0 1 0 0 0];
Coupling = 1;

G=A1;

N=size(A1,1);
%%%%%%%%%%%%%%%%%%

Grow = G*ones(size(G,2),1);
        if norm(mod(G,1))==0
        disp('integer adjacency matrix')
        [Grows,IA,IC] = unique(Grow);
        else % non-integer row sum
            if nargin~=3
            disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
            return
            else
        epsilon=varargin{1};
        tolu=epsilon/max(abs(Grow(:)));
        [Grows,IA,IC] = unique(Grow,tolu);
            end
        end
        A=[0];
        for k=2:N
            B=[];
            for l=1:size(A,1)
                index=(IC(k)==IC(1:k-1))
                Arow=A(l,:);
                Arow=Arow';
                Poss=unique(Arow(index));
                Poss=[Poss;max(Arow)+1];
                for m=1:length(Poss)
                    B=[B;A(l,:) Poss(m)];
                end
            end
            A=B;
        end
        partitions=A;