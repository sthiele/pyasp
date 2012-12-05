% gene(G)
% redge(R1,R2).
% cat(G,R).
% pmax(U)

% this logic program yields only one answerset
% it converts into the logic format used by ksip.lp
% redge/3
% vertex/1




% genes that catalyse reactions
cat_gene(X) :- gene(X), cat(X,_).


last(L) :- gene(L), not gene(L+1).

% setup redge format
% count distance clockwise
redge(p(G1,R1),p(G2,R2),W1) :- redge(R1,R2), cat(G1,R1), cat(G2,R2), G2 >G1,
                               last(L), pmax(U), W1 = G2-G1, W2 = (L-G2)+G1, W1<=W2, W1 <=U.

redge(p(G1,R1),p(G2,R2),W1) :- redge(R1,R2), cat(G1,R1), cat(G2,R2), G1 >G2, 
                               last(L), pmax(U), W1 = G1-G2, W2 = (L-G1)+G2, W1<=W2, W1 <=U.

% count distance counter clockwise
redge(p(G1,R1),p(G2,R2),W2) :- redge(R1,R2), cat(G1,R1), cat(G2,R2), G2 >G1, 
                               last(L), pmax(U), W1 = G2-G1, W2 = (L-G2)+G1, W2<W1, W2 <=U.
redge(p(G1,R1),p(G2,R2),W2) :- redge(R1,R2), cat(G1,R1), cat(G2,R2), G1 >G2, 
                               last(L), pmax(U), W1 = G1-G2, W2 = (L-G1)+G2, W2<W1, W2 <=U.

redge(p(G,R1),p(G,R2),0) :- redge(R1,R2), cat(G,R1), cat(G,R2).

#hide.
#show redge(X,Y,W).

#show pmax/1.
