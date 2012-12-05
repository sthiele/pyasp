% redge/3
% pair/2
% pmax(U)


start(p(G,S),S) :-  redge(p(G,S),Y,W), pair(S,Z).

fpath(X,0,S) :- start(X,S).
fpath(p(GY,RY),L+W,S) :- redge(p(GX,RX),p(GY,RY),W), fpath(p(GX,RX),L,S),
                       pmax(N), L+W <= N.

vpair(S,Z)     :- fpath(p(G,Z),L,S), pair(S,Z).

#hide.
#show vpair/2.
