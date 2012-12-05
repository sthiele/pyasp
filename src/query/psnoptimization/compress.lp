
keep(X) :- stimuli(X).
keep(X) :- inhibitor(X).
keep(X) :- readout(X).

%end of input

%marked for reduction
marked(V) :- vertex(V), not keep(V).

%Max number of reductions.
reductions(R) :- R={marked(V)}.

%Only one vertex could be reduced in each reduction.
1{reduction(I,V) : marked(V)}1 :- reductions(R), I=1..R.
%Each vertex could be reduced just once.
:- marked(V), 2{reduction(_,V)}.

%Force an order of reductions related to the vertex id's. (to obtain a unique compression)
% rel(V,R) :- reduction(R,V), vertex(V).
:-reduction(I1,V1), reduction(I2,V2),I1<I2, V1>V2.
% :- rel(V,I,R), I+1<R. %I+1 because I starts form 0
% :- rel(V1,I1,R1), rel(V2,I2,R2), I1>I2, R1<R2, I1+1>=R1, I2+1>=R2.

%Each readuction could be applied or not.
{app(I)} :- reductions(R), I=1..R.
app(1).

%Read the input network
tmp_vertex(1,V) :- vertex(V).
tmp_edge(1,V,W,S) :- edge(V,W,S).

tmp_vertex(I,V) :- keep(V), reductions(R), I=2..R+1.

%From reduction R to reduction R+1, we keep one vertex if it's not the one who we are
%reducing or if the reduction it's not applied.
tmp_vertex(R+1,V) :- tmp_vertex(R,V), not reduction(R,V), reductions(RM),R<=RM.
tmp_vertex(R+1,V) :- tmp_vertex(R,V), not app(R), reductions(RM), R<=RM.

%From reduction R to reduction R+1, we keep an edge if none of it's vertices are the
%one who we are reducing or if the reduction it's not applied.
tmp_edge(R+1,V,W,S) :- tmp_edge(R,V,W,S), reductions(RM), R<=RM, not reduction(R,V),
		       not reduction(R,W).


tmp_edge(R+1,V,W,S) :- tmp_edge(R,V,W,S), not app(R), reductions(RM), R<=RM.

%When a reduction R is applied, we need to update the network by adding new edges in
%order to keep the previous paths.
tmp_edge(R+1,T,W,S1*S2) :- tmp_edge(R,T,V,S1), tmp_edge(R,V,W,S2), reduction(R,V),
app(R), reductions(RM), R<=RM.

%Force to apply reductions when:
%in(V) = 1 and out(V) = 1
 :- 1{tmp_edge(R,V,_,_)}1, 1{tmp_edge(R,_,V,_)}1, reduction(R,V), not app(R).
%in(V) > 1 and out(V) = 1
:- 1{tmp_edge(R,V,_,_)}1, 2{tmp_edge(R,_,V,_)}, reduction(R,V), not app(R).
%in(V) = 1 and out(V) > 1
:- 2{tmp_edge(R,V,_,_)}, 1{tmp_edge(R,_,V,_)}1, reduction(R,V), not app(R).

%Force to NOT apply a reduction when in(V) > 1 and out(V) > 1.
 :- 2{tmp_edge(R,V,_,_)}, 2{tmp_edge(R,_,V,_)}, reduction(R,V), app(R).

%Reduction result
compress_vertex(V) :- tmp_vertex(R+1,V), reductions(R).
compress_edge(V,W,S) :- tmp_edge(R+1,V,W,S), reductions(R).

#hide.
#show compress_vertex/1.
#show compress_edge/3.