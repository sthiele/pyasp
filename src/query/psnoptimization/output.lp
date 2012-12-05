%pre-format easier for parsing.
i(si(U,1,end),a,V) :- clause(subset(U,1,end),1,V).
i(si(U,1,end),i,V) :- clause(subset(U,-1,end),1,V).
i(and(F,Z,S),a,V) :- clause(subset(F,Z,S),L,V), L>1.
i(si(U,1,end),a,and(S)):- in(U,1,S), i(and(S),a,V).
i(si(U,1,end),i,and(S)):- in(U,-1,S), i(and(S),a,V).

#hide.
#show i/3.
#show clause/3.
