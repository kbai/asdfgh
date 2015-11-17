function hw7p1()
G1=[1 1 1 1; 1 -3 4 5]';
G2=[1 1 1 1; -0.1 0.3 -0.4 0.5]';
G3=[1 1 1 1; 101 97 104 105]';
G4=[1 1 1 1; 3 3 3 3]';

R1=pseduoInverse(G1)
R2=pseduoInverse(G2)
R3=pseduoInverse(G3)
R4=pseduoInverse(G4)

% least square inverse
L1=LLSInverse(G1)
L2=LLSInverse(G2)
L3=LLSInverse(G3)
L4=LLSInverse(G4)

function r=LLSInverse(G)
r = inv(G'*G)*G';

function r=pseduoInverse(G)
[U,S,V]=svd(G); % check if U*S*V' == G
eps=1E-16;
for i=1:min(size(S))
    if(S(i,i) < eps)
        S(i,i) = 0;
    else
        S(i,i) = 1.0/S(i,i);
    end
end
r=V*S'*U';
