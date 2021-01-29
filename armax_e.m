% VISAN ALEXANDRU 342 B2

function [Mid]=armax_e(D, si)
    % armax_e implementeaza MCMMPE pentru sisteme ARMAX
    % si=[na, nb, nc, nk]


    na = si(1);
    nb = si(2);
    nc = si(3);
%     nk = si(4); % intarzierea modelului

    if(nargin < 1 || isempty(D))
        A = [1 -1.5 0.7];
        B = [0 1 0.5];
        C = [1 -1 0.2];
        P = idpoly(A,B,C,1,1,1);
        D = gen_data(P,250,1,1);
    end

    n_alfa = max([na,nb,nc])*2; 
    n_beta = n_alfa;

    Mid = arx(D, [n_alfa, n_beta, 1]);

    %zgomotul estimat
    e = pe(Mid,D) ;% Noise estimation. 
    e = e.y;
    
%     N = 250
    N = length(e);
    suma = na+nb+nc;
    phi = zeros(N,suma);
    
    % calcul phi - input, output si zgomot 
    for i = 1:N
        
        for j = 1:na
            if (i-j>0)
                phi(i, j) = -D.y(i-j);
            else
                phi(i, j) = 0;
            end
        end
        
        for j=1:nb
            if (i-j>0)
                phi(i, na+j)=D.u(i-j);
            else
                phi(i, na+j)=0;
            end
        end
    
        for j=1:nc
            if (i-j>0)
                phi(i, na+nb+j) = e(i-j);
            else
                phi(i, na+nb+j) = 0;
            end
        end
    end
    
    theta = phi\D.y;
    A = theta(1:na);
    B = theta(na+1:na+nb);
    C = theta(na+nb+1:na+nb+nc);

    Mid.a = [1;A]';
    Mid.b = [0;B]';
    Mid.c = [1;C]';
%   Mid = idpoly(A,[zeros(1,nk) B],C,1,1,Ts);

    e_final = pe(Mid,D);
    Mid.NoiseVariance = norm(e_final.y)/sqrt(N-4);
end
