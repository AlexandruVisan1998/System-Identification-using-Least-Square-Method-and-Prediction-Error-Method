% VISAN ALEXANDRU 342 B2

function [Mid] = bj_e(Did,si)
% adaptat dupa armax_e - folosit in ISLAB_6C_bj

%indici structurali
    na = si(1);
    nb = si(2);
    nc = si(3);
    nd = si(4);
    nf = si(5);
%     nk = si(6); % intarziere

    if(nargin < 1 || isempty(Did))
        F = [1 1.5 0.7];
        B = [0 1 0.5];
        C = [1 -1 0.2];
        D = [1 -1.5 0.7];
        P = idpoly(conv(D, F),conv([zeros(1, 1) B], D),conv(C, F),1,1,1);
        Did = gen_data(P,250,1,1);
    end

    n_alfa = max([nf,nb,nc,nd])*2; 
    n_beta = n_alfa;

    Mid = armax_e(Did, [nf+nd, nb+nd, nc+nf]);
%   Mid = arx(Did, [n_alfa, n_beta, 1]);
    
    
    %zgomotul estimat
    e = pe(Mid, Did);
    e = e.y;
    
%   N = 250
    N = length(e);
    suma = nd+nf+nb+nc+nd+nf;
    phi = zeros(N,suma);
    
    % calcul phi - input, output si zgomot 
    for i = 1:N
        
        for j = 1:(nd+nf)
            if (i-j)>0
                phi(i,j) = -Did.y(i-j);
            else
                phi(i,j) = 0;
            end
        end
        
        for j = 1:(nb+nd)
            if (i-j)>0
                phi(i, nd+nf+j) = Did.u(i-j);
            else
                phi(i, nd+nf+j) = 0;
            end
        end
        
        for j = 1:(nc+nf)
            if (i-j)>0
                phi(i,nd+nf+nb+nd+j) = e(i-j);
            else
                phi(i,nd+nf+nb+nd+j) = 0;
            end
        end
        
    end

    theta = phi\Did.y;
    B = theta(1:nb);
    C = theta(1+nb:nb+nc);
    D = theta(1+nb+nc:nb+nc+nd);
    F = theta(1+nb+nc+nd:nb+nc+nd+nf);

    Mid.a = [1];
    Mid.b = [0;B]';
    Mid.c = [1;C]';
    Mid.d = [1;D]';
    Mid.f = [1;F]';
%   Mid = idpoly(1,[zeros(1,nk) B],C,D,F,Ts);
%   Mid = bj(Did,[nb nc nd nf nk]) 
    
    e_final = pe(Mid,Did);
    Mid.NoiseVariance = norm(e_final.y)/sqrt(N-4);

end