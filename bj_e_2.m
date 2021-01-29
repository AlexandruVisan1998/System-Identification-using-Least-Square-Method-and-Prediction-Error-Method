% VISAN ALEXANDRU 342 B2

function [Mid] = bj_e_2(D,si)
% implementare bazata pe radacinile comune ale polinoamelor
% si pe functia ce converteste radacinile in polinoame - poly
%indici structurali
    na = si(1);
    nb = si(2);
    nc = si(3);
    nd = si(4);
    nf = si(5);
%     nk = si(6); % intarziere
    Ts=1;

    D2 = armax_e(D,[nf+nd, nb+nd, nc+nf]);
    aa = D2.A;
    bb = D2.B;
    cc = D2.C;
    
    % Roots
    roots_A = roots(aa);
    roots_B = roots(bb);
    roots_C = roots(cc);
    roots_D = [];
    k=1;

    for i = 1:length(roots_A)
        
        for j = 1:length(roots_B)
        if(roots_A(i) == roots_B(j)) 
            roots_D(k) = roots_A(i);
            k=k+1;
        end
        
    end
    D = poly(roots_D);
    B = bb/D;
    F = aa/D;
    C = cc/F;

    Mid.b = [0;B]';
    Mid.c = [1;C]';
    Mid.d = [1;D]';
    Mid.f = [1;F]';
end