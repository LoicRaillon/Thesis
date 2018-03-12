function Y = fct(tfm,par,ts,X,T)

Np                          = length(X);
Y                           = zeros(1,Np);
switch tfm
    case 'rcss' 
        switch par
            case 'a11'
                Rci         = X(:,1);
                Rv          = X(:,2);
                Rw          = X(:,3);
                Ca          = X(:,4);
                Y(1,:)      = -(2.*Rci+2.*Rv+Rw)./(Ca.*Rv.*(2.*Rci+Rw));
            case 'a12' 
                Rci         = X(:,1);
                Rw          = X(:,2);
                Ca          = X(:,3);
                Y(1,:)      = 2./(Ca.*(2.*Rci+Rw));
            case 'a21'
                Rci         = X(:,1);
                Rw          = X(:,2);
                Cw          = X(:,3);
                Y(1,:)      = 2./(Cw.*(2.*Rci+Rw));
            case 'a22'
                Rco         = X(:,1);
                Rci         = X(:,2);
                Rw          = X(:,3);
                Cw          = X(:,4);
                Y(1,:)      = -4.*(Rci+Rco+Rw)./(Cw.*(2.*Rco+Rw).*(2.*Rci+Rw));
            case 'b11'
                Rv          = X(:,1);
                Ca          = X(:,2);
                Y(1,:)      = 1./(Rv.*Ca);
            case 'b13'
                Ca          = X(:,1);
                Y(1,:)      = 1./Ca;
            case 'b21'
                Rco         = X(:,1);
                Rw          = X(:,2);
                Cw          = X(:,3);
                Y(1,:)      = 2./(Cw.*(2.*Rco+Rw));
            case 'b22'
                Rco         = X(:,1);
                Rw          = X(:,2);
                Cw          = X(:,3);
                Y(1,:)      = (2.*Rco)./(Cw.*(2.*Rco+Rw));
        end
    case 'sstf'
        a11d                = 1+X(:,1).*ts;
        a12d                = X(:,2).*ts;
        a21d                = X(:,3).*ts;
        a22d                = 1+X(:,4).*ts;
        b11d                = X(:,5).*ts;
        b13d                = X(:,6).*ts;
        b21d                = X(:,7).*ts;
        b22d                = X(:,8).*ts;
        
        for jj = 1:Np
            Ad              = [a11d(jj) a12d(jj) ; a21d(jj) a22d(jj)];
            Bd              = [b11d(jj) 0 b13d(jj) ; b21d(jj) b22d(jj) 0];
            T               = [1 0 ; -a11d(jj)/a12d(jj) 1/a12d(jj)];
            Acd             = T\Ad*T;
            Beta            = T\Bd;
            Bcd             = [1 0 ; -Acd(2,2) 1]*Beta;
            switch par
                case 'n11'
                    Y(1,jj) = Bcd(1,1);
                case 'n12'
                    Y(1,jj) = Bcd(2,1);
                case 'n22'
                    Y(1,jj) = Bcd(2,2);
                case 'n31'
                    Y(1,jj) = Bcd(1,3);
                case 'n32'
                    Y(1,jj) = Bcd(2,3);
                case 'd1'
                    Y(1,jj) = -Acd(2,2);
                case 'd2'
                    Y(1,jj) = -Acd(2,1);
            end
        end
    case 'rctf'
        Rco                 = X(:,1);
        Rci                 = X(:,2);
        Rv                  = X(:,3);
        Rw                  = X(:,4);
        Cw                  = X(:,5);
        Ca                  = X(:,6);
        a11c                = -(2.*Rci+2.*Rv+Rw)./(Ca.*Rv.*(2.*Rci+Rw));
        a12c                = 2./(Ca.*(2.*Rci+Rw));
        a21c                = 2./(Cw.*(2.*Rci+Rw));
        a22c                = -4.*(Rci+Rco+Rw)./(Cw.*(2.*Rco+Rw).*(2.*Rci+Rw));
        b11c                = 1./(Rv.*Ca);
        b13c                = 1./Ca;
        b21c                = 2./(Cw.*(2.*Rco+Rw));
        b22c                = (2.*Rco)./(Cw.*(2.*Rco+Rw));
        a11d                = 1+a11c.*ts;
        a12d                = a12c.*ts;
        a21d                = a21c.*ts;
        a22d                = 1+a22c.*ts;
        b11d                = b11c.*ts;
        b13d                = b13c.*ts;
        b21d                = b21c.*ts;
        b22d                = b22c.*ts;
        for jj = 1:Np
            Ad              = [a11d(jj) a12d(jj) ; a21d(jj) a22d(jj)];
            Bd              = [b11d(jj) 0 b13d(jj) ; b21d(jj) b22d(jj) 0];
            T               = [1 0 ; -a11d(jj)/a12d(jj) 1/a12d(jj)];
            Acd             = T\Ad*T;
            Beta            = T\Bd;
            Bcd             = [1 0 ; -Acd(2,2) 1]*Beta;
            switch par
                case 'n11'
                    Y(1,jj) = Bcd(1,1);
                case 'n12'
                    Y(1,jj) = Bcd(2,1);
                case 'n22'
                    Y(1,jj) = Bcd(2,2);
                case 'n31'
                    Y(1,jj) = Bcd(1,3);
                case 'n32'
                    Y(1,jj) = Bcd(2,3);
                case 'd1'
                    Y(1,jj) = -Acd(2,2);
                case 'd2'
                    Y(1,jj) = -Acd(2,1);
            end
        end
    case 'tfss'
        n11                 = X(:,1);
        n12                 = X(:,2);
        n22                 = X(:,3);
        n31                 = X(:,4);
        n32                 = X(:,5);
        d1                  = X(:,6);
        d2                  = X(:,7);
        B1                  = [n11 zeros(length(n11),1) n31];
        B2                  = [n12 n22 n32];
        for jj = 1:Np
            Ad              = [0 1 ; -d2(jj) -d1(jj)];
            invZ            = [1 0 ; -d1(jj) 1];
            Bd              = invZ*[B1(jj,:) ; B2(jj,:)];
            AdT             = T*Ad/T;
            BdT             = T*Bd;
            AcT             = (AdT-eye(2))./ts;
            BcT             = BdT./ts;
            switch par
                case 'a21'
                    Y(1,jj) = AcT(2,1);
                case 'a22'
                    Y(1,jj) = AcT(2,2);
                case 'b11'
                    Y(1,jj) = BcT(1,1);
                case 'b13'
                    Y(1,jj) = BcT(1,3);
                case 'b21'
                    Y(1,jj) = BcT(2,1);
                case 'b22'
                    Y(1,jj) = BcT(2,2);
            end
        end
    case 'ssrc'
        a11                 = X(:,1);
%         a12                 = X(:,2);
        a21                 = X(:,3);
        a22                 = X(:,4);
        b11                 = X(:,5);
        b13                 = X(:,6);
%         b21                 = X(:,7);
        b22                 = X(:,8);
        switch par
            case 'Rco'
                Y(1,:)      = -b22./(a22+a21);
            case 'Rci'
                Y(1,:)      = -((b11+a11).*b22+(a22+2.*a21).*b13)./((a22+a21).*b11+a11.*a22+a11.*a21);
            case 'Rv'
                Y(1,:)      = b13./b11;
            case 'Rw'
                Y(1,:)      = ((2.*b11+2.*a11).*b22+2.*a21.*b13)./((a22+a21).*b11+a11.*a22+a11.*a21); % without eq1 and eq4
            case 'Cw'
                Y(1,:)      = -(b11+a11)./(a21.*b13);
            case 'Ca'
                Y(1,:)      = 1./b13;
        end
    case 'tfrc'
        n11                 = X(:,1);
        n12                 = X(:,2);
        n22                 = X(:,3);
        n31                 = X(:,4);
        n32                 = X(:,5);
        d1                  = X(:,6);
        d2                  = X(:,7);
        B1                  = [n11 zeros(length(n11),1) n31];
        B2                  = [n12 n22 n32];
        for jj = 1:Np
            Ad              = [0 1 ; -d2(jj) -d1(jj)];
            invZ            = [1 0 ; -d1(jj) 1];
            Bd              = invZ*[B1(jj,:) ; B2(jj,:)];
            AdT             = T*Ad/T;
            BdT             = T*Bd;
            AcT             = (AdT-eye(2))./ts;
            BcT             = BdT./ts;
            a11             = AcT(1,1);
%             a12             = AcT(1,2);
            a21             = AcT(2,1);
            a22             = AcT(2,2);
            b11             = BcT(1,1);
            b13             = BcT(1,3);
%             b21             = BcT(2,1);
            b22             = BcT(2,2);
            switch par
                case 'Rco'
                    Y(1,jj) = -b22/(a22+a21);
                case 'Rci'
                    Y(1,jj) = -((b11+a11)*b22+(a22+2*a21)*b13)/((a22+a21)*b11+a11*a22+a11*a21);
                case 'Rv'
                    Y(1,jj) = b13/b11;
                case 'Rw'
                    Y(1,jj) = ((2*b11+2*a11)*b22+2*a21*b13)/((a22+a21)*b11+a11*a22+a11*a21);
                case 'Cw'
                    Y(1,jj) = -(b11+a11)/(a21*b13);
                case 'Ca'
                    Y(1,jj) = 1/b13;
            end
        end
end