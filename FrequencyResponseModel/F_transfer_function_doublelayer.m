function Result = F_transfer_function_doublelayer(freq,dl1,dl2,D1,D2,hf,tau)

    dl2 = dl1 + dl2;

    %% complex constants
    img   = 1i;
    quath = 0.5 + 0.5i;
    quat  = 1.0 + 1.0i;

    %% numerical parameters
    nmax = 1000;   % spatial discretization points

    %% preallocation
    C = zeros(nmax+1,1);
    I = zeros(nmax+1,1);
    F = zeros(nmax+1,1);

    Gain  = zeros(length(freq),1);
    Phase = zeros(length(freq),1);

    %% loop over frequency
    for ii = 1:length(freq)

        frequency = freq(ii);
        omega = 2*pi*frequency;

        dk1 = sqrt(2*omega/D1);
        dk2 = sqrt(2*omega/D2);

        %% depth-wise calculation
        for n = 1:nmax+1

            z = (n-1)/nmax * dl1;

            % ===== exponent arguments =====
            A1 = quat * dl1 * dk2;
            A2 = quat * dl2 * dk2;

            B1 = quat * dl1 * dk1;
            B2 = quat * z   * dk1;

            % ===== exponential shifting =====
            M1 = max([real(A1), real(A2)]);
            M2 = max([real(B1), real(B2), 0]);   % 0 corresponds to exp(0)

            EA1 = exp(A1 - M1);
            EA2 = exp(A2 - M1);

            EB1 = exp(B1 - M2);
            EB2 = exp(B2 - M2);
            E1  = exp(-M2);   % shifted exp(0)

            % ===== numerator & denominator =====
            num = exp(-quath*z*dk1) .* ( ...
                    sqrt(D1)*(EA1 + EA2).*(EB1 + EB2) ...
                  - sqrt(D2)*(EA1 - EA2).*(EB1 - EB2) );

            den = sqrt(D1)*(EA1 + EA2).*(EB1 + E1) ...
                - sqrt(D2)*(EA1 - EA2).*(EB1 - E1);

            C(n) = num ./ den;

            I(n) = C(n) ./ (1 + img*tau*omega);
            F(n) = exp(-z*hf);
        end

        %% trapezoidal integration
        TotalIref = 0;
        TotalI    = 0;

        for n = 1:nmax
            TotalIref = TotalIref + 0.5*(F(n)     + F(n+1));
            TotalI    = TotalI    + 0.5*(F(n)*I(n)+ F(n+1)*I(n+1));
        end

        %% frequency response
        H = TotalI / TotalIref;

        Gain(ii)  = 20*log10(abs(H));
        Phase(ii) = atan2(imag(H),real(H)) * 180/pi;

    end

    Result = [Gain, Phase];

end
