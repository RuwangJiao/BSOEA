function radius = ReduceRadius(R, k, MaxK)
    %% The shrink of the dynamic niching radius
    cp = 5;
    z        = 1e-8;
    Nearzero = 1e-15;
    B        = MaxK./power(log((R + z)./z), 1.0./cp);
    B(B==0)  = B(B==0) + Nearzero;
    f        = R.* exp( -(k./B).^cp );
    tmp      = find(abs(f-z) < Nearzero);
    f(tmp)   = f(tmp).*0 + z;
    radius     = f - z;
    radius(radius<=0) = 0;
end