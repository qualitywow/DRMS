function val = mod2idx(input, n)
    val = mod(input, n);
    if (val == 0)
        val = n;
    end
end