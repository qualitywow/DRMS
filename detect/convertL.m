function l = convertL(Nzc, Bcs, zcIndex, Ncs, simPara)
    l = (zcIndex - 1) * length(Bcs) + find(Bcs == Ncs);
    % fprintf(fid, "l = %d\n", l);
    % go through |U|* |Bcs|, `l` still has a chance to be zero
    if (l > simPara.numActive) l = 0; end;
end