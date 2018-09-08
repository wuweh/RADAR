function y = qr_FUNC(A,COV)
    [~,Q]=qr([A sqrt(COV)]',0);
    y = Q;
end
