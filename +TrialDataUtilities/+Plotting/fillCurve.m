function h = fillCurve(x, y, args, patch_args)
    % fills the area under a curve 
    arguments
        x (:, 1)
        y (:, 1)
        args.baseline (1, 1) = 0;
        patch_args.?matlab.graphics.primitive.Patch;
    end

    xpatch = [x(1); x; x(end)];
    ypatch = [args.baseline; y; args.baseline];

    patch_args_c = namedargs2cell(patch_args);
    patch(x, y, 'b', patch_args)
end