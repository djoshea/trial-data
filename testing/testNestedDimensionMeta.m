sz = struct();
sz.A = 1;
sz.B = 2;
sz.C = 3;
sz.D = 4;
sz.E = 5;
sz.F = 6;
sz.G = 7;

ndm = NestedDimensionMeta({{'A', 'B', 'C'}, {'D', 'E'}, {'F', 'G'}});

data = ndm.emptyFromSizes(sz);
vals = ndm.sizeValuesFromData(ndm.filterByName(data, 'G', 1:4))
assert(vals.G == 4);



