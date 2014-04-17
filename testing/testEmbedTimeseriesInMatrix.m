N = 5;
delta = 0.001;

RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 0));

for i = 1:5
    time = i:delta:2*i;
    timeCell{i} = time + 0.01*randn(size(time));
    dataCell{i} = time;
end
    
[mat, tvec] = embedTimeseriesInMatrix(dataCell, timeCell, 'timeDelta', delta);

%pmat(mat);