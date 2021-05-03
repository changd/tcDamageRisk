function dist = pos2distVec(lag1,lon1,lag2,lon2,method)
    %% post2dist function but for lat/long vectors
    
    N=length(lag2);
    dist=zeros(size(lag2));
    for i=1:numel(lag2)
        dist(i)=pos2dist(lag1,lon1,lag2(i),lon2(i),method);
    end