function neighbors(dim,radius,xc,n)
    node_ind=1:n
    dist=zeros(n,1)
    for i=1:dim
      dist=dist+xc[:,i].^2
    end
    dist=sqrt.(dist)
    contribute=node_ind[dist[:,1].<radius[:,1]]
    len=length(contribute)
    return contribute,len
end
