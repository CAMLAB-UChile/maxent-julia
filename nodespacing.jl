using LinearAlgebra

function nodespacing(dim,n,ncoord)
  """
 Purpose
 =======
 Compute characteristic nodal spacing to use in basis function support-width

 Input variables
 ===============
 dim = dimension of the domain of discretization
 n = total number of nodes in the domain
 ncoord = vector of nodal coordinates
          1d: ncoord = [x1;x2;...;xn]
          2d: ncoord = [x1 y1;x2 y2;...;xn yn]
          3d: ncoord = [x1 y1 z1;x2 y2 z2;...;x3 y3 z3]

 Return
 ======
 h_node = [h1; h2;...; hn]
"""
  h_node=zeros(n,1)
  distances=zeros(n,n)
  k=1
  @inbounds for i=1:n
    distances[i,:]=sqrt.(sum((ones(n,1)*ncoord[i,:]'-ncoord).^2, dims=2))'
    distances[i,i]=-1
    distances[i,:]=sort(distances[i,:])
    if dim==1
      if n<3
        h_node[i,1]=ncoord[2,1]-ncoord[1,1]
        k=k+1
      else
        h_node[i,1]=distances[i,2]
        k=k+1
      end
    elseif dim==2
      h_node[i,1]=distances[i,3]
      k=k+1
    elseif dim==3
      h_node[i,1]=distances[i,4]
      k=k+1
    else
      error("Fatal error! Dimension not yet coded.")
    end
  end
  return h_node

end
