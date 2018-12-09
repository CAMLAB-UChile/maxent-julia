using LinearAlgebra

function phidermatrices(dim,w,wder,s,xc,phi)

  H=zeros(dim)
  A=zeros(dim)
  MC=zeros(1,dim)
  # compute the hessian matrix
  if (dim == 1)
      H[1,1]=dot(phi,xc[1:s,1].^2)
  elseif (dim == 2)
    H=zeros(dim,dim)
    H[1,1]=sum(phi.*(xc[1:s,1].^2))
    H[1,2]=sum(phi.*(xc[1:s,1].*xc[1:s,2]))
    H[2,1]=H[1,2]
    H[2,2]=sum(phi.*(xc[1:s,2].^2))
  elseif (dim == 3)
    H=zeros(dim,dim)
    H[1,1]=dot(phi,xc[1:s,1].^2)
    H[1,2]=dot(phi,xc[1:s,1].*xc[1:s,2])
    H[1,3]=dot(phi,xc[1:s,1].*xc[1:s,3])
    H[2,1]=H[1,2]
    H[2,2]=dot(phi,xc[1:s,2].^2)
    H[2,3]=dot(phi,xc[1:s,2].*xc[1:s,3])
    H[3,1]=H[1,3]
    H[3,2]=H[2,3]
    H[3,3]=dot(phi,xc[1:s,3].^2)
  else
    error("Fatal error! Dimension not coded yet.")
  end

  # compute MA matrix
  MA=zeros(s,dim)
  @inbounds for i=1:s
      MA[i,1:dim]=1/w[i]*wder[i,1:dim]
  end

  # compute A matrix
  if (dim == 1)
    A[1,1]=dot(phi.*MA[1:s,1],xc[1:s,1])
  elseif (dim == 2)
    A=zeros(dim,dim)
    A[1,1]=dot(phi.*MA[1:s,1],xc[1:s,1])
    A[1,2]=dot(phi.*MA[1:s,2],xc[1:s,1])
    A[2,1]=dot(phi.*MA[1:s,1],xc[1:s,2])
    A[2,2]=dot(phi.*MA[1:s,2],xc[1:s,2])
  elseif (dim == 3)
    A=zeros(dim,dim)
    A[1,1]=dot(phi.*MA[1:s,1],xc[1:s,1])
    A[1,2]=dot(phi.*MA[1:s,2],xc[1:s,1])
    A[1,3]=dot(phi.*MA[1:s,3],xc[1:s,1])
    A[2,1]=dot(phi.*MA[1:s,1],xc[1:s,2])
    A[2,2]=dot(phi.*MA[1:s,2],xc[1:s,2])
    A[2,3]=dot(phi.*MA[1:s,3],xc[1:s,2])
    A[3,1]=dot(phi.*MA[1:s,1],xc[1:s,3])
    A[3,2]=dot(phi.*MA[1:s,2],xc[1:s,3])
    A[3,3]=dot(phi.*MA[1:s,3],xc[1:s,3])
  else
    error("Fatal error! Dimension not coded yet.")
  end

  # compute MC matrix
  if (dim == 1)
    MC[1,1]=dot(phi,MA[1:s,1])
  elseif (dim == 2)
    MC[1,1]=dot(phi,MA[1:s,1])
    MC[1,2]=dot(phi,MA[1:s,2])
  elseif (dim == 3)
    MC[1,1]=dot(phi,MA[1:s,1])
    MC[1,2]=dot(phi,MA[1:s,2])
    MC[1,3]=dot(phi,MA[1:s,3])
  else
    error("Fatal error! Dimension not coded yet.")
  end
  return H,MA,A,MC

end
