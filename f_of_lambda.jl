function f_of_lambda(dim,w,xc,lambda,len)

  if dim==1
    Z_I=w.*exp.(-(xc*lambda))
    Z=sum(Z_I)
    phi=Z_I/Z # contribuitng basis function evaluated at xc = xi - x
    gam=log(Z)
    temp1=-xc.*phi
    dgam=sum(temp1[:]) # gradient of log(Z)
    hgam = sum( phi.*(xc).^2) - dgam*dgam
  else
    Z_I=w.*exp.(-(xc*lambda))
    Z=sum(Z_I)
    phi=Z_I/Z # basis function
    gam=log(Z)
    dgam=zeros(dim,1)
    hgam=zeros(dim,dim)
    @inbounds for id=1:dim
      a=sum(-(xc[:,id]).*phi[:])
      dgam[id]=sum((-xc[:,id]).*phi[:]) # gradient of log(Z)
    end
    @inbounds for id=1:dim
      @inbounds for jd=1:dim
          hgam[id,jd]=sum( phi[:].*( -xc[:,id] ).*(-xc[:,jd]) )- dgam[id]*dgam[jd]
      end
    end
  end
  return gam,dgam,hgam,phi

end
