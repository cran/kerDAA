kerdaa = function(X1, Y1, X2, Y2, perm=0) {
  X = rbind(X1, X2)
  Y = rbind(Y1, Y2)
  
  m = nrow(X1)
  n = nrow(X2)
  N = m + n
  
  Dx = dist(X)^2
  sigma = median(Dx)
  if (sigma == 0) {
    sigma = sigma + 0.1
  }
  KX_g = exp(-as.matrix(Dx)/sigma/2)
    
  Dy = dist(Y)^2
  sigma = median(Dy)
  if (sigma == 0) {
    sigma = sigma + 0.1
  }
  KY_g = exp(-as.matrix(Dy)/sigma/2)
  
  KX_l = X%*%t(X)
  KY_l = Y%*%t(Y)
  
  I.m = diag(1,m)
  I.1 = rep(1,m)
  H1 = I.m-I.1%*%t(I.1)/m

  I.n = diag(1,n)
  I.1 = rep(1,n)
  H2 = I.n-I.1%*%t(I.1)/n

  H = rbind( cbind(H1, matrix(0, m, n)), cbind(matrix(0,n,m), H2) )
  
  Kx_g = H%*%KX_g%*%H
  Ky_g = H%*%KY_g%*%H
  
  Kx_l = H%*%KX_l%*%H
  Ky_l = H%*%KY_l%*%H
  
  K_g = Kx_g*Ky_g
  diag(K_g) = 0
  
  K_l = Kx_l*Ky_l
  diag(K_l) = 0
  
  K1_g = sum(K_g[1:m,1:m])/m/(m-1)
  K2_g = sum(K_g[(m+1):N,(m+1):N])/n/(n-1)
  
  K1_l = sum(K_l[1:m,1:m])/m/(m-1)
  K2_l = sum(K_l[(m+1):N,(m+1):N])/n/(n-1)
  
  p1 = m*(m-1)/N/(N-1)
  p2 = p1*(m-2)/(N-2)
  p3 = p2*(m-3)/(N-3)
  
  q1 = n*(n-1)/N/(N-1)
  q2 = q1*(n-2)/(N-2)
  q3 = q2*(n-3)/(N-3)
  
  mu_g = sum(K_g)/N/(N-1)
  
  A_g = sum(K_g^2)
  B_g = sum(rowSums(K_g)^2) - A_g
  C_g = sum(K_g)^2 - 2*A_g - 4*B_g
  
  var_K1_g = (2*A_g*p1 + 4*B_g*p2 + C_g*p3)/m/m/(m-1)/(m-1) - mu_g^2
  var_K2_g = (2*A_g*q1 + 4*B_g*q2 + C_g*q3)/n/n/(n-1)/(n-1) - mu_g^2
  cov_g = C_g/N/(N-1)/(N-2)/(N-3) - mu_g^2
  
  var_g = var_K1_g + var_K2_g - 2*cov_g
  
  stat_g = (K1_g - K2_g)/sqrt(var_g)
  
  mu_l = sum(K_l)/N/(N-1)
  
  A_l = sum(K_l^2)
  B_l = sum(rowSums(K_l)^2) - A_l
  C_l = sum(K_l)^2 - 2*A_l - 4*B_l
  
  var_K1_l = (2*A_l*p1 + 4*B_l*p2 + C_l*p3)/m/m/(m-1)/(m-1) - mu_l^2
  var_K2_l = (2*A_l*q1 + 4*B_l*q2 + C_l*q3)/n/n/(n-1)/(n-1) - mu_l^2
  cov_l = C_l/N/(N-1)/(N-2)/(N-3) - mu_l^2
  
  var_l = var_K1_l + var_K2_l - 2*cov_l
  
  stat_l = (K1_l - K2_l)/sqrt(var_l)
  
  pval_g = 2*pnorm(-abs(stat_g))
  pval_l = 2*pnorm(-abs(stat_l))
  
  x = c(pval_g, pval_l)
  cauchy.t = sum(tan((0.5-pmin(x,0.99))*pi))/length(x)
  pval_o = 1 - pcauchy(cauchy.t)
  
  res = list()
  res$stat_g = stat_g
  res$stat_l = stat_l
  res$pval = pval_o
  
  
  if (perm>0) {
    temp1 = temp2 = rep(0, perm)
    for (i in 1:perm) {
      id = sample(N, replace = FALSE)
      K_g_i = K_g[id,id]
      K_l_i = K_l[id,id]
      
      K1_g_i = sum(K_g_i[1:m,1:m])/m/(m-1)
      K2_g_i = sum(K_g_i[(m+1):N,(m+1):N])/n/(n-1)
      
      K1_l_i = sum(K_l_i[1:m,1:m])/m/(m-1)
      K2_l_i = sum(K_l_i[(m+1):N,(m+1):N])/n/(n-1)
      
      stat_g_i = (K1_g_i - K2_g_i)/sqrt(var_g)
      stat_l_i = (K1_l_i - K2_l_i)/sqrt(var_l)
      
      temp1[i] = stat_g_i
      temp2[i] = stat_l_i
    }
    pval_g_perm = 2*length(which(temp1>=abs(stat_g)))/perm
    pval_l_perm = 2*length(which(temp2>=abs(stat_l)))/perm
  }
  
  x = c(pval_g_perm, pval_l_perm)
  cauchy.t = sum(tan((0.5-pmin(x,0.99))*pi))/length(x)
  pval_o_perm = 1 - pcauchy(cauchy.t)
  
  res$pval_perm = pval_o_perm
  
  
  return(res)
}

















