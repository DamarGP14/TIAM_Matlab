function n_f=fermi(E,beta,lambda_0)
lambda_0=0;
n_f=1./(exp(beta.*(E+lambda_0))+1);

end