function fermi_factor=fermi_factor1(omega,mesh,beta)

if omega < 0
      
     fermi_factor=(exp(beta*(omega))+1)./(exp(beta*(omega))+exp(beta*(-mesh))+exp(beta*(mesh+omega))+1);

else
fermi_factor=(exp(beta*(-omega))+1)./(exp(beta*(-omega))+exp(beta*(mesh))+exp(beta*(-mesh-omega))+1);

end
end