function fermi_factor=fermi_factor1(omega,mesh,beta,lambda_0)

 lambda_0=0;
if omega < 0
      
     fermi_factor=(exp(beta*(omega+lambda_0))+1)./(exp(beta*(omega+lambda_0))+exp(beta*(-mesh))+exp(beta*(mesh+omega+lambda_0))+1);
%     %lambda_0=lambda_0;
%     %fermi_factor=(exp(beta*(omega-lambda_0))+1)./(exp(beta*lambda_0)+exp(beta*(mesh+omega))+exp(-beta*mesh)+exp(beta*(omega-lambda_0)));
% %     
%     omega=-omega;
%     
%     fermi_factor_1=(exp(beta*(-omega+lambda_0))+1)./(exp(beta*(-omega+lambda_0))+exp(beta*(mesh-lambda_0))+exp(beta*(-mesh-omega+lambda_0))+1);
%     
%     fermi_factor=zeros(1,length(fermi_factor_1));
%     
%     for i=1:length(fermi_factor)
%         
%         fermi_factor(i)=fermi_factor_1(length(fermi_factor)+1-i);
%     end
else

%fermi_factor=(exp(beta*(-omega+lambda_0))+1)./(exp(beta*(-omega+2*lambda_0))+exp(beta*(mesh+lambda_0))+exp(beta*(-mesh-omega+lambda_0))+1);
fermi_factor=(exp(beta*(-omega-lambda_0))+1)./(exp(beta*(-omega-lambda_0))+exp(beta*(mesh))+exp(beta*(-mesh-omega-lambda_0))+1);


end
end