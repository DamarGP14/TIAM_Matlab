function Im_s_down=Im_s_downf_tilda(omega,mesh,hyb,spec_b,spec_d,beta,lambda_0)

% Im_s_down=fermi_factor1(omega,mesh-sign(omega)*omega,beta,lambda_0).*hyb.*(spec_b+spec_d)/(pi);
Im_s_down=fermi_factor1(omega,mesh-omega,beta,lambda_0).*hyb.*(spec_b+spec_d)/(pi);

end