function Im_d=Im_df_tilda(omega,mesh,hyb,spec_up,spec_down,beta,lambda_0)

% Im_d=fermi_factor1(omega,mesh-sign(omega)*omega,beta,lambda_0).*hyb.*(spec_up+spec_down)/(pi);
Im_d=fermi_factor1(omega,mesh-omega,beta,lambda_0).*hyb.*(spec_up+spec_down)/(pi);

end