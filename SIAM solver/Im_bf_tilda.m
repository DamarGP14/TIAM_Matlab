function Im_b=Im_bf_tilda(omega,mesh,hyb,spec_up,spec_down,beta,lambda_0)

Im_b=fermi_factor1(omega,mesh-omega,beta,lambda_0).*hyb.*(spec_up+spec_down)/(pi);

end