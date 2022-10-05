function spec_d=spec_df(omega,lambda_0,Re_d,Im_d,Im_d_tilda,E_0,U)

spec_d=Im_d_tilda./(((omega+lambda_0-2*E_0-U-Re_d).*(omega+lambda_0-2*E_0-U-Re_d)+Im_d.*Im_d));

end