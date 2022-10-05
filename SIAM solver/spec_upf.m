function spec_up=spec_upf(omega,lambda_0,Re_up,Im_up,Im_up_tilda,E_0)

spec_up=Im_up_tilda./(((omega+lambda_0-E_0-Re_up).*(omega+lambda_0-E_0-Re_up)+Im_up.*Im_up));

end