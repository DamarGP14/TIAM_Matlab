function spec_down=spec_downf(omega,lambda_0,Re_down,Im_down,Im_down_tilda,E_0)

spec_down=Im_down_tilda./(((omega+lambda_0-E_0-Re_down).*(omega+lambda_0-E_0-Re_down)+Im_down.*Im_down));

end