function spec_b=spec_bf(omega,lambda_0,Re_b,Im_b,Im_b_tilda)

spec_b=Im_b_tilda./(((omega+lambda_0-Re_b).*(omega+lambda_0-Re_b)+Im_b.*Im_b));

end