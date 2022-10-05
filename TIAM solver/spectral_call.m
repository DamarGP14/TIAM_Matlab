function spectral_call=spectral_call(mesh,lambda_0,D,Re_E,Im_E,Im_E_tilda)

spectral_call=Im_E_tilda./((mesh+lambda_0-D-Re_E).*(mesh+lambda_0-D-Re_E)+Im_E.*Im_E);

end