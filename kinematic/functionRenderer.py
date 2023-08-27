import mako.template

mytemplate = mako.template.Template(filename='py_functions_template.py')
varList=['th','qv_curr', 'qc_curr', 'qr_curr', 'qi_curr', 'qs_curr', 'qg_curr',\
    'qni_curr', 'qns_curr', 'qnr_curr', 'qng_curr', 'rho', 'pi_phy', 'p', 'dt','dz8w']
dimens=['nz','ny','nx']
renderedF=mytemplate.render(\
    varList=varList,dimens=dimens)

with open('wrf2python_interface.py','w') as f:
    f.write(renderedF)

mytemplate_f90 = mako.template.Template(filename='call_emulator_template.py')

with open('call_emulator.f90','w') as f:
    f.write(mytemplate_f90.render(varList=varList,dimens=dimens))