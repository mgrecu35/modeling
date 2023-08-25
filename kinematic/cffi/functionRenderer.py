import mako.template

mytemplate = mako.template.Template(filename='py_functions_template.py')

renderedF=mytemplate.render(\
    varList=['th','qv_curr', 'qc_curr', 'qr_curr', 'qi_curr', 'qs_curr', 'qg_curr',\
    'qni_curr', 'qns_curr', 'qnr_curr', 'qng_curr', 'rho', 'pi_phy', 'dt','dz8w'],\
    dimens=['nz','ny','nx'])

with open('wrf2python_interface.py','w') as f:
    f.write(renderedF)