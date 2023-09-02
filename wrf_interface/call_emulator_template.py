<%
argnames=''
c_argnames=''
for i,var in enumerate(varList):
    if var!='dt':
        argnames+=var+'(its:ite,kts:kte,jts:jte),'
    else:
        argnames+=var+','
    if i%3==0:
        argnames+='&\n    '
    #c_argnames+='float *'+var+','
argnames+='nx_py,ny_py,nz_py,ireturn_py'
c_argnames+='int *nx,int *ny,int *nz'                
%>

nx_py=ite-its+1
ny_py=jte-jts+1
nz_py=kte-kts+1
call emulator_interface(${argnames})
