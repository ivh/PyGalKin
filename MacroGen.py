"""
 MacroGen.py

 For documentation see the major file PyCigale.py
 
 This file contains the functions that create ADHOC
 macro files (.adm) and distribute them in the data
 directories
"""
from PyGalKin import *

## WHERE WILL ADHOC FIND THE DATA?
pathprefix='C:\\data\\CIGALE-2002\\'



def doreductions():
    """ distribute the reduct.adm files for all objects
    and register each macro into a global macro
    """
    
    globalmacro=open('global_reduct_autogen.adb','w')
    globalcounter=0
    reductname='reduct.adm'
    
    for i in os.listdir('.'):
        if os.path.isdir(i):
            os.chdir(i)    

            if 'p800.adp' in os.listdir('.'): paramname='p800.adp'
            elif 'p2000.adp' in os.listdir('.'): paramname='p2000.adp'
            else: print("There is no parameter file.")

            globalmacro.seek(globalcounter*256*3)
            globalmacro.write(pathprefix + i)
            globalmacro.seek((globalcounter*256*3)+256)
            globalmacro.write(paramname)
            globalmacro.seek((globalcounter*256*3)+512)
            globalmacro.write(reductname)
            globalcounter+=1

            par=C.ADP(paramname)
            workingdir= pathprefix + i + bs
            print(workingdir)
            reduc(reductname,workingdir,par)
            
            os.chdir('..')
        else: pass

    globalmacro.seek(196607)
    globalmacro.write(null)
    globalmacro.close()



def docubes():
    """ generate the cube.adm macros for all objects
    and register them in a global macro
    """

    globalmacro=open('global_cube_autogen.adb','w')
    globalcounter=0
    cubename='cube.adm'

    for i in os.listdir('.'):
        if os.path.isdir(i):
            os.chdir(i)
            cubefile=open(cubename,'w')
            if 'p800.adp' in os.listdir('.'): paramname='p800.adp'
            elif 'p2000.adp' in os.listdir('.'): paramname='p2000.adp'
            else: print("There is no parameter file.")

            globalmacro.seek(globalcounter*256*3)
            globalmacro.write(pathprefix + i)
            globalmacro.seek((globalcounter*256*3)+256)
            globalmacro.write(paramname)
            globalmacro.seek((globalcounter*256*3)+512)
            globalmacro.write(cubename)
            globalcounter+=1
            
            for j in os.listdir('.'):
                if os.path.isdir(j):
                    os.chdir(j)
                    workingdir= pathprefix + i + bs + j + bs
                    print(workingdir)
                    
                    
                    
                    if 'DATA' in os.listdir('.'):   ### DATA DIRECTORIES

                        #numbercycles=createADP('../' + paramname, i, j)
                        print(numbercycles)
                        cube(cubefile,workingdir,j,str(numbercycles))   ### WRITE THE LINE INTO THE CUBE-MACRO
                        cubefile.write(nl)
                        
                        if not os.path.islink('../neb_org.ad3'):
                            os.symlink(j + '/org.ad3','../neb_org.ad3')


                    else:                            ### CALIB DIRECTORIES
                        if not os.path.islink('../eta_org.ad3'):
                            os.symlink(j + '/' + j + '.AD3','../eta_org.ad3')


                    os.chdir('..')
                else: pass
                
            #cubefile.write(nl)
            cubefile.close()
            os.chdir('..')
        else: pass
    globalmacro.seek(196607)
    globalmacro.write(null)
    globalmacro.close()




def reduc(filename,workingdir,par):
    """ making an .ADM file for the data reduction """

    ## NAME OF THE SKY BOXES ADZ FILE
    sky='sky.adz'

    outfile=open(filename,'w')

    # SPECTRAL SMOOTH OF THE CALIBRATION FILE
    outfile.write('SMOOTHZ fin=' + workingdir + 'eta_org.ad3 fout=' + workingdir + 'scratch1.ad3 gauss=2 ;'+nl+nl)

    # PHASE COMPUTATION
    # ALL IMPORTANT FRQUENCES AND PARAMETERS GO IN HERE
    outfile.write('PHASE fin=' + workingdir + 'scratch1.ad3 ph_bru=' + workingdir + 'Phas_bru.AD2 ph_prb=' + workingdir + 'Phas_prb.AD2 cal_sum=' + workingdir + 'Sum_cal.AD2 cal_rv=' + workingdir + 'RV_cal.AD2 xl1='+ str(par.xl1) +' xc='+ str(par.xc) +' yc='+ str(par.yc) +' ellipt=' + str(par.ellipt) + ' center=0 elliptCCD=0 p_pf='+ str(par.p) +' lmbd_p='+ str(par.xlp) +' lmbd_cal='+ str(par.xleta) +' lmbd_scancal='+ str(par.xlbeta) +' lmbd_obj='+ str(par.xlneb) +' lmbd_scanobj='+ str(par.xlbneb) +' minangle=' + str(par.min_a) + ' maxangle=' + str(par.max_a) + ' minradius=' + str(par.min_r) + ' maxradius=' + str(par.max_r) + ' algorithm=1 chg_chnls=2 ;'+nl+nl)

    # SCALE ORIGINAL DATA UP
    outfile.write('ARITHM fin=' + workingdir + 'neb_org.ad3 fout=' + workingdir + 'scratch1.ad3 oper=0 a=10 b=0;'+nl+nl)

    # SPECTRAL SMOOTH OF THE DATA 
    outfile.write('SMOOTHZ fin=' + workingdir + 'scratch1.ad3 fout=' + workingdir + 'scratch2.ad3 gauss=2 ;'+nl+nl)

    # SPATIAL SMOOTH
    outfile.write('SMOOTHXY fin=' + workingdir + 'scratch2.ad3 fout=' + workingdir + 'scratch1.ad3 gauss=2 ;'+nl+nl)

    # APPLY THE PHASE
    outfile.write('LAMBDA fin=' + workingdir + 'scratch1.ad3 ph_prb=' + workingdir + 'Phas_prb.AD2 fout=' + workingdir + 'scratch2.ad3 ;'+nl+nl)

    # SKY SUBTRACTION
    outfile.write('SUBSTRACT_OH fin=' + workingdir + 'scratch2.ad3 fout=' + workingdir + 'l_oh.ad3 fzones=' + workingdir + sky + ' fbmp=' + workingdir + 'OH.BMP ftxt=' + workingdir + 'OH.TXT mode_ref=0 nbpass=1 fact1=1.000000 mode_cont=0 mode_prof=0;'+nl+nl)

    # FLIPPING X AN Y DIMENSION - SEMMS NOT TO BE NECESSARY
    #outfile.write('GEOMETRY oper=1 fin=' + workingdir + 'scratch1.ad3 fout=' + workingdir + 'l_oh.ad3 inversX=1 inversY=1 inversZ=0 ;'+nl+nl)

    # CALCULATE RV, MONO, WIDTH
    outfile.write('MONORV fin=' + workingdir + 'l_oh.ad3 fmono=' + workingdir + 'Mono_ORG.AD2 fcont=' + workingdir + 'Cont_ORG.AD2 frv=' + workingdir + 'RV_ORG.AD2 fwidths=' + workingdir + 'Width_ORG.AD2 mode_emiss=0 algorithm=1 chg_chnls=2 mode_cont=2 ;'+nl+nl)

    # MORE SPECTRAL SMOOTH
    outfile.write('SMOOTHZ fin=' + workingdir + 'l_oh.ad3 fout=' + workingdir + 'scratch1.ad3 gauss=2 ;'+nl+nl)

    # MORE SPATIAL SMOOTH 
    outfile.write('SMOOTHXY fin=' + workingdir + 'scratch1.ad3 fout=' + workingdir + 'l_5.ad3 gauss=5 ;'+nl+nl)

    # CALCULATE RV, MONO, WIDTH WITH SLOPE CHANGE METHOD
    outfile.write('MONORV fin=' + workingdir + 'l_5.ad3 fmono=' + workingdir + 'Mono_5_sl.AD2 fcont=' + workingdir + 'Cont_5_sl.AD2 frv=' + workingdir + 'RV_5_sl.AD2 fwidths=' + workingdir + 'Width_5_sl.AD2 mode_emiss=0 algorithm=1 chg_chnls=2 mode_cont=2 ;'+nl+nl)

    # CALCULATE RV, MONO, WIDTH WITH HISTOGRAM METHOD
    outfile.write('MONORV fin=' + workingdir + 'l_5.ad3 fmono=' + workingdir + 'Mono_5_fl.AD2 fcont=' + workingdir + 'Cont_5_fl.AD2 frv=' + workingdir + 'RV_5_fl.AD2 fwidths=' + workingdir + 'Width_5_fl.AD2 mode_emiss=0 algorithm=2 pc_mono=35 mode_cont=1 pc_cont=10 ;'+nl+nl)

    #MASK
    outfile.write('ARITHM2 fin1=' + workingdir + 'RV_5_sl.AD2 fin2=' + workingdir + 'Mono_5_sl.AD2 fout=' + workingdir + 'RV_5_sl_m.AD2 oper=7 a=4 b=3.1E+38 c=0 ;'+nl+nl)

    #MASK
    outfile.write('ARITHM2 fin1=' + workingdir + 'RV_5_fl.AD2 fin2=' + workingdir + 'Mono_5_fl.AD2 fout=' + workingdir + 'RV_5_fl_m.AD2 oper=7 a=4 b=3.1E+38 c=0 ;'+nl+nl)

    outfile.close()






def cube(outfile,workingdir,rawdir,cycle_end):
    """ generating a macro for cube integration"""
    
    tmpfile='adw_temp_1.tmp'
    if not os.path.exists(tmpfile): os.system('touch '+tmpfile)
    outfile.write('INTEGR_ADA fin=' + workingdir + rawdir + '.ADT' + ' fout=' +  workingdir + 'org.ad3' + ' ftmp=' + workingdir[:-5] + tmpfile + ' cycle_dep=1 cycle_end='+ cycle_end +' norm=1 reman=0 rempix=0 parasit=1 sigma_parasit=5 ;')







