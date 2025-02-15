# for generating fatband plot py file

#####################  Custom Variables  #######################
modeofplot = 'o'
spins_mode = "[0,1]" # spin plot mode, "[0,1]" for spin considered DFT, band spin polarization can be showed 
Emin = '-15'
Emax = '30'
cmap = "\'hot_r\'" # the pattern of color map, see more in "https://zhuanlan.zhihu.com/p/114420786",spins can use seismic, orbitals can use hot_r
targetelements = ['Fe', 'Fe', 'Fe', 'Pt', 'Pt', 'Pt']
targetorbitals = ['s', 'p', 'd', 's', 'p', 'd']
fermienergy = '7.1567'

#####################  Custom Variables  #######################



orbitals_name= {'s':"[0]",  'py':"[1]",    'pz':"[2]",     'px':"[3]",    'dxy':"[4]",   'dyz':"[5]",    'dz2':"[6]",    'dxz':"[7]",  'dx2-y2':"[8]", 'p':"[1,2,3]", 'd':"[4,5,6,7,8]", 'f':"[9,10,11,12,13,14,15]"}

#####################  From POSCAR Get Atoms_dict  #######################
names = " "
numbers = " "
with open("POSCAR",mode="r") as poscar:
    for i in range(5) :
        poscar.readline()
    names = poscar.readline()
    numbers = poscar.readline()

names = names.rstrip("\n")
numbers = numbers.rstrip("\n")

names = names.split(" ")
numbers = numbers.split(" ")
#print(len(names))
i = 0
j = 0
atoms_index = []
while True :
    if (i+1) > len(names) :
        break
    tmp = []
    for k in range(int(numbers[i])):
        tmp.append(j)
        j = j + 1
    atoms_index.append(str(tmp))
    i = i + 1

atom_dict = dict(zip(names, atoms_index))
#####################  From POSCAR Get Atoms_dict  #######################


#####################  From KPOINTS Get Knames And Ksticks  #######################



knames = [] # name of each k labels
kticks = [] # positions of k labels

with open('KPOINTS', mode='r') as kpoints :
    kpoints.readline()
    k_interval = int(kpoints.readline())
    kpoints.readline()
    kpoints.readline()
    i = 0
    while True :
        tmp_k = kpoints.readline()
        if not tmp_k :
            knames.append(end_tmp_k[len(end_tmp_k)-1])
            kticks.append(i)
            break
        tmp_k = tmp_k.rstrip('\n')
        tmp_k = tmp_k.rstrip()
        tmp_k = tmp_k.split(' ')
        if i == 0 :
            knames.append(tmp_k[len(tmp_k)-1])
        else :    
            if (tmp_k[len(tmp_k)-1] ==  end_tmp_k[len(end_tmp_k)-1]) :
                knames.append(tmp_k[len(tmp_k)-1])
            else :
                knames.append( end_tmp_k[len(end_tmp_k)-1]+ "|" + tmp_k[len(tmp_k)-1] )
            
        kticks.append(i)
        i = i + k_interval

        end_tmp_k = kpoints.readline()
        end_tmp_k = end_tmp_k.rstrip('\n')
        end_tmp_k = end_tmp_k.rstrip()
        end_tmp_k = end_tmp_k.split(' ')
        
        kpoints.readline()

knames = str(knames)
kticks = str(kticks)
 # set knames and kticks manually
     
#knames= '[\"$\Gamma$\",\"T\",\'H_2|H_0\',\"L\",\"$\Gamma$\",\"S_0|S_2\",\"F\",\"$\Gamma$\"]' 
#kticks= '[0,  19,  39,   59,   79,   99,   119,   139]' 
#####################  From KPOINTS Get Knames And Ksticks  #######################



#####################  Define Function to Write Specific Plot (Orbital Mode)  ####################### 
def write_specific_orbital(ele, orb) :
    f.write("#-------------------------------------#")
    f.write('\n')
    f.write('#'+ele+'\n')
    f.write("atoms="+atom_dict[ele]+'\n')
    f.write("title=\'"+ele+"_"+orb+"\'\n")
    f.write("orbitals="+orbitals_name[orb]+'\n')
    f.write("pyprocar.bandsplot(dirname=\'.\',code=\'vasp\',"+"title=title,"+"fermi="+fermienergy+','+"elimit=[Emin,Emax],"+"cmap=cmap"+","+"mode=\'parametric\',"+"atoms=atoms"+","+"orbitals=orbitals"+","+"knames=knames"+","+"kticks=kticks"+","+"savefig=\'"+ele+"_"+orb+".png\')"+'\n')
      
#####################  Define Function to Write Specific Plot (Orbital Mode)  #######################       
      
#####################  Define Function to Write Specific Plot (Spin Mode)  #######################
def write_specific_spin() :
    f.write("#-------------------------------------#")
    f.write('\n')
    #f.write('#'+ele+'\n')
    #f.write("atoms="+atom_dict[ele]+'\n')
    f.write("title=\'"+"band"+"_"+"spin"+"\'\n")    
    #f.write("orbitals="+orbitals_name[orb]+'\n')
    #f.write("pyprocar.bandsplot(\'PROCAR-repaired\',outcar=\'OUTCAR-scf\',"+"elimit=[Emin,Emax],"+'cmap=cmap'+","+"mode=\'parametric\',"+"atoms=atoms"+","+"orbitals=orbitals"+","+"knames=knames"+","+"repair=False,"+"kticks=kticks"+","+"discontinuities=discontinuities,"+"savefig=\'"+ele+"_"+orb+".png\')"+'\n') 
    #f.write("pyprocar.bandsplot(dirname=\'.\',code=\'vasp\',"+"title=title,"+"fermi="+fermi_energy+','+"elimit=[Emin,Emax],"+"color=\'red\'"+","+"mode=\'plain\',"+"spins=[0],"+"knames=knames"+","+"kticks=kticks"+","+")"+'\n')
    f.write("pyprocar.bandsplot(dirname=\'.\',code=\'vasp\',"+"title=title,"+"fermi="+fermienergy+','+"elimit=[Emin,Emax],"+"color=\'blue\'"+","+"mode=\'plain\',"+"spins="+spins_mode+","+"knames=knames"+","+"kticks=kticks"+","+"savefig=\'"+'band'+"_"+'spin'+".png\')"+'\n')  
#####################  Define Function to Write Specific Plot (Spin Mode)  #######################        
       

with open('fatband.py' , mode='w') as f :
   
    f.write('#we need to use pyprocar, make sure that you have this module before running this file'+'\n')
    f.write('#we need to use POSCAR and KPOINTS to generate this file, you can also set it by yourself'+'\n')
    f.write('#we have POSCAR under the standard of Materials Project, it could be wrong if not under that standard'+'\n')
    f.write('#we have KPOINTS under the standard generated by vaspkit, it could be wrong if not under that standard'+'\n')
    f.write('#the basic frame was written by my senior \phq/(thanks), hope this file can help you ^_^'+'\n')     
    f.write("#-------------------------------------#")
    f.write('\n')
    f.write('import pyprocar '+'\n')
    
    f.write("#-------------------------------------#")
    f.write('\n')
    
    f.write('Emin='+Emin+'\n')
    f.write('Emax='+Emax+'\n')
    f.write('cmap='+cmap+'\n')
    f.write('knames='+knames+'\n')
    f.write('kticks='+kticks+'\n')
    f.write('discontinuities=[]'+'\n')
    
    f.write("#-------------------------------------#")
    f.write('\n')
    
    f.write("pyprocar.repair(\'PROCAR\',\'PROCAR-band\')"+'\n')
    
    if modeofplot == "o" :        
        for i in range(len(targetelements)) :
            write_specific_orbital(targetelements[i], targetorbitals[i])
    else :
        write_specific_spin()
    



