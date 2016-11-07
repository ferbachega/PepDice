import os

cwd     = os.getcwd()
folders = os.listdir('.')


for folder in folders:
    try:
        path = os.path.join(cwd, folder)
        os.chdir(os.path.join(cwd, folder))
        #os.system('rm *decoy*.crd')
        #os.system('rm *decoy*.seq')
        #os.system('rm *decoy*.top')
        #os.system('rm *decoy*A.pdb')

        os.system('rm *AMBER.pdb')
        
		#os.chdir('original')
        #os.system('rm *decoy*.pdb')
        #os.system('mv list.txt '+ path)
        #os.system('mv rst.dat '+ path)
    
    except:
        
        print 'The path is not a folder:' ,os.path.join(cwd, folder)
        
