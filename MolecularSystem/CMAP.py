from Geometry    import  distance_ab
from Molecule    import Molecule    
import numpy as np


class CMAP:
    """ Class doc """
    
    def __init__ (self, pdb  = False, cutoff = 6.0, log = False):
        """ Class initialiser """
        
        self.cutoff = 0
        self.cmap   = None
        self._type  = 'simple'
        self.number_of_contacts = 0
        
        if pdb:
            self.import_contact_map_from_pdb (  pdb    = pdb    , 
                                                cutoff = cutoff ,
                                                log    = log)
        
    
    def import_contact_map_from_pdb (self           , 
                                     pdb    = None  , 
                                     cutoff = 6.0   ,
                                     log    = False):
        
        """ Gera o mapa de contato a patir de um PDB (estrutura de referencia) """
        
       
        system = Molecule() 
        system.load_PDB_to_system (filename = pdb   )   
        cmap = 0*np.random.rand(len(system.residues),len(system.residues))
        
        self.number_of_contacts =  0
        
        for index_i in range(0, len(system.residues)):
            for index_j in range(index_i+2, len(system.residues)):
                
                name_i = system.residues[index_i].name
                
                for atom in system.residues[index_i].atoms:
                    if atom.name == 'CA':
                        atom_i    = atom  
                name_j = system.residues[index_j].name
                
                for atom in system.residues[index_j].atoms:
                    if atom.name == 'CA':
                        atom_j = atom  
                R_ab = distance_ab (atom_i, atom_j)
                
                
                if R_ab <= cutoff:
                    cmap[index_i][index_j] = 1
                    self.number_of_contacts += 1
                    if log:
                        print index_i, name_i, index_j, name_j, R_ab, 'contact'
                else:
                    pass
        
        if log:
            print '\ncmap matrix:'
            print cmap
            
            print '\n---------------------------------------------'
            print 'total number of residues = ', len(system.residues)
            print 'Cutoff size              = ', cutoff
            print 'number of contacts       = ', self.number_of_contacts
            print '---------------------------------------------\n'
        
        self.cmap   = cmap
        self.cutoff = cutoff
