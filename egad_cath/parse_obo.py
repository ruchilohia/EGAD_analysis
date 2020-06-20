import traceback
import sys
import os
import numpy as np
from scipy import sparse
import pickle
GOASPECTMAP= { 'biological_process' : 'bp',
               'cellular_component' : 'cc',
               'molecular_function' : 'mf',
               'external'           : 'ex'
            }

UPASPECTMAP = { 'C': 'cc',
                'F': 'mf',
                'P': 'bp'
              }
def parse_obo():
    """
creates dict of dicts. key is goterm, contents is dict of 
       
       goterm  ""
       goname  ""
       goasp   ""
       godef   ""
       goisa   [] of goterms
       gohasa  [] of goterms
    
    
    """
    obofile = "/home/lohia/egad_trial/go.obo"
    filehandle = open(obofile)
    godict = {}
    altids = {}
    current = None
    try:
        for line in filehandle:
            if line.startswith("[Typedef]"):
                godict[current['goterm']]= current
                break
            elif line.startswith("[Term]"):     
                if current is not None:
                    godict[current['goterm']]= current
                # create new item...
                current = {}
                current['is_a'] = []
                current['part_of'] = []
                
            elif line.startswith("id: "):
                current['goterm'] = line[4:].strip()
                
            elif line.startswith("name: "):
                current['goname'] = line[6:].strip()
            
            elif line.startswith("namespace: "):
                asp = line[11:].strip()
                
                current['goasp'] = GOASPECTMAP[asp]

            # must create a separate mapping from alt_ids that maps to
            # primary, so it can be added to gomatrix properly
            elif line.startswith("alt_id: "):
                #current['alt_id'] = line[8:18].strip()
                #current['goasp'] = GOASPECTMAP[asp]            
                altid = line[8:18].strip()
                altids[altid] = current['goterm'] 
            
            #elif line.startswith("def: "):
            #    current['godef'] = line[5:].strip()

            #elif line.startswith("synonym: "):
            #    current.synonym.append(line[9:].strip())

            elif line.startswith("is_a: "):
                current['is_a'].append(line[6:16].strip())
            
            elif line.startswith("relationship"):
                if "part_of" in line:
                    current['part_of'].append(line[22:32])
                                 
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    ##print(altids)   
    return (godict, altids)

def converge_sparse(matrix):
    #logging.debug(f"starting matrix: \n{matrix_info(matrix)}")
    ##logging.debug(f"{#print_square(matrix.todense(), GOTERMLIST)}")
    oldval = 0
    #logging.debug("Summing inbound matrix...")
    newval = matrix.sum()
    #print(newval)
    #logging.debug("Beginning convergence loop.")
    icount = 0
    while oldval != newval:
        ##logging.debug(f"Inbound matrix:\n{matrix_info(matrix)}")
        ##logging.debug(f"oldval={oldval} newval={newval}")
        oldval = newval
        if not isinstance(matrix,  sparse.lil.lil_matrix): 
            #logging.debug(f"{type(matrix)} is not scipy.sparse.lil.lil_matrix, converting... ")
            matrix = sparse.lil_matrix(matrix, dtype=np.bool)
        else:
            pass
            ##logging.debug("matrix already lil_matrix...")
        ##logging.debug("Multiplying...")
        mult = matrix @ matrix
        ##logging.debug("Adding back original...")
        matrix = mult + matrix
        ##logging.debug("Getting new sum...")
        newval = matrix.sum()
        ##logging.debug(f"New matrix {icount}:\n{matrix.todense()}")
        ##logging.debug(f"{#print_square(matrix.todense(), GOTERMLIST)}")
        icount += 1
    #logging.debug(f"Convergence complete. {icount} iterations. matrix:\n{matrix_info(matrix)}")
    return matrix


def build_ontology(usecache):
    """
    obofile=~/data/go/go.obo
    cachedir = ~/play/cafa4      
    
    from parse_obo():
    { 'GO:2001315': 
         {'is_a': ['GO:0009226', 'GO:0046349', 'GO:2001313'], 
         'part_of': [], 
         'goterm': 'GO:2001315', 
         'goname': 'UDP-4-deoxy-4-formamido-beta-L-arabinopyranose biosynthetic process', 
         'goasp': 'bp', 
       ...  
    }
 
    result:  Numpy boolean matrix of all relations in ontology, ordered by sorted goterm name.  
       
    """
    #global GOMATRIX
    #global GOTERMIDX
    #global ALTIDDICT
    #global GOTERMLIST
    # def __init__(self, gomatrix, gotermidx, altidx):
    
    #logging.debug(f"usecache={usecache}")
    cachedir = "/home/lohia/egad_trial"
    ontologyfile = f"{cachedir}/ontology.npy"
    termidxfile = f"{cachedir}/gotermindex.pickle"
    altiddictfile = f"{cachedir}/altiddict.pickle"
    include_partof = "FALSE"
    
    gomatrix = None
    
    if os.path.exists(ontologyfile) and usecache:
        #logging.debug("Cache hit. Using existing matrix...")
        gomatrix = np.load(ontologyfile)
        #logging.debug(f"loaded matrix: {matrix_info(gomatrix)} from {ontologyfile}")
        
        f = open(termidxfile, 'rb')
        gotermidx = pickle.load(f)
        f.close()

        f = open(altiddictfile, 'rb')
        altiddict = pickle.load(f)
        f.close()
                
        #logging.debug(f"goterm index, e.g. : \n{list(gotermidx)[0]} :  {gotermidx[list(gotermidx)[0]]} ")
    
    else:
        (godict, altiddict) = parse_obo()
        
        # get keys from dict
        gotermlist = list(godict)
        ##print(gotermlist)
        #logging.debug(f"parsed obo with {len(gotermlist)} entries. ")
        #logging.debug(f"example entry:\n{gotermlist[0]}")
        #logging.debug("sorting goterms")
        gotermlist.sort()

        #logging.debug(f"sorted: e.g. {gotermlist[0:5]} ")
        #logging.debug("creating goterm index dict.")
        #
        i = 0
        gotermidx = {}
        for gt in gotermlist:
            gotermidx[gt] = i
            i = i + 1

              
        #logging.debug(f"creating zero matrix of dimension {len(gotermlist)}")
        shape = (len(gotermlist), len(gotermlist))
        gomatrix = np.zeros( shape, dtype=bool )
        
        #logging.debug(f"filling in parent matrix for all goterms...")
        for gt in godict.keys():
            for parent in godict[gt]['is_a']:
                    gomatrix[gotermidx[gt]][gotermidx[parent]] = True

        if include_partof:
            #logging.debug("Including part_of relationships as is_a")
            for gt in godict.keys():
                for parent in godict[gt]['part_of']:
                        gomatrix[gotermidx[gt]][gotermidx[parent]] = True
        
        ##logging.debug(f"initial matrix:\n{#print_square(gomatrix, GOTERMLIST)}")
        #logging.debug("Calculating sparsity...")
        #logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")
        #logging.debug("converting to sparse matrix.")
        #print(gomatrix)
        gomatrix = sparse.lil_matrix(gomatrix, dtype=bool)
        
        #logging.debug(f"converging matrix: {matrix_info(gomatrix)}")
        gomatrix = converge_sparse(gomatrix)
        #logging.info(f"got converged matrix:\n{matrix_info(gomatrix)} ")
        #logging.debug(f"converged matrix sum={gomatrix.sum()}")
        ##logging.debug("Calculating sparsity...")
        #sparsity = 1.0 - np.count_nonzero(gomatrix) / gomatrix.size
        ##logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")        
        gomatrix = gomatrix.todense()

        gomatrix = np.asarray(gomatrix)
        #print(gomatrix[0])
            
        #logging.debug(f"Caching all values/indexes...")
        #logging.debug(f"Saving matrix: {matrix_info(gomatrix)} to {ontologyfile}")
        np.save(ontologyfile, gomatrix)

        #logging.debug(f"Saving gotermidx {len(gotermidx)} items to {termidxfile}")
        f = open(termidxfile, 'wb')   
        pickle.dump(gotermidx, f )
        f.close() 

        f = open(altiddictfile, 'wb')
        pickle.dump(altiddict, f)
        f.close()
        
        #logging.debug("Done constructing input for Ontology().")
build_ontology("")
  print 
