# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

import labeling


# TODO: This needs to be in a class, and eventually refer to an exernal data source.
def getFragDictBasic():
    """
    Function to create basic dictionary of fragments
    As explained above the inputs for the GCMS fragment go as:
        
    Ala159 =  GCMSfragment('alaL_2_3'              ,'Alanine,159'      ,'Ala159',2,'C6H20NSi'    ,0,317,'C15H35NO2Si2')

    alaL_2_3      = EMU name
    Alanine,159   = fragment name
    Ala159        = abbreviation
    2             = number of carbons in fragment
    C6H20NSi      = formula for fragment without carbon backbone
    0             = antoniewicz approved?
    317           = weight for derivatized fragment
    C15H35NO2Si2  = formula for derivatized fragment
    """

    # Defining basic GCMS Dictionary
    #                          EMU name                 fragment name      abbreviation, number of carbons, formula minus backbone, approved by antoniewicz, weight of derivatized fragment, formula of derivatized fragment
    Ala159 =  labeling.GCMSfragment('alaL_2_3'              ,'Alanine,159'      ,'Ala159',2,'C6H20NSi'    ,0,317,'C15H35NO2Si2')
    Ala85  =  labeling.GCMSfragment('alaL_2_3'              ,'Alanine,85'       ,'Ala85' ,2,'C8H26NOSi2'  ,1,317,'C15H35NO2Si2')
    Ala57  =  labeling.GCMSfragment('alaL_1_2_3'            ,'Alanine,57'       ,'Ala57' ,3,'C8H26NO2Si2' ,1,317,'C15H35NO2Si2')
    Gly159 =  labeling.GCMSfragment('gly_2'                 ,'Glycine,159'      ,'Gly159',1,'C6H18NSi'    ,0,303,'C14H33NO2Si2')
    Gly85  =  labeling.GCMSfragment('gly_2'                 ,'Glycine,85'       ,'Gly85' ,1,'C8H24NOSi2'  ,1,303,'C14H33NO2Si2')
    Gly57  =  labeling.GCMSfragment('gly_1_2'               ,'Glycine,57'       ,'Gly57' ,2,'C8H24NO2Si2' ,1,303,'C14H33NO2Si2')
    Val85  =  labeling.GCMSfragment('valL_2_3_4_5'          ,'Valine,85'        ,'Val85' ,4,'C8H30NOSi2'  ,1,345,'C17H39NO2Si2')
    Val159 =  labeling.GCMSfragment('valL_2_3_4_5'          ,'Valine,159'       ,'Val159',4,'C6H24NSi'    ,0,345,'C17H39NO2Si2')
    Val57  =  labeling.GCMSfragment('valL_1_2_3_4_5'        ,'Valine,57'        ,'Val57' ,5,'C8H30NO2Si2' ,1,345,'C17H39NO2Si2')
    Leu159 =  labeling.GCMSfragment('leuL_2_3_4_5_6'        ,'Leucine,159'      ,'Leu159',5,'C6H26NSi'    ,0,359,'C18H41NO2Si2')
    Leu85  =  labeling.GCMSfragment('leuL_2_3_4_5_6'        ,'Leucine,85'       ,'Leu85' ,5,'C8H32NOSi2'  ,1,359,'C18H41NO2Si2')
    Leu57  =  labeling.GCMSfragment('leuL_1_2_3_4_5_6'      ,'Leucine,57'       ,'Leu57' ,6,'C8H32NO2Si2' ,0,359,'C18H41NO2Si2')
    Ile159 =  labeling.GCMSfragment('ileL_2_3_4_5_6'        ,'Isoleucine,159'   ,'Ile159',5,'C6H26NSi'    ,0,359,'C18H41NO2Si2')
    Ile85  =  labeling.GCMSfragment('ileL_2_3_4_5_6'        ,'Isoleucine,85'    ,'Ile85' ,5,'C8H32NOSi2'  ,1,359,'C18H41NO2Si2')
    Ile57  =  labeling.GCMSfragment('ileL_1_2_3_4_5_6'      ,'Isoleucine,57'    ,'Ile57' ,6,'C8H32NO2Si2' ,0,359,'C18H41NO2Si2')
    Ser159 =  labeling.GCMSfragment('serL_2_3'              ,'Serine,159'       ,'Ser159',2,'C12H34NOSi2' ,1,447,'C21H49NO3Si3')
    Ser85  =  labeling.GCMSfragment('serL_2_3'              ,'Serine,85'        ,'Ser85' ,2,'C14H40NO2Si3',1,447,'C21H49NO3Si3')
    Ser57  =  labeling.GCMSfragment('serL_1_2_3'            ,'Serine,57'        ,'Ser57' ,3,'C14H40NO3Si3',1,447,'C21H49NO3Si3')
    Thr85  =  labeling.GCMSfragment('thrL_2_3_4'            ,'Threonine,85'     ,'Thr85' ,3,'C14H42NO2Si3',1,461,'C22H51NO3Si3')
    Thr57  =  labeling.GCMSfragment('thrL_1_2_3_4'          ,'Threonine,57'     ,'Thr57' ,4,'C14H42NO3Si3',1,461,'C22H51NO3Si3')
    Met159 =  labeling.GCMSfragment('metL_2_3_4_5'          ,'Methionine,159'   ,'Met159',4,'C6H24NSSi'   ,1,377,'C17H39NO2SSi2')
    Met85  =  labeling.GCMSfragment('metL_2_3_4_5'          ,'Methionine,85'    ,'Met85' ,4,'C8H30NOSSi2' ,1,377,'C17H39NO2SSi2')
    Met57  =  labeling.GCMSfragment('metL_1_2_3_4_5'        ,'Methionine,57'    ,'Met57' ,5,'C8H30NO2SSi2',1,377,'C17H39NO2SSi2')
    Phe159 =  labeling.GCMSfragment('pheL_2_3_4_5_6_7_8_9'  ,'Phenylalanine,159','Phe159',8,'C6H24NSi'    ,1,393,'C21H39NO2Si2')
    Phe85  =  labeling.GCMSfragment('pheL_2_3_4_5_6_7_8_9'  ,'Phenylalanine,85' ,'Phe85' ,8,'C8H30NOSi2'  ,1,393,'C21H39NO2Si2')
    Phe57  =  labeling.GCMSfragment('pheL_1_2_3_4_5_6_7_8_9','Phenylalanine,57' ,'Phe57' ,9,'C8H30NO2Si2' ,1,393,'C21H39NO2Si2')
    Phe91  =  labeling.GCMSfragment('pheL_1_2'              ,'Phenylalanine,91' ,'Phe91' ,2,'C12H33NO2Si2',1,393,'C21H39NO2Si2')
    Asp159 =  labeling.GCMSfragment('aspL_2_3_4'            ,'Aspartate,159'    ,'Asp159',3,'C12H34NO2Si2',0,475,'C22H49NO4Si3')
    Asp85  =  labeling.GCMSfragment('aspL_2_3_4'            ,'Aspartate,85'     ,'Asp85' ,3,'C14H40NO3Si3',1,475,'C22H49NO4Si3')
    Asp57  =  labeling.GCMSfragment('aspL_1_2_3_4'          ,'Aspartate,57'     ,'Asp57' ,4,'C14H40NO4Si3',1,475,'C22H49NO4Si3')
    Asp173 =  labeling.GCMSfragment('aspL_1_2'              ,'Aspartate,173'    ,'Asp173',2,'C12H32NO2Si2',1,475,'C22H49NO4Si3')
    Asp99  =  labeling.GCMSfragment('aspL_1_2'              ,'Aspartate,99'     ,'Asp99' ,2,'C14H38NO3Si3',1,475,'C22H49NO4Si3')
    Glu159 =  labeling.GCMSfragment('gluL_2_3_4_5'          ,'Glutamate,159'    ,'Glu159',4,'C12H36NO2Si2',1,489,'C23H51NO4Si3')
    Glu85  =  labeling.GCMSfragment('gluL_2_3_4_5'          ,'Glutamate,85'     ,'Glu85' ,4,'C14H42NO3Si3',1,489,'C23H51NO4Si3')
    Glu57  =  labeling.GCMSfragment('gluL_1_2_3_4_5'        ,'Glutamate,57'     ,'Glu57' ,5,'C14H42NO4Si3',1,489,'C23H51NO4Si3')
    Tyr57  =  labeling.GCMSfragment('tyrL_1_2_3_4_5_6_7_8_9','Tyrosine,57'      ,'Tyr57' ,9,'C14H44NO3Si3',0,523,'C27H53NO3Si3')
    Tyr159 =  labeling.GCMSfragment('tyrL_2_3_4_5_6_7_8_9'  ,'Tyrosine,159'     ,'Tyr159',8,'C12H38NOSi2' ,0,523,'C27H53NO3Si3')
    Tyr221 =  labeling.GCMSfragment('tyrL_1_2'              ,'Tyrosine,221'     ,'Tyr221',2,'C12H32NO2Si2',1,523,'C27H53NO3Si3')
    Pro159 =  labeling.GCMSfragment('proL_2_3_4_5'          ,'Proline,159'      ,'Pro159',4,'C6H22NSi'    ,0,343,'C17H37NO2Si2')

    Ala   =  labeling.LCMSfragment('alaL_1_2_3'            ,'Ala,0'             ,'Ala' ,3,'H5NO')
    Gly   =  labeling.LCMSfragment('gly_1_2'               ,'Gly,0'             ,'Gly' ,2,'H3NO')
    Val   =  labeling.LCMSfragment('valL_1_2_3_4_5'        ,'Val,0'             ,'Val' ,5,'H9NO')
    Leu   =  labeling.LCMSfragment('leuL_1_2_3_4_5_6'      ,'Leu,0'             ,'Leu' ,6,'H11NO')
    Ile   =  labeling.LCMSfragment('ileL_1_2_3_4_5_6'      ,'Ile,0'             ,'Ile' ,6,'H11NO')
    Ser   =  labeling.LCMSfragment('serL_1_2_3'            ,'Ser,0'             ,'Ser' ,3,'H5NO2')
    Thr   =  labeling.LCMSfragment('thrL_1_2_3_4'          ,'Thr,0'             ,'Thr' ,4,'H7NO2')
    Met   =  labeling.LCMSfragment('metL_1_2_3_4_5'        ,'Met,0'             ,'Met' ,5,'H9NOS')
    Phe   =  labeling.LCMSfragment('pheL_1_2_3_4_5_6_7_8_9','Phe,0'             ,'Phe' ,9,'H9NO')
    Asp   =  labeling.LCMSfragment('aspL_1_2_3_4'          ,'Asp,0'             ,'Asp' ,4,'H5NO3')
    Glu   =  labeling.LCMSfragment('gluL_1_2_3_4_5'        ,'Glu,0'             ,'Glu' ,5,'H7NO3')
    Tyr   =  labeling.LCMSfragment('tyrL_1_2_3_4_5_6_7_8_9','Tyr,0'             ,'Tyr' ,9,'H9NO2')
    Pro   =  labeling.LCMSfragment('proL_1_2_3_4_5'        ,'Pro,0'             ,'Pro' ,5,'H7NO')
    Lys   =  labeling.LCMSfragment('lysL_1_2_3_4_5_6'      ,'Lys,0'             ,'Lys' ,6,'H12N2O')
    Arg   =  labeling.LCMSfragment('argL_1_2_3_4_5_6'      ,'Arg,0'             ,'Arg' ,6,'H12N4O')
    His   =  labeling.LCMSfragment('hisL_1_2_3_4_5_6'      ,'His,0'             ,'His' ,6,'H7N3O')
    Trp   =  labeling.LCMSfragment('trpL_1_2_3_4_5_6_7_8_9_10_11','Trp,0'      ,'Trp'  ,11,'H10N2O')
    Cys   =  labeling.LCMSfragment('cysL_1_2_3'             ,'Cys,0'            ,'Cys' ,3,'H5NOS')
    Asn   =  labeling.LCMSfragment('asnL_1_2_3_4'           ,'Asn,0'            ,'Asn' ,4,'H6N2O2')
    Gln   =  labeling.LCMSfragment('glnL_1_2_3_4_5'         ,'Gln,0'            ,'Gln' ,5,'H8N2O2')

    Suc   =  labeling.LCMSfragment('succ_1_2_3_4'          ,'succ,0'            ,'succ'  ,4,'H4O4')    
    Mal   =  labeling.LCMSfragment('mal-L_1_2_3_4'         ,'mal-L,0'           ,'mal-L' ,4,'H4O5')
    Cit   =  labeling.LCMSfragment('cit_1_2_3_4_5_6'       ,'cit,0'             ,'cit'   ,6,'H5O7')    
    Pyr   =  labeling.LCMSfragment('pyr_1_2_3'             ,'pyr,0'             ,'pyr'   ,3,'H3O3')  
    Mev   =  labeling.LCMSfragment('mev-R_1_2_3_4_5_6'     ,'mev-R,0'           ,'mev-R' ,6,'H11O4')  
    Malcoa=  labeling.LCMSfragment('malcoa_1_2_3'          ,'malcoa,0'          ,'malcoa',3,'HN7O3')  
    Akg   =  labeling.LCMSfragment('akg_1_2_3_4_5'         ,'akg,0'             ,'akg'   ,5,'H4O5')  
    Oaa   =  labeling.LCMSfragment('oaa_1_2_3_4'           ,'oaa,0'             ,'oaa'   ,4,'H2O5')  
    G6p   =  labeling.LCMSfragment('g6p_1_2_3_4_5_6'       ,'g6p,0'             ,'g6p'   ,6,'H11O9P')  
    Glx   =  labeling.LCMSfragment('glx_1_2'               ,'glx,0'             ,'glx'   ,6,'H1O3')  
    octa  =  labeling.LCMSfragment('octa_1_2_3_4_5_6_7_8'  ,'octa,0'            ,'octa'  ,8,'H15O2')  
    dca   =  labeling.LCMSfragment('dca_1_2_3_4_5_6_7_8_9_10','dca,0'           ,'dca'  ,10,'H20O2')  
    ddca  =  labeling.LCMSfragment('ddca_1_2_3_4_5_6_7_8_9_10_11_12','mev-R,0'  ,'ddca' ,12,'H23O2') 

  
    fdp    =  labeling.LCMSfragment('fdp_1_2_3_4_5_6'     ,'fdp,0'            ,'fdp'   ,6,'H10O2P' )
    dhap   =  labeling.LCMSfragment('dhap_1_2_3'          ,'dhap,0'           ,'dhap'  ,3,'H5O6P'  )
    M_3pg  =  labeling.LCMSfragment('3pg_1_2_3'           ,'3pg,0'            ,'3pg'   ,3,'H4O7P'  )  # We can't have a variable starting with a number
    pep    =  labeling.LCMSfragment('pep_1_2_3'           ,'pep,0'            ,'pep'   ,3,'H2O6P'  )
    pyr    =  labeling.LCMSfragment('pyr_1_2_3'           ,'pyr,0'            ,'pyr'   ,3,'H3O3'   )
    ru5p   =  labeling.LCMSfragment('ru5p-D_1_2_3_4_5'    ,'ru5p-D,0'         ,'ru5p-D',5,'H9O8P'  )
    r5p    =  labeling.LCMSfragment('r5p_1_2_3_4_5'       ,'r5p,0'            ,'r5p'   ,5,'H9O8P'  )
    s7p    =  labeling.LCMSfragment('s7p_1_2_3_4_5_6_7'   ,'s7p,0'            ,'s7p'   ,7,'H13O10P')
    mal    =  labeling.LCMSfragment('mal-L_1_2_3_4'       ,'mal-L,0'          ,'mal-L' ,4,'H4O5'   )

    allFragments =   [Ala159,Ala85,Ala57,Gly159,Gly85,Gly57,Val85,Val159,Val57,Leu159,Leu85,Leu57,Ile159,Ile85,Ile57,
                      Ser159,Ser85,Ser57,Thr85,Thr57,Met159,Met85,Met57,Phe159,Phe85,Phe57,Phe91,Asp159,Asp85,Asp57,
                      Asp173,Asp99,Glu159,Glu85,Glu57,Tyr57,Tyr159,Tyr221,Pro159]

    allFragments.extend([Ala,Gly,Val,Leu,Ile,Ser,Thr,Met,Phe,Asp,Glu,Tyr,Pro,Lys,Arg,His,Trp,Cys,Asn,Gln])

    allFragments.extend([Suc,Mal,Cit,Pyr,Mev,Malcoa,Akg,Oaa,G6p,Glx,octa,dca,ddca])

    allFragments.extend([fdp,dhap,M_3pg,pep,pyr,ru5p,r5p,s7p,mal])

    return allFragments

        
def metAliases():
    """ Correspondance between labeling name and aliases """
     
    # TODO(Hector): Eliminate this completely, or substitute by aliases
    MetDict   = {
                'Glucose'   :['glc-D[e]','glcDE','glc_D_e','glcDE_c'],
                'Acetate'   :['ac[e]','ac_e'],
                'Ethanol'   :['etoh[e]','etoh_e'],
                'AcetylCoa' :['accoa','accoa[e]','accoa_c'],
                'Aspartate' :['Asp','asp','asp_c'],
                'Pyruvate'  :['Pyr','pyr','pyr_c'],
                'XX'        :['XX','xx','XX_c','xx_c']
                     }   
                     
    return MetDict


def getCompartmentDict():
    "A reference for all compartments for SBML models, e.g. c = cytosolic, m = mitochondrial, etc"
    # Note: Compartment short name 'b' is not in BiGG.
    # Nevertheless, it appears in some SBML, and is usually paired with the 'boundaryCondition = true' flag when it does appear. 
    compartments = {
        'c': 'cytosol',
        'e': 'extracellular space',
        'f': 'flagellum',
        'g': 'golgi apparatus',
        'h': 'chloroplast',
        'i': 'intracellular',
        'l': 'lysosome',
        'm': 'mitochondria',
        'n': 'nucleus',
        'p': 'periplasm',
        'r': 'endoplasmic reticulum',
        's': 'eyespot',
        'u': 'thylakoid',
        'v': 'vacuole',
        'x': 'peroxisome/glyoxysome',
        'b': 'b'    # See note above about compartment 'b'
        }
    return compartments
