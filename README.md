# stellar-environments

to do    :: vhone.f90 - add Dr. B's cooling algorithm
            vhone.f90 - fix accretion issues


02/19/18 :: updated vhone.f90
             - shock finding algorithm only executes in stage 3
             - bugfix accretion algorithm for stampede
            updated init.f90
             - disabled noise algorithm, rand() conflict on stampede
            updated make/stampede
             - now includes read.f90 
            
   tests :: ran full sim on stampede (branched)
            ran full sim on kepler (merged)


02/05/18 :: updated vhone.f90 
             - added simple accretion algorithm
             - removed annoying status printout
            updated init.f90
             - removed annoying status printout



01/29/18 :: updated read.f90, vhone.f90, init.f90 to use icons
            updated icons - new read procedure
            discont icons2 - unused
            updated setup.py - now copies icons instead of icons2

   tests :: ran full sim with GM2=0
