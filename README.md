# IVO-plugin
A plugin of psi4 for improved virtual orbital (STEX)

Usage: compile by typing

    psi4 --plugin-compile
    make
    
Include ivo directory in $PATHONPATH or state sys.path in psi4 input file.

Then in input file, add

    import ivo
    
    set ivo {
        print 1 # 1 for regular printing, 2 for debug
        hole  0 # 0 for the first MO. Using 1, 2 for other MOs.
    }
    
    energy('ivo')
    
Currently, it only works for close shell molecules.
