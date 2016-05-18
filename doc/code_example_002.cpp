"""
    source/
          |--> mat/  (GPU matrix operation library).
          |       |
          |       |--> include/ (include path for the compiler).
          |       |--> main.cpp (entry testing file). 
          |       |--> src/ (source code). 
          |       |--> ... 
          |
          |--> lift/  (Lattice Interface For Tight binding).
                   |
                   |--> bin/ (place for the python code).
                   |       |
                   |       |--> lift.py (python script
                   |       |              for processing the lattice file).
                   |       |--> vlift.py (python script
                   |       |              for visualize input/output file).
                   |       |--> wxBand.py (python script
                   |       |              for visualize calculated bandstructure).
                   |       |--> mhg (simple script for compile the GPU based code).
                   |       |--> m16 (compile environment settings).
                   |
                   |--> include/ (include path for the compiler).
                   |--> src/ (TBM^cube source code).
                   |--> examples/ (several example code of the TBM^cube). 
                   |            |
                   |            |--> main_cubic.cpp (the first example). 
                   |            |--> ...
                   |
                   |--> doc/ (the place for this document). 
                   |--> ... 
"""
