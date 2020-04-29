Potential Issues and Solutions
==============================

1.  Occasionally one of the AN spike train files, which are stored in the directory `cnmodel/an_model/cache`, become locked. This can occur if the calling routines (e.g., simulation runs) are aborted (^C, ^Z) in the middle of a transaction accessing the cache file, or perhaps during when parallel processing is enabled and a routine fails or is aborted. In this case, a file with the extension ``".lock"`` exists, which prevents the an_model code from accessing the file. The ``".lock"`` file needs to be deleted from the cache directory. Because the cache directory contains an hierarchical arrangement of subdirectories, and can be populated with thousands of files after a few runs requiring many auditory nerve datasets, finding the lock file can be somewhat tedious. The following should help under Unix:
    
  *  First, print a list of the locked files::
      
          - find /path/to/cache -name '*.lock'
    
  * Where /path/to/cache may be something like `cnmodel/an_model/cache`. 
    There is most likely only one such file in the diretory.

  * Next, to delete the files::
  
      - find /path/to/cache -name '*.lock' -delete
       
  * Under Windows (and other OS's), you should be able do accomplish the same thing
    with the File Explorer/Finder, limiting the files by extension.
    
  * An alternative (for any OS) is to take advantage of Python's pathlib module. The glob search is 
    remarkably fast (on my system, it takes under a minute to search through more than 3.5 million
    cached AN spike trains)::
    
            >>python
            > from pathlib import Path
            > gl = Path('.').rglob('*.lock')
            > locks = list(gl) # (could do this in the next line)
            > # print(locks)  # see the lock files
            > for g in locks:  # now remove the lock files
            >    g.unlink()
            >
   